#include "Linear-Algebra-C/linear-algebra.h"
#include "polylib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

// --- Global Variables ---
double time_int = 10.0;
int time_points = 6000;
double time_delta; 


double x_0 = 0.0, x_f = 9575.0;     // Dimensions for the Marmousi
double y_0 = 0.0, y_f = 4025.0;
double delta_x, delta_y;
double Lx, Ly;
int n_x = 20, n_y = 20;            // Number of intervals per axis
int n_elements;
int n_nodes = 6;      // Nodes per element per direction
int total_points;     // Will be calculated in main, total nodes in Omega
matrix* xy_points = NULL;           // In row n, coordinates of point n in global ID
matrix* center_element = NULL;     // Row n contains the center of element n
matrix** local_K = NULL;
vector* elasticity = NULL;

int** connectivity = NULL;

// Initialization of the GLL points and Derivatives
vector* gll_points = NULL;
vector* gll_weights = NULL;
matrix* D = NULL;
matrix* Dt = NULL;
vector* Ku = NULL;

// --- Function Prototypes ---
double density(double x, double y);
void calc_elasticity(matrix *velocity, vector *elasticity, int rows, int cols);
double f(double x, double y, double t);
void F_e(double xi, double nu, int e, double *output);
void init(int **connectivity, matrix* xy);          // Initializes the connectivity matrix and the points of interpolation
double T(int m, int n, int i, int j, int l, int k);
void initial_conditions(matrix* u, vector* u_vel);
double K(int i, int j, int m, int n, int e);
void comp_Ku(vector* Ku, matrix* u);                // Computes the product Ku
void saveStepToFile(FILE *fp, matrix* u, int total_points);
double analytical_sol(double x, double y, double t);
matrix* load_from_xyz(const char* filename, int* out_rows, int* out_cols, double* out_dx, double* out_dy);
double calc_velocity(double x, double y, int rows, int cols, matrix *velocity_grid);

int main()
{
    // Initialize distances and intervals
    Lx = x_f - x_0;
    Ly = y_f- y_0;
    delta_x = Lx/n_x;
    delta_y = Ly/n_y;
    n_elements = n_x * n_y;
    time_delta = time_int / (time_points - 1);

    // File to save results
    FILE *results_fp = fopen("simulation_results_2d.csv", "w");

    // Calculate the center of each element
    center_element = zeroMatrix(n_elements, 2);
    int e = 0;
    for(int j = 0; j<n_y; j++){
        for(int i = 0; i<n_x; i++){
            double center_x = x_0 + delta_x/2 + i*delta_x;
            double center_y = y_f - delta_y/2 - j*delta_y;
            setMatrixElement(center_element, e, 0, center_x);
            setMatrixElement(center_element, e, 1, center_y);
            e++;
        }
    }
    // printMatrix(center_element, true);

    // Allocate space for the connectivity matrix. It was not done with the linear algebra library as it stores doubles and not int.
    int rows = n_nodes * n_nodes;
    int cols = n_elements;
    connectivity = (int**)malloc(rows * sizeof(int*));
    if (connectivity == NULL) return 1;
    for (int i = 0; i < rows; i++) {
        connectivity[i] = (int*)malloc(cols * sizeof(int));
        if (connectivity[i] == NULL) return 1;
    }

    // Initialization of nodes weights and points
    total_points = (n_x*(n_nodes - 1) + 1)*(n_y*(n_nodes - 1) + 1);     
    gll_points = nullVector(n_nodes);
    gll_weights = nullVector(n_nodes);
    xy_points = zeroMatrix(total_points, 2);

    zwgll(gll_points->data, gll_weights->data, n_nodes);
    init(connectivity, xy_points);

    // --- Save X Coordinates ---
    int nodes_printed = 0;
    for (int i = 0; i < total_points; i++) {
        fprintf(results_fp, "%e%s", getMatrixElement(xy_points, i, 0), (i == total_points - 1) ? "" : ",");
        nodes_printed++;
    }
    fprintf(results_fp, "\n");

    // --- Save Y Coordinates ---
    for (int i = 0; i < total_points; i++) {
        fprintf(results_fp, "%e%s", getMatrixElement(xy_points, i, 1), (i == total_points - 1) ? "" : ",");
    }
    fprintf(results_fp, "\n");
    printf("Saving %d points to CSV\n", nodes_printed);


    // Calculation of the elasticity using the marmousi file
    int rows_marmousi, cols_marmousi;
    double dx, dy;

    matrix* marmousi = load_from_xyz("marmousi.xyz", &rows_marmousi, &cols_marmousi, &dx, &dy);
    elasticity = zeroVector(total_points);
    calc_elasticity(marmousi, elasticity, rows_marmousi, cols_marmousi);
    printf("Elasticity at random point %lf\n", getVectorElement(elasticity, 500));

    // Initialization of M and f
    vector* M = zeroVector(total_points);
    vector* vector_f = zeroVector(total_points);
    Ku = zeroVector(total_points);

    // Compute the derivatives of the Legendre Polynomials
    D = zeroMatrix(n_nodes, n_nodes);
    Dt = zeroMatrix(n_nodes, n_nodes);
    Dglj(D->data, Dt->data, gll_points->data, n_nodes, 0.0, 0.0);


    printf("DEBUG: Starting M assembly\n");
    // Computation of M
    for (int e=0; e<n_elements; e++){
        double jacobian = delta_x*delta_y/4;      // For equidistant nodes
        double density_ij;
        for(int i=0; i<n_nodes; i++){
            for(int j=0; j<n_nodes; j++){
                int node_index = connectivity[i*n_nodes+j][e];
                double w_i = getVectorElement(gll_weights, i);
                double w_j = getVectorElement(gll_weights, j);
                double x_coord = getMatrixElement(xy_points, node_index, 0);
                double y_coord = getMatrixElement(xy_points, node_index, 1);
                density_ij = density(x_coord, y_coord);
                
                // Calculation accounts for the possible overlapping

                double M_ij = getVectorElement(M, node_index);
                double calc = density_ij*w_i*w_j*jacobian;       
                setVectorElement(M, node_index, M_ij+calc);
            }
        }
    }
    // Computation of the inverse of M
    vector* M_inv = zeroVector(total_points);
    for(int i = 0; i<total_points; i++){
        double M_inv_i = 1./getVectorElement(M, i);
        setVectorElement(M_inv, i, M_inv_i);
    }


    // we store the past, present and future (how poetic) and initial velocity 
    matrix* u[3];
    u[0] = zeroMatrix(1, total_points);
    u[1] = zeroMatrix(1, total_points);
    u[2] = zeroMatrix(1, total_points);
    vector* u_vel = zeroVector(total_points);
    initial_conditions(u[0], u_vel);

    // Computation of local K's
    local_K = (matrix**)malloc(n_elements * sizeof(matrix*));
    printf("DEBUG: Starting K assembly\n");
    for (int e = 0; e < n_elements; e++) {
        // Each element gets its own local K
        local_K[e] = zeroMatrix(n_nodes*n_nodes, n_nodes*n_nodes); 

        for(int m = 0; m<n_nodes; m++){
            for(int n = 0; n<n_nodes; n++){

                for (int i = 0; i < n_nodes; i++) {
                    for (int j = 0; j < n_nodes; j++) {
                    
                        double val = K(i, j, m, n, e);
                        setMatrixElement(local_K[e], i*n_nodes+j, m*n_nodes+n, val);
                    }
                }
            }
        }
    }

    vector* Ku = zeroVector(total_points);
    // Time loop 
    printf("Starting time steps\n");
    for(int t_i=0; t_i<time_points-1; t_i++){
        double t = t_i*time_delta;
        // printf("DEBUG: Entering Time Loop t_i=%d\n", t_i);
        fillVector(&vector_f, 0.0);
        fillVector(&Ku, 0.0);

        double jacobian = delta_x*delta_y/4;      // For equidistant nodes

        // Computation of vector f
        for(int e = 0; e<n_elements; e++){
            
            for(int i=0; i<n_nodes; i++){
                for(int j=0; j<n_nodes; j++){
                        double xi_i = getVectorElement(gll_points, i);
                        double w_i = getVectorElement(gll_weights, i);
                        double nu_j = getVectorElement(gll_points, j);
                        double w_j = getVectorElement(gll_weights, j);

                        int global_id = connectivity[i*n_nodes+j][e];
                        double x = getMatrixElement(xy_points, global_id, 0);
                        double y = getMatrixElement(xy_points, global_id, 1);
                        
                        double f_eval_ij = f(x, y, t);
                        // Calculation accounts for the possible overlapping

                        double vector_f_ij = getVectorElement(vector_f, global_id);
                        double calc = f_eval_ij*jacobian*w_i*w_j;       
                        setVectorElement(vector_f, global_id, vector_f_ij+calc);
                }
            }
        }
        

        //  Case u_1
        // We first compute u_1 (first instant after 0) as it is done with a different approximation to \ddot u
        if(t_i == 0){
            // Computation of Ku
            comp_Ku(Ku, u[0]);
            for(int m=0; m<total_points; m++){
                double x = getMatrixElement(xy_points, m, 0);
                double y = getMatrixElement(xy_points, m, 1);

                // Do not calculate if the node is in the edges (Neumann conditions)
                double eps = 1e-9; // Small tolerance
                if (x <= x_0 + eps || x >= x_f - eps || y <= y_0 + eps || y >= y_f - eps) {
                    setMatrixElement(u[1], 0, m, 0.0); 
                }
                else{
                    double Ku_m = getVectorElement(Ku, m);
                    double M_inv_m = getVectorElement(M_inv, m);
                    double vector_f_m = getVectorElement(vector_f, m);
                    double a_m = (vector_f_m - Ku_m) * M_inv_m;     

                    double sum = getMatrixElement(u[0], 0, m) + time_delta*getVectorElement(u_vel, m) + time_delta*time_delta/2*a_m;
                    setMatrixElement(u[1], 0, m, sum);
                }
            }
            saveStepToFile(results_fp, u[1], total_points);
        }

        // Other cases
        else{
            // Computation of Ku
            comp_Ku(Ku, u[1]);
            for (int d = 0; d < total_points; d++) {
                double x = getMatrixElement(xy_points, d, 0);
                double y = getMatrixElement(xy_points, d, 1);

                double eps = 1e-9; // Small tolerance
                if (x <= x_0 + eps || x >= x_f - eps || y <= y_0 + eps || y >= y_f - eps) {
                    setMatrixElement(u[2], 0, d, 0.0); 
                }
                else{
                    double Ku_d = getVectorElement(Ku, d);

                    // Now use Ku_d in the central difference formula
                    double M_inv_d = getVectorElement(M_inv, d);
                    double F_d = getVectorElement(vector_f, d);
                    double u_curr = getMatrixElement(u[1], 0, d);
                    double u_old  = getMatrixElement(u[0], 0, d);

                    double accel = (F_d - Ku_d) * M_inv_d;
                    double u_next = (time_delta * time_delta * accel) + (2.0 * u_curr) - u_old;

                    setMatrixElement(u[2], 0, d, u_next);
                }
            }
            // Save results and swap the generations
            if(t_i % 50 == 0) saveStepToFile(results_fp, u[2], total_points);
            matrix* temp = u[0];
            u[0] = u[1];
            u[1] = u[2];
            u[2] = temp;
        }
        if(t_i % 1000 == 0) printf("Time step %d completed\n", t_i);
    }


    
    printf("Simulation done!\n");

    // Free matrices, vectors, etc.
    for (int i = 0; i < rows; i++) {
        free(connectivity[i]);
    }
    for (int e = 0; e < n_elements; e++) {
        deleteMatrix(local_K[e]); 
    }  
    for(int i = 0; i<3; i++){
        deleteMatrix(u[i]);
    }
    free(connectivity);
    deleteMatrix(center_element);
    deleteMatrix(xy_points);
    deleteVector(gll_points);
    deleteVector(gll_weights);
    deleteVector(M);
    deleteMatrix(D);
    deleteMatrix(Dt);
    deleteVector(vector_f);
    deleteVector(Ku);
    deleteVector(M_inv);
    free(local_K);
    return 0;
}

double density(double x, double y) {
    return 1.0;
}

// double density(double x, double y) {
//     double cx = 5.0; // Center X
//     double cy = 3.0; // Center Y
//     double sigma = 2.0;
//     // Creates a bell-curve density centered at (5,3)
//     return 1.0 + 5.0 * exp(-((x - cx) * (x - cx) + (y - cy) * (y - cy)) / (2 * sigma * sigma));
// }

void calc_elasticity(matrix *velocity, vector *elasticity, int rows, int cols) {

    for(int i = 0; i<total_points; i++) {
        double x_coord = getMatrixElement(xy_points, i, 0);
        double y_coord = getMatrixElement(xy_points, i, 1);
        double vel_i = calc_velocity(x_coord, y_coord, rows, cols, velocity);
        double elasticity_i = vel_i*vel_i;
        setVectorElement(elasticity, i, elasticity_i);
    }
}

double f(double x, double y, double t) {
    double x_source = 4787.0; // Middle of the Marmousi model
    double y_source = 0.0;    // Surface
    double f0 = 10.0;         // Central frequency in Hz
    double t0 = 1.0 / f0;     // Time delay
    
    // Distance from the source point
    double dist_sq = (x - x_source)*(x - x_source) + (y - y_source)*(y - y_source);
    
    // Ricker Wavelet formula
    double tau = t - t0;
    double ricker = (1.0 - 2.0 * M_PI * M_PI * f0 * f0 * tau * tau) * exp(-M_PI * M_PI * f0 * f0 * tau * tau);
    
    // Apply spatial distribution (Gaussian) so it's not just a single point
    return ricker * exp(-dist_sq / 100.0); 
}

void F_e(double xi, double nu, int e, double *output) {     // Assuming equidistant elements in x and y direction

    output[0] = xi*delta_x/2+getMatrixElement(center_element, e, 0);
    output[1] = nu*delta_y/2+getMatrixElement(center_element, e, 1);
}


// Initiates the connectivity matrix and the coordinates of the nodes in \Omega 
void init(int **connectivity, matrix *xy_points) {

    int grid_rows = n_y * (n_nodes - 1) + 1;
    int grid_cols = n_x * (n_nodes - 1) + 1;
    // printf(" rows = %d cols =%d \n", grid_rows, grid_cols);
    // Internal mapping grid
    int** auxiliary = (int**)malloc(grid_rows * sizeof(int*));
    for (int i = 0; i < grid_rows; i++) {
        auxiliary[i] = (int*)malloc(grid_cols * sizeof(int));
        for (int j = 0; j < grid_cols; j++) auxiliary[i][j] = -1;
    }

    // Counter nodes counts how many nodes there are, floor the vertical position of the element and door the horiziontal
    // For example, in a grid with n_x=3, n_y=2, element 5 is in floor 1 and door 2 (multiplied by n_nodes-1 to account for row and column)
    int counter = 0;
    int counter_nodes = 0;
    int floor = 0;
    int door = 0;

    // The ordering of the elements is from right to left, top to bottom
    // -------------------------- 0  1  2 ----------------------------
    // -------------------------- 3  4  5 ----------------------------
    // And inside each element the order of the nodes is from top to bottom, left to right
    // ------------------ 0=(0,0)  3=(1,0)  6=(2,0) ------------------
    // ------------------ 1=(0,1)  4=(1,1)  7=(2,1) ------------------
    // ------------------ 2=(0,2)  5=(1,2)  8=(2,2) ------------------

    for(int e = 0; e<n_elements; e++){
        if(counter % n_x == 0 && counter != 0){
            floor = floor + n_nodes - 1;
            door = 0;
        }
        for(int i = 0; i<n_nodes; i++){
            for(int j = 0; j<n_nodes; j++){
                int index_x = j + floor;
                int index_y = i + door;
                if(auxiliary[index_x][index_y] == -1){
                    auxiliary[index_x][index_y] = counter_nodes;
                    counter_nodes++;
                }
                connectivity[i*n_nodes+j][e] = auxiliary[index_x][index_y];
                
                double xi = getVectorElement(gll_points, i);
                double nu = getVectorElement(gll_points, n_nodes - j - 1);
                double point[2];
                F_e(xi, nu, e, point);
                // Check if nodes that share position also share coordinates (Debug)
                // if(fabs(point[0]-getMatrixElement(xy_points, auxiliary[index_x][index_y], 0))>1e-9 || fabs(point[1]-getMatrixElement(xy_points, auxiliary[index_x][index_y], 1))>1e-9){
                //     printf("   it changed!   ");
                // }
                setMatrixElement(xy_points, auxiliary[index_x][index_y], 0, point[0]);
                setMatrixElement(xy_points, auxiliary[index_x][index_y], 1, point[1]);
                // printf("index (%d, %d) element %d with node id %d and coordinates (%lf, %lf)\n "
                //         ,i, j, e, auxiliary[index_x][index_y], point[0], point[1]);
            }
        }
        door = door + n_nodes - 1;
        counter++;
    }

    // Debug 
    // for(int i = 0; i<grid_rows; i++){
    //     for(int j = 0; j<grid_cols; j++){
    //         printf("  %d  ", auxiliary[i][j]);
    //     }
    //     printf("\n");
    // }

    for (int i = 0; i < grid_rows; i++) free(auxiliary[i]);
    free(auxiliary);
}


double T(int m, int n, int i, int j, int l, int k) {

    double G11, G12, G22;
    G11 = 4/(delta_x*delta_x);
    G22 = 4/(delta_y*delta_y);
    G12 = 0;
    // Delta kroeneeckers 
    double d_il = (i == l) ? 1.0 : 0.0;
    double d_jk = (j == k) ? 1.0 : 0.0;
    double d_ml = (m == l) ? 1.0 : 0.0;
    double d_nk = (n == k) ? 1.0 : 0.0;

    // Derivative values from the D matrix: D[evaluation_point][basis_index]
    // l'_i(xi_l) is the slope of basis i at point l
    double lp_il = getMatrixElement(D, l, i);
    double lp_jk = getMatrixElement(D, k, j);
    double lp_ml = getMatrixElement(D, l, m);
    double lp_nk = getMatrixElement(D, k, n);

    // Term 1: G11 * l'_i(xi_l) * l_j(nu_k) * l'_m(xi_l) * l_n(nu_k)
    double term1 = G11 * (lp_il * d_jk * lp_ml * d_nk);

    // Term 2: G12 * (l'_i * l_j * l_m * l'_n + l_i * l'_j * l'_m * l_n)
    double term2 = G12 * ( (lp_il * d_jk * d_ml * lp_nk) + 
                           (d_il * lp_jk * lp_ml * d_nk) );

    // Term 3: G22 * l_i * l'_j * l_m * l'_n
    double term3 = G22 * (d_il * lp_jk * d_ml * lp_nk);

    return term1 + term2 + term3;
}

// void initial_conditions(matrix* u, vector* u_vel) {

//     for (int i=0; i<total_points; i++){
//         double x = getMatrixElement(xy_points, i, 0);
//         double y = getMatrixElement(xy_points, i, 1);
//         double init_u = 0;
//         double init_u_vel = 0;
//         setMatrixElement(u, 0, i, init_u);
//         setVectorElement(u_vel, i, init_u_vel);
//     }
//     return;
// }

// ------------------------------------------Droplet--------------------------------------
void initial_conditions(matrix* u, vector* u_vel) {
    double cx = 5.0, cy = 5.0; // Center of the pulse
    double sigma = 0.3;        // Width of the pulse
    double amplitude = 1.0;

    for (int i = 0; i < total_points; i++) {
        double x = getMatrixElement(xy_points, i, 0);
        double y = getMatrixElement(xy_points, i, 1);
        
        // Gaussian pulse formula: A * exp(-(dist^2)/(2*sigma^2))
        double dist_sq = (x - cx) * (x - cx) + (y - cy) * (y - cy);
        double val = amplitude * exp(-dist_sq / (2 * sigma * sigma));
        
        setMatrixElement(u, 0, i, 0);    // Initial displacement (u_0)
        // setMatrixElement(u, 0, i, val);    // Assume u_1 is same for zero velocity start
        setVectorElement(u_vel, i, 0.0);   // Initial velocity is zero
    }
}
// ---------------------------------------With known analytical sol-------------------------------------------------
// void initial_conditions(matrix* u, vector* u_vel) {
//     for (int i = 0; i < total_points; i++) {
//         double x = getMatrixElement(xy_points, i, 0);
//         double y = getMatrixElement(xy_points, i, 1);
//         double val = analytical_sol(x, y, 0.0);
//         setMatrixElement(u, 0, i, val);
//         setVectorElement(u_vel, i, 0.0); // Velocity is 0 at t=0 for cos(wt)
//     }
// }

double analytical_sol(double x, double y, double t) {
    double kx = M_PI / 10.0;
    double ky = M_PI / 6.0;
    double omega = sqrt(kx*kx + ky*ky); // for c=1
    return sin(kx * x) * sin(ky * y) * cos(omega * t);
}

double K(int i, int j, int m, int n, int e) {

    double sum_ij = 0.0;
    
    double jacobian = (delta_x * delta_y) / 4.0;
    double K_ij_mn_e = 0.0;
    for (int l = 0; l < n_nodes; l++) {
        for (int k = 0; k < n_nodes; k++) {
            
            // We only compute if the result isn't obviously zero
            // Term 1 needs k == j == n
            // Term 3 needs l == i == m
            // Term 2 (if G12 != 0) needs cross-matches
            
            double val_T = T(m, n, i, j, l, k);
            
            if (val_T != 0.0) {
                double w_lk = getVectorElement(gll_weights, l) * getVectorElement(gll_weights, k);
                int global_id_lk = connectivity[l*n_nodes+k][e];
                double x = getMatrixElement(xy_points, global_id_lk, 0);
                double y = getMatrixElement(xy_points, global_id_lk, 1);
                double mu_lk = getVectorElement(elasticity, global_id_lk);
                K_ij_mn_e += w_lk * jacobian * mu_lk * val_T;
            }
        }
    }
    return K_ij_mn_e;
}

// Computes the product of Ku
void comp_Ku(vector* Ku, matrix* u) {
    for(int e = 0; e<n_elements; e++){
        for(int i = 0; i<n_nodes; i++){
            for(int j = 0; j<n_nodes; j++){

                double sum_ij = 0;
                for(int m = 0; m<n_nodes; m++){
                    for(int n = 0; n<n_nodes; n++){

                        int global_id_mn = connectivity[m*n_nodes+n][e];
                        double u_mn = getMatrixElement(u, 0, global_id_mn);
                        double K_ij_mn_e = getMatrixElement(local_K[e], i*n_nodes+j, m*n_nodes+n);
                        sum_ij = sum_ij + K_ij_mn_e*u_mn;
                    }
                }
            int global_id_ij = connectivity[i*n_nodes+j][e];
            double K_old = getVectorElement(Ku, global_id_ij);
            setVectorElement(Ku, global_id_ij, K_old + sum_ij);
            }
        }
    }
    return;
}

void saveStepToFile(FILE *fp, matrix* u, int total_points) {
    if (fp == NULL) return;

    for (int j = 0; j < total_points; j++) {
        // Use %e for scientific notation (better for small wave amplitudes)
        fprintf(fp, "%.7e%s", getMatrixElement(u, 0, j), (j == total_points- 1) ? "" : ",");
    }
    fprintf(fp, "\n"); // New line for the next time step
    return; 
}


// To load Marmousi data

matrix* load_from_xyz(const char* filename, int* out_rows, int* out_cols, double* out_dx, double* out_dy) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open %s\n", filename);
        return NULL;
    }

    double x, y, v;
    double first_x, first_y, second_y, second_x = -1.0;
    int rows = 0;
    int total_points = 0;
    bool second_x_found = false;

    // --- Phase 1: Pre-scan for Dimensions and Spacing ---
    // Read first point
    if (fscanf(fp, "%lf %lf %lf", &first_x, &first_y, &v) != 3) return NULL;
    rows = 1;
    total_points = 1;

    // Read second point to get dy
    if (fscanf(fp, "%lf %lf %lf", &x, &y, &v) == 3) {
        *out_dy = fabs(y - first_y);
        total_points++;
        if (x == first_x) rows++;
    }

    // Continue scanning to find when X changes (to get rows) and total points
    while (fscanf(fp, "%lf %lf %lf", &x, &y, &v) == 3) {
        total_points++;
        if (!second_x_found && x != first_x) {
            second_x = x;
            *out_dx = fabs(second_x - first_x);
            second_x_found = true;
        }
        if (x == first_x) {
            rows++;
        }
    }

    int cols = total_points / rows;
    *out_rows = rows;
    *out_cols = cols;

    // --- Phase 2: Allocate and Load Data ---
    rewind(fp);
    matrix* vel_matrix = zeroMatrix(rows, cols); // From linear-algebra.h

    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (fscanf(fp, "%lf %lf %lf", &x, &y, &v) == 3) {
                // Mapping: i is row (Y/depth), j is column (X)
                // Accessing raw data pointer for efficiency as discussed
                vel_matrix->data[i * cols + j] = v;
            }
        }
    }

    fclose(fp);
    return vel_matrix;
}

// Calculates velocity using bilinear interpolation
double calc_velocity(double x, double y, int rows, int cols, matrix *velocity_grid) {

    int col_id = (int)(x/25);
    int row_id = (int)(y/25);

    // printf("Row id %d col id %d\n", row_id, col_id);
    double a = x/25 - col_id;
    double b = y/25 - row_id;

    // printf("a is %lf and b is %lf\n", a, b);
    double v_bottomleft = 0;
    double v_topright = 0;
    double v_bottomright = 0;

    int check = 0;
    if (row_id != rows - 1) {
        v_bottomleft = getMatrixElement(velocity_grid, row_id + 1, col_id);
        check++;
    }

    if (col_id != cols - 1) {
        v_topright = getMatrixElement(velocity_grid, row_id, col_id + 1);
        check++;
    }

    if (check == 2) {
        v_bottomright = getMatrixElement(velocity_grid, row_id + 1, col_id + 1);
    }

    double v_topleft = getMatrixElement(velocity_grid, row_id, col_id);
    // printf("vtl %lf, vtr %lf, vbl %lf, vbr %lf\n", v_topleft, v_topright, v_bottomleft, v_bottomright);

    return v_topleft*(1 - a)*(1 - b) + v_bottomleft*(1 - a)*b + v_topright*a*(1 - b) + v_bottomright*a*b;
}