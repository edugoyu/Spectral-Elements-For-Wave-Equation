#include "Linear-Algebra-C/linear-algebra.h"
#include "polylib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
// #define NDEBUG

// --- Global Variables ---
double time_int = 20.0;
int time_points = 5000;
double time_delta; 
int MESH_MODE = 0; // 0 = Manual Rectangle, 1 = Read from .dat
int PML = 1; // PML can only be used in MESH_MODE = 0, it is present if PML = 1, and not present if PML = 0
int PML_layers = 8;
int boundary_conditions = 0;   // If 1, Dirichlet, if 0, Neumann
double c = 1.;

double x_0 = 0.0, x_f = 100.0;
double y_0 = 0.0, y_f = 100.0;
double delta_x, delta_y;
double Lx, Ly;
int n_x = 30, n_y = 30;            // Number of intervals per axis
int n_elements;
int n_nodes = 6;      // Nodes per element per direction
int total_points;     // Will be calculated in main, total nodes in Omega
matrix* xy_points = NULL;           // In row n, coordinates of point n in global ID
matrix* center_element = NULL;     // Row n contains the center of element n
matrix** local_K = NULL;
int* boundary_nodes = NULL;    // Lookup array: boundary_nodes[i] = 1 if boundary, 0 if not
int num_boundary_nodes = 0;    // Total count of boundary nodes

int** connectivity = NULL;

// Initialization of the GLL points and Derivatives
vector* gll_points = NULL;
vector* gll_weights = NULL;
matrix* D = NULL;
matrix* Dt = NULL;
vector* Ku = NULL;

// CPML Global Variables
matrix* psi_x = NULL;
matrix* psi_y = NULL;
matrix* zeta_x = NULL;
matrix* zeta_y = NULL;

vector* a_x = NULL;
vector* b_x = NULL;
vector* a_y = NULL;
vector* b_y = NULL;

vector* du_dx = NULL;
vector* du_dy = NULL;
vector* d2u_dx2 = NULL;
vector* d2u_dy2 = NULL;
vector* dpsi_x_dx = NULL; 
vector* dpsi_y_dy = NULL;
vector* node_multiplicity = NULL; 
vector* F_cpml_vec = NULL; // Holds the integrated weak-form CPML forces


// --- Function Prototypes ---
double density(double x, double y);
double elasticity(double x, double y);
double f(double x, double y, double t);
void F_e(double xi, double nu, int e, double *output);
void init(int **connectivity, matrix* xy);          // Initializes the connectivity matrix and the points of interpolation
void init_from_file(int **connectivity, matrix *xy_points, FILE *fp);
double T(int m, int n, int i, int j, int l, int k, int e);
void initial_conditions(matrix* u, vector* u_vel);
double K(int i, int j, int m, int n, int e);
void comp_Ku(vector* Ku, matrix* u);                // Computes the product Ku
void saveStepToFile(FILE *fp, matrix* u, int total_points);
double analytical_sol(double x, double y, double t);
double jacobian_calc(int e, int i, int j);
int is_boundary(int global_id);
double calculate_pml_sigma(double distance, double thickness, double wave_velocity);
void init_cpml_profiles();
void comp_F_cpml(vector* F_cpml);
void compute_spatial_derivatives(matrix* field, vector* d_dx, vector* d_dy, vector* d2_dx2, vector* d2_dy2, int compute_second_order);

int main()
{
    time_delta = time_int / (time_points - 1);
    FILE *results_fp = fopen("simulation_results_2d.csv", "w");

    // 1D GLL Arrays (Needed for both modes)
    gll_points = nullVector(n_nodes);
    gll_weights = nullVector(n_nodes);
    zwgll(gll_points->data, gll_weights->data, n_nodes);

    if (MESH_MODE == 0) {
        printf("Mode: Manual Grid Generation\n");

        Lx = x_f - x_0; Ly = y_f - y_0;
        delta_x = Lx / n_x; delta_y = Ly / n_y;

        if(PML == 1) {
            x_f = x_f + delta_x * PML_layers;
            x_0 = x_0 - delta_x * PML_layers;
            y_f = y_f + delta_y * PML_layers;
            n_x = n_x + 2 * PML_layers;
            n_y = n_y + PML_layers;
        }

        n_elements = n_x * n_y;
        total_points = (n_x * (n_nodes - 1) + 1) * (n_y * (n_nodes - 1) + 1);

        // Center element calculation (Manual only)
        center_element = zeroMatrix(n_elements, 2);
        int e = 0;
        for(int j = 0; j < n_y; j++){
            for(int i = 0; i < n_x; i++){
                setMatrixElement(center_element, e, 0, x_0 + delta_x/2 + i*delta_x);
                setMatrixElement(center_element, e, 1, y_f - delta_y/2 - j*delta_y);
                e++;
            }
        }

        // Allocate connectivity and points. It was not done with the linear algebra library as it stores doubles and not int.
        connectivity = (int**)malloc(n_nodes * n_nodes * sizeof(int*));
        for (int i = 0; i < n_nodes * n_nodes; i++) 
            connectivity[i] = (int*)malloc(n_elements * sizeof(int));
        xy_points = zeroMatrix(total_points, 2);

        init(connectivity, xy_points);

    } else {
        printf("Mode: Importing Mesh from .dat\n");
        FILE *fp = fopen("mesh/high_order_mesh.dat", "r");
        if (!fp) { printf("Error opening mesh file.\n"); return 1; }

        if (fscanf(fp, "%d %d %d", &n_elements, &total_points, &n_nodes) != 3) {
            fprintf(stderr, "Error: Failed to read mesh header\n");
            return 1;
        }
        rewind(fp);
        printf("Number of elements->%d, number of points->%d, number of nodes per direction->%d\n", n_elements, total_points, n_nodes);

        // Allocate connectivity and points based on file data
        connectivity = (int**)malloc(n_nodes * n_nodes * sizeof(int*));
        for (int i = 0; i < n_nodes * n_nodes; i++) 
            connectivity[i] = (int*)malloc(n_elements * sizeof(int));
        xy_points = zeroMatrix(total_points, 2);

        init_from_file(connectivity, xy_points, fp);
        fclose(fp);
    }

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
        double density_ij;
        for(int i=0; i<n_nodes; i++){
            for(int j=0; j<n_nodes; j++){
                double jacobian = jacobian_calc(e, i, j);
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

    if (PML == 1) {
        printf("Initializing Weak-Form CPML Profiles...\n");
        init_cpml_profiles();
    }

    // Time loop 
    printf("Starting time steps\n");
    for(int t_i=0; t_i<time_points-1; t_i++){
        double t = t_i*time_delta;
        // printf("DEBUG: Entering Time Loop t_i=%d\n", t_i);
        fillVector(&vector_f, 0.0);
        fillVector(&Ku, 0.0);     

        // Computation of vector f
        for(int e = 0; e<n_elements; e++){
            
            for(int i=0; i<n_nodes; i++){
                for(int j=0; j<n_nodes; j++){
                        double jacobian = jacobian_calc(e, i, j);
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

                // Boundary conditions
                if (is_boundary(m) * boundary_conditions == 1) {
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
        else {
            // Compute Ku
            comp_Ku(Ku, u[1]);

            if (PML == 1) {
                // Calculate inner domain limits once for the loop
                double L_pml_x = PML_layers * delta_x;
                double L_pml_y = PML_layers * delta_y;
                double x_inner_min = x_0 + L_pml_x;
                double x_inner_max = x_f - L_pml_x;
                double y_inner_max = y_f - L_pml_y;

                compute_spatial_derivatives(u[1], du_dx, du_dy, d2u_dx2, d2u_dy2, 1);

                //  Update Psi
                for (int m = 0; m < total_points; m++) {
                    double x = getMatrixElement(xy_points, m, 0);
                    double y = getMatrixElement(xy_points, m, 1);
                    
                    // Skip inner domain nodes completely
                    if (x >= x_inner_min && x <= x_inner_max && y <= y_inner_max) continue;

                    double px_new = getVectorElement(a_x, m) * getMatrixElement(psi_x, 0, m) + getVectorElement(b_x, m) * getVectorElement(du_dx, m);
                    double py_new = getVectorElement(a_y, m) * getMatrixElement(psi_y, 0, m) + getVectorElement(b_y, m) * getVectorElement(du_dy, m);
                    
                    setMatrixElement(psi_x, 0, m, px_new);
                    setMatrixElement(psi_y, 0, m, py_new);
                }

                // Get the spatial derivative of Psi 
                compute_spatial_derivatives(psi_x, dpsi_x_dx, du_dy, d2u_dx2, d2u_dy2, 0); 
                compute_spatial_derivatives(psi_y, du_dx, dpsi_y_dy, d2u_dx2, d2u_dy2, 0);

                // 5. Update Zeta
                for (int d = 0; d < total_points; d++) {
                    double x = getMatrixElement(xy_points, d, 0);
                    double y = getMatrixElement(xy_points, d, 1);
                    
                    // OPTIMIZATION: Skip inner domain nodes completely
                    if (x >= x_inner_min && x <= x_inner_max && y <= y_inner_max) continue;

                    double zx_new = getVectorElement(a_x, d) * getMatrixElement(zeta_x, 0, d) + 
                                    getVectorElement(b_x, d) * (getVectorElement(d2u_dx2, d) + getVectorElement(dpsi_x_dx, d));
                    double zy_new = getVectorElement(a_y, d) * getMatrixElement(zeta_y, 0, d) + 
                                    getVectorElement(b_y, d) * (getVectorElement(d2u_dy2, d) + getVectorElement(dpsi_y_dy, d));
                    
                    setMatrixElement(zeta_x, 0, d, zx_new);
                    setMatrixElement(zeta_y, 0, d, zy_new);
                }

                // Integrate CPML forces into the Weak Form 
                comp_F_cpml(F_cpml_vec);

                // Calculate Acceleration and Advance Time
                for (int d = 0; d < total_points; d++) {
                    if (is_boundary(d) * boundary_conditions == 1) {
                        setMatrixElement(u[2], 0, d, 0.0); 
                    } else {
                        double Ku_d = getVectorElement(Ku, d);
                        double F_cpml_d = getVectorElement(F_cpml_vec, d); 
                        double M_inv_d = getVectorElement(M_inv, d);
                        double F_d = getVectorElement(vector_f, d);
                        
                        double accel = (F_d - Ku_d + F_cpml_d) * M_inv_d;
                        
                        double u_curr = getMatrixElement(u[1], 0, d);
                        double u_old  = getMatrixElement(u[0], 0, d);
                        double u_next = (time_delta * time_delta * accel) + (2.0 * u_curr) - u_old;

                        setMatrixElement(u[2], 0, d, u_next);
                    }
                }
            } 
            else {
                // --- STANDARD SOLVER (PML OFF) ---
                for (int d = 0; d < total_points; d++) {
                    if (is_boundary(d) * boundary_conditions == 1) {
                        setMatrixElement(u[2], 0, d, 0.0); 
                    } else {
                        double Ku_d = getVectorElement(Ku, d);
                        double M_inv_d = getVectorElement(M_inv, d);
                        double F_d = getVectorElement(vector_f, d);
                        double u_curr = getMatrixElement(u[1], 0, d);
                        double u_old  = getMatrixElement(u[0], 0, d);

                        double accel = (F_d - Ku_d) * M_inv_d;
                        double u_next = (time_delta * time_delta * accel) + (2.0 * u_curr) - u_old;

                        setMatrixElement(u[2], 0, d, u_next);
                    }
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


    fclose(results_fp);
    printf("Simulation done!\n");

    // Free matrices, vectors, etc.
    for (int i = 0; i < n_nodes * n_nodes; i++) {
        free(connectivity[i]);
    }
    for (int e = 0; e < n_elements; e++) {
        deleteMatrix(local_K[e]); 
    }  
    for(int i = 0; i<3; i++){
        deleteMatrix(u[i]);
    }
    free(connectivity);
    if (center_element != NULL) {
        deleteMatrix(center_element);
    }
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
    free(boundary_nodes);  
    if (PML == 1) {
        deleteMatrix(psi_x);
        deleteMatrix(psi_y);
        deleteMatrix(zeta_x);
        deleteMatrix(zeta_y);

        deleteVector(a_x);
        deleteVector(b_x);
        deleteVector(a_y);
        deleteVector(b_y);
        
        deleteVector(du_dx);
        deleteVector(du_dy);
        deleteVector(d2u_dx2);
        deleteVector(d2u_dy2);
        
        deleteVector(dpsi_x_dx);
        deleteVector(dpsi_y_dy);
        
        deleteVector(node_multiplicity);
        deleteVector(F_cpml_vec);
    }

    return 0;
}

double density(double x, double y) {
    return 1.0;
}


double elasticity(double x, double y) {
    return 1.0;
}

// Ricker Wavelet
double f(double x, double y, double t) {
    // --- 1. Source Parameters ---
    double f_p = 25.0;            // Peak frequency of the Ricker wavelet (Hz)
    double t_d = 1.2 / f_p;       // Time delay to make the source causal (starts at 0)
    double amplitude = 100.0;     // Overall strength of the source

    // --- 2. Target Location ---
    // x_s: Exactly halfway between the left and right boundaries
    // y_s: Placed at the surface (y_0)
    double x_s = (x_0 + x_f) / 2.0; 
    double y_s = y_0 + 0.1 * delta_y;             
    
    // Spatial width of the source (sigma)
    double sigma = delta_x;       

    // --- 3. Time Component (Ricker Wavelet) ---
    // Formula: (1 - 2 * (pi * f_p * (t - t_d))^2) * exp(-(pi * f_p * (t - t_d))^2)
    double arg = M_PI * f_p * (t - t_d);
    double arg_sq = arg * arg;
    double ricker_time = (1.0 - 2.0 * arg_sq) * exp(-arg_sq*10);

    // --- 4. Spatial Component (Gaussian Distribution) ---
    // Distributes the source energy over a small local area to prevent ringing
    double dist_sq = (x - x_s) * (x - x_s) + (y - y_s) * (y - y_s);
    double spatial_dist = exp(-dist_sq / (2.0 * sigma * sigma));

    // Combine space, time, and amplitude
    return amplitude * spatial_dist * ricker_time;
}

void F_e(double xi, double nu, int e, double *output) {     // Assuming equidistant elements in x and y direction

    output[0] = xi*delta_x/2+getMatrixElement(center_element, e, 0);
    output[1] = nu*delta_y/2+getMatrixElement(center_element, e, 1);
}


double jacobian_calc(int e, int i, int j) {

    if(MESH_MODE == 0) {
        return delta_x * delta_y / 4;
    }

    double xi = getVectorElement(gll_points, i);
    double nu = getVectorElement(gll_points, j);

    double l_1_xi = (xi + 1)/2;
    double l_n1_xi = -(xi - 1)/2;
    double l_1_nu = (nu + 1)/2;
    double l_n1_nu = -(nu - 1)/2;

    int global_id1 = connectivity[0][e];
    int global_id2 = connectivity[n_nodes - 1][e];
    int global_id3 = connectivity[n_nodes * n_nodes -1][e];
    int global_id4 = connectivity[n_nodes*(n_nodes - 1)][e];

    double v1[2], v2[2], v3[2], v4[2];
    v1[0] = getMatrixElement(xy_points, global_id1, 0), v1[1] = getMatrixElement(xy_points, global_id1, 1);
    v2[0] = getMatrixElement(xy_points, global_id2, 0), v2[1] = getMatrixElement(xy_points, global_id2, 1);
    v3[0] = getMatrixElement(xy_points, global_id3, 0), v3[1] = getMatrixElement(xy_points, global_id3, 1);
    v4[0] = getMatrixElement(xy_points, global_id4, 0), v4[1] = getMatrixElement(xy_points, global_id4, 1);

    double dx_dxi = -0.5 * l_n1_nu * v1[0] - 0.5  * l_1_nu * v2[0] + 0.5  * l_1_nu * v3[0]+ 1./2 * l_n1_nu * v4[0];
    double dy_dxi = -0.5  * l_n1_nu * v1[1] - 0.5  * l_1_nu * v2[1] + 0.5  * l_1_nu * v3[1]+ 0.5  * l_n1_nu * v4[1];

    double dx_dnu = -0.5  * l_n1_xi * v1[0] + 0.5  * l_n1_xi * v2[0] + 0.5  * l_1_xi * v3[0]- 0.5  * l_1_xi * v4[0];
    double dy_dnu = -0.5 * l_n1_xi * v1[1] + 0.5  * l_n1_xi * v2[1] + 0.5  * l_1_xi * v3[1]- 0.5  * l_1_xi * v4[1];

    double det_j = dx_dxi * dy_dnu - dy_dxi * dx_dnu;
    return fabs(det_j);

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

// Initiates the connectivity matrix and the coordinates of the nodes in \Omega from a file 
// Reads first number of nodes, 
void init_from_file(int **connectivity, matrix *xy_points, FILE *fp) {

    if (fp == NULL) {
        printf("Error: File pointer is null.\n");
        return;
    }

    int num_elements, num_nodes;
    
    if (fscanf(fp, "%d %d %d", &num_elements, &num_nodes, &n_nodes) != 3) {
        printf("Error: Could not read the header.\n");
        return;
    }

    // Verify the pre-allocated matrix matches the file data 
    if (xy_points->rows != num_nodes || xy_points->cols != 2) {
        printf("Warning: xy_points matrix dimensions (%dx%d) do not match file data (%dx2).\n", 
               xy_points->rows, xy_points->cols, num_nodes);
    }

    // Read the coordinates into the linear-algebra matrix
    // The file stores them as: X Y
    for (int i = 0; i < num_nodes; i++) {
        double x, y;
        if (fscanf(fp, "%lf %lf", &x, &y) != 2) {
                fprintf(stderr, "Error: Failed to read coordinates at node %d\n", i);
                // handle error
            }
        
        setMatrixElement(xy_points, i, 0, x);
        setMatrixElement(xy_points, i, 1, y);
    }

    // Read the Connectivity Matrix

    for (int e = 0; e < num_elements; e++) {
        for (int local_id = 0; local_id < n_nodes * n_nodes; local_id++) {
            int global_id;

            if (fscanf(fp, "%d", &global_id) != 1) {
                printf("Error: Could not read global id.\n");
                return;
            }
            
            connectivity[local_id][e] = global_id;
        }
    }

    char label[50];
    // Check if the file has the BOUNDARY_NODES header
    if (fscanf(fp, "%s %d", label, &num_boundary_nodes) == 2) {
        if (strcmp(label, "BOUNDARY_NODES") == 0) {
            // Allocate the lookup array initialized to 0
            boundary_nodes = (int*)calloc(num_nodes, sizeof(int));
            
            for (int i = 0; i < num_boundary_nodes; i++) {
                int b_node;

                if (fscanf(fp, "%d", &b_node) != 1) {
                        printf("Error: Could not read global id.\n");
                        return;
                    }
            
                boundary_nodes[b_node] = 1; // Flag this specific ID as a boundary
            }
            printf("Successfully loaded %d boundary nodes.\n", num_boundary_nodes);
        }
    } else {
        printf("Warning: No boundary nodes found in mesh file.\n");
    }
    
    printf("Successfully initialized %d nodes and %d elements from file.\n", num_nodes, num_elements);
}


double T(int m, int n, int i, int j, int l, int k, int e) {

    // Get the reference coordinates for the current integration point (l, k)
    double xi = getVectorElement(gll_points, l);
    double nu = getVectorElement(gll_points, k);

    // Fetch the 4 corners of the element to calculate the components of the jacobian
    int id1 = connectivity[0][e];
    int id2 = connectivity[n_nodes - 1][e];
    int id3 = connectivity[n_nodes * n_nodes - 1][e];
    int id4 = connectivity[n_nodes * (n_nodes - 1)][e];

    double x1 = getMatrixElement(xy_points, id1, 0), y1 = getMatrixElement(xy_points, id1, 1);
    double x2 = getMatrixElement(xy_points, id2, 0), y2 = getMatrixElement(xy_points, id2, 1);
    double x3 = getMatrixElement(xy_points, id3, 0), y3 = getMatrixElement(xy_points, id3, 1);
    double x4 = getMatrixElement(xy_points, id4, 0), y4 = getMatrixElement(xy_points, id4, 1);

    double l_1_xi = (xi + 1.0) / 2.0;   double l_n1_xi = -(xi - 1.0) / 2.0;
    double l_1_nu = (nu + 1.0) / 2.0;   double l_n1_nu = -(nu - 1.0) / 2.0;

    double dx_dxi = -0.5 * l_n1_nu * x1 - 0.5 * l_1_nu * x2 + 0.5 * l_1_nu * x3 + 0.5 * l_n1_nu * x4;
    double dy_dxi = -0.5 * l_n1_nu * y1 - 0.5 * l_1_nu * y2 + 0.5 * l_1_nu * y3 + 0.5 * l_n1_nu * y4;
    double dx_dnu = -0.5 * l_n1_xi * x1 + 0.5 * l_n1_xi * x2 + 0.5 * l_1_xi * x3 - 0.5 * l_1_xi * x4;
    double dy_dnu = -0.5 * l_n1_xi * y1 + 0.5 * l_n1_xi * y2 + 0.5 * l_1_xi * y3 - 0.5 * l_1_xi * y4;

    // Compute the signed determinant
    double signed_det = dx_dxi * dy_dnu - dy_dxi * dx_dnu;

    double xi_x = dy_dnu / signed_det;
    double xi_y = -dx_dnu / signed_det;
    double nu_x = -dy_dxi / signed_det;
    double nu_y = dx_dxi / signed_det;

    // Metric Tensor Components
    double G11 = xi_x * xi_x + xi_y * xi_y;
    double G22 = nu_x * nu_x + nu_y * nu_y;
    double G12 = xi_x * nu_x + xi_y * nu_y;


    //  T FUNCTION MATH

    double d_il = (i == l) ? 1.0 : 0.0;
    double d_jk = (j == k) ? 1.0 : 0.0;
    double d_ml = (m == l) ? 1.0 : 0.0;
    double d_nk = (n == k) ? 1.0 : 0.0;

    double lp_il = getMatrixElement(D, l, i);
    double lp_jk = getMatrixElement(D, k, j);
    double lp_ml = getMatrixElement(D, l, m);
    double lp_nk = getMatrixElement(D, k, n);

    double term1 = G11 * (lp_il * d_jk * lp_ml * d_nk);
    double term2 = G12 * ( (lp_il * d_jk * d_ml * lp_nk) + (d_il * lp_jk * lp_ml * d_nk) );
    double term3 = G22 * (d_il * lp_jk * d_ml * lp_nk);

    return term1 + term2 + term3;
}


void initial_conditions(matrix* u, vector* u_vel) {
        
    setMatrixElement(u, 0, i, 0);    // Initial displacement (u_0)
    setVectorElement(u_vel, i, 0.0);   // Initial velocity is zero
}

double K(int i, int j, int m, int n, int e) {

    double sum_ij = 0.0;
    
    
    double K_ij_mn_e = 0.0;
    for (int l = 0; l < n_nodes; l++) {
        for (int k = 0; k < n_nodes; k++) {
            double jacobian = jacobian_calc(e, l, k);
            // We only compute if the result isn't obviously zero
            // Term 1 needs k == j == n
            // Term 3 needs l == i == m
            // Term 2 (if G12 != 0) needs cross-matches
            
            double val_T = T(m, n, i, j, l, k, e);
            
            if (val_T != 0.0) {
                double w_lk = getVectorElement(gll_weights, l) * getVectorElement(gll_weights, k);
                int global_id_lk = connectivity[l*n_nodes+k][e];
                double x = getMatrixElement(xy_points, global_id_lk, 0);
                double y = getMatrixElement(xy_points, global_id_lk, 1);
                double mu_lk = elasticity(x, y);
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
        fprintf(fp, "%.7e%s", getMatrixElement(u, 0, j), (j == total_points- 1) ? "" : ",");
    }
    fprintf(fp, "\n"); // New line for the next time step
    return; 
}

// Returns 1 if the node is on the boundary, 0 otherwise
int is_boundary(int global_id) {
    if (MESH_MODE == 0) {
        // Mode 0: Manual Rectangle limits
        double x = getMatrixElement(xy_points, global_id, 0);
        double y = getMatrixElement(xy_points, global_id, 1);
        double eps = 1e-9; 
        
        if (x <= x_0 + eps || x >= x_f - eps || y <= y_0 + eps || y >= y_f - eps) {
            return 1;
        }
        return 0;
    } else {
        // Mode 1: Check the lookup table populated from the .dat file
        if (boundary_nodes != NULL && boundary_nodes[global_id] == 1) {
            return 1;
        }
        return 0;
    }
}


// ----------------------------------------------   CPML FUNCTIONS   ----------------------------------------------------

// Damping Profile
double calculate_pml_sigma(double distance, double thickness, double wave_velocity) {
    double R = 1e-5; 
    double m = 2.0;  

    if (thickness <= 0.0 || distance <= 0.0) return 0.0;
    if (distance > thickness) distance = thickness;

    double sigma_max = (3.0 * wave_velocity) / (2.0 * thickness) * log(1.0 / R);
    return sigma_max * pow(distance / thickness, m);
}


// Initialization of CPML vectors and calculation of a_x and b_x
void init_cpml_profiles() {
    psi_x = zeroMatrix(1, total_points);   psi_y = zeroMatrix(1, total_points);
    zeta_x = zeroMatrix(1, total_points);  zeta_y = zeroMatrix(1, total_points);
    a_x = zeroVector(total_points);     b_x = zeroVector(total_points);
    a_y = zeroVector(total_points);     b_y = zeroVector(total_points);
    
    du_dx = zeroVector(total_points);     du_dy = zeroVector(total_points);
    d2u_dx2 = zeroVector(total_points);   d2u_dy2 = zeroVector(total_points);
    dpsi_x_dx = zeroVector(total_points); dpsi_y_dy = zeroVector(total_points);
    node_multiplicity = zeroVector(total_points);
    F_cpml_vec = zeroVector(total_points);

    double L_pml_x = PML_layers * delta_x;
    double L_pml_y = PML_layers * delta_y;
    double x_inner_min = x_0 + L_pml_x;
    double x_inner_max = x_f - L_pml_x;
    double y_inner_max = y_f - L_pml_y;
    
    double v_p = sqrt(60.0 / 1.0);

    double alpha = 0.5; 

    for (int i = 0; i < total_points; i++) {
        double x = getMatrixElement(xy_points, i, 0);
        double y = getMatrixElement(xy_points, i, 1);
        
        double sig_x = 0.0, sig_y = 0.0;

        // X is active on both sides
        if (x < x_inner_min) sig_x = calculate_pml_sigma(x_inner_min - x, L_pml_x, v_p);
        else if (x > x_inner_max) sig_x = calculate_pml_sigma(x - x_inner_max, L_pml_x, v_p);

        // Y is only active at the upper limit (y_f)
        if (y > y_inner_max) sig_y = calculate_pml_sigma(y - y_inner_max, L_pml_y, v_p);

        double alp_x = (sig_x > 0) ? alpha : 0.0;
        double ax_val = exp(-(sig_x + alp_x) * time_delta);
        double bx_val = (sig_x > 0) ? (sig_x / (sig_x + alp_x)) * (ax_val - 1.0) : 0.0;
        setVectorElement(a_x, i, ax_val);
        setVectorElement(b_x, i, bx_val);

        double alp_y = (sig_y > 0) ? alpha : 0.0;
        double ay_val = exp(-(sig_y + alp_y) * time_delta);
        double by_val = (sig_y > 0) ? (sig_y / (sig_y + alp_y)) * (ay_val - 1.0) : 0.0;
        setVectorElement(a_y, i, ay_val);
        setVectorElement(b_y, i, by_val);
    }
}

// Weak formulation of the CPML function
void comp_F_cpml(vector* F_cpml) {
    fillVector(&F_cpml, 0.0);
    
    // Calculate boundaries of the domain
    double L_pml_x = PML_layers * delta_x;
    double L_pml_y = PML_layers * delta_y;
    double x_inner_min = x_0 + L_pml_x;
    double x_inner_max = x_f - L_pml_x;
    double y_inner_max = y_f - L_pml_y;
    
    for(int e = 0; e < n_elements; e++){

        // Get the coordinates of the first and last node to find the element's bounding box
        int id_first = connectivity[0][e];
        int id_last = connectivity[n_nodes * n_nodes - 1][e];
        
        double x_min = fmin(getMatrixElement(xy_points, id_first, 0), getMatrixElement(xy_points, id_last, 0));
        double x_max = fmax(getMatrixElement(xy_points, id_first, 0), getMatrixElement(xy_points, id_last, 0));
        double y_max = fmax(getMatrixElement(xy_points, id_first, 1), getMatrixElement(xy_points, id_last, 1));
        
        // If the entire element is safely inside the inner domain, skip it
        if (x_min >= x_inner_min && x_max <= x_inner_max && y_max <= y_inner_max) {
            continue; 
        }
        // If not, apply the math
        for(int i = 0; i < n_nodes; i++){
            for(int j = 0; j < n_nodes; j++){
                int global_id_ij = connectivity[i*n_nodes+j][e];
                double w_i = getVectorElement(gll_weights, i);
                double w_j = getVectorElement(gll_weights, j);
                double jacobian_ij = jacobian_calc(e, i, j);
                
                double x_ij = getMatrixElement(xy_points, global_id_ij, 0);
                double y_ij = getMatrixElement(xy_points, global_id_ij, 1);
                double mu_ij = elasticity(x_ij, y_ij); 
                
                double zx = getMatrixElement(zeta_x, 0, global_id_ij);
                double zy = getMatrixElement(zeta_y, 0, global_id_ij);
                
                double term1 = w_i * w_j * jacobian_ij * mu_ij * (zx + zy);
                
                double term2_x = 0.0;
                for (int l = 0; l < n_nodes; l++) {
                    int global_id_lj = connectivity[l*n_nodes+j][e];
                    double w_l = getVectorElement(gll_weights, l);
                    double jac_lj = jacobian_calc(e, l, j);
                    double x_lj = getMatrixElement(xy_points, global_id_lj, 0);
                    double y_lj = getMatrixElement(xy_points, global_id_lj, 1);
                    double mu_lj = elasticity(x_lj, y_lj);
                    double psi_x_lj = getMatrixElement(psi_x, 0, global_id_lj);
                    
                    term2_x += w_l * w_j * jac_lj * mu_lj * getMatrixElement(D, l, i) * (2.0 / delta_x) * psi_x_lj;
                }
                
                double term2_y = 0.0;
                for (int k = 0; k < n_nodes; k++) {
                    int global_id_ik = connectivity[i*n_nodes+k][e];
                    double w_k = getVectorElement(gll_weights, k);
                    double jac_ik = jacobian_calc(e, i, k);
                    double x_ik = getMatrixElement(xy_points, global_id_ik, 0);
                    double y_ik = getMatrixElement(xy_points, global_id_ik, 1);
                    double mu_ik = elasticity(x_ik, y_ik);
                    double psi_y_ik = getMatrixElement(psi_y, 0, global_id_ik);
                    
                    term2_y += w_i * w_k * jac_ik * mu_ik * getMatrixElement(D, k, j) * (2.0 / delta_y) * psi_y_ik;
                }
                
                // Assemble to global node
                double current_F = getVectorElement(F_cpml, global_id_ij);
                setVectorElement(F_cpml, global_id_ij, current_F + term1 - term2_x - term2_y);
            }
        }
    }
}

// Derivative Engine
void compute_spatial_derivatives(matrix* field, vector* d_dx, vector* d_dy, vector* d2_dx2, vector* d2_dy2, int compute_second_order) {
    fillVector(&d_dx, 0.0); fillVector(&d_dy, 0.0);
    if (compute_second_order) { fillVector(&d2_dx2, 0.0); fillVector(&d2_dy2, 0.0); }
    fillVector(&node_multiplicity, 0.0);

    for (int e = 0; e < n_elements; e++) {
        double xi_x = 2.0 / delta_x;
        double nu_y = 2.0 / delta_y; 

        // Boundary box optimization: skip inner domain
        int id_first = connectivity[0][e];
        int id_last = connectivity[n_nodes * n_nodes - 1][e];
        double x_inner_min = x_0 + (PML_layers * delta_x);
        double x_inner_max = x_f - (PML_layers * delta_x);
        double y_inner_max = y_f - (PML_layers * delta_y);
        
        double x_min = fmin(getMatrixElement(xy_points, id_first, 0), getMatrixElement(xy_points, id_last, 0));
        double x_max = fmax(getMatrixElement(xy_points, id_first, 0), getMatrixElement(xy_points, id_last, 0));
        double y_max = fmax(getMatrixElement(xy_points, id_first, 1), getMatrixElement(xy_points, id_last, 1));
        
        if (x_min >= x_inner_min && x_max <= x_inner_max && y_max <= y_inner_max) continue; 

        for (int i = 0; i < n_nodes; i++) {
            for (int j = 0; j < n_nodes; j++) {
                int global_id = connectivity[i * n_nodes + j][e];
                double local_d_dxi = 0.0, local_d_dnu = 0.0;
                double local_d2_dxi2 = 0.0, local_d2_dnu2 = 0.0;

                for (int m = 0; m < n_nodes; m++) {
                    double val_x = getMatrixElement(field, 0, connectivity[m * n_nodes + j][e]);
                    local_d_dxi += getMatrixElement(D, i, m) * val_x;

                    double val_y = getMatrixElement(field, 0, connectivity[i * n_nodes + m][e]);
                    local_d_dnu += getMatrixElement(D, j, m) * val_y;

                    if (compute_second_order) {
                        for (int k = 0; k < n_nodes; k++) {
                            // 2nd derivatives aligned with their respective axes
                            local_d2_dxi2 += getMatrixElement(D, i, m) * getMatrixElement(D, m, k) * getMatrixElement(field, 0, connectivity[k * n_nodes + j][e]);
                            local_d2_dnu2 += getMatrixElement(D, j, m) * getMatrixElement(D, m, k) * getMatrixElement(field, 0, connectivity[i * n_nodes + k][e]);
                        }
                    }
                }
                
                setVectorElement(d_dx, global_id, getVectorElement(d_dx, global_id) + local_d_dxi * xi_x);
                setVectorElement(d_dy, global_id, getVectorElement(d_dy, global_id) + local_d_dnu * nu_y);
                if (compute_second_order) {
                    setVectorElement(d2_dx2, global_id, getVectorElement(d2_dx2, global_id) + local_d2_dxi2 * (xi_x * xi_x));
                    setVectorElement(d2_dy2, global_id, getVectorElement(d2_dy2, global_id) + local_d2_dnu2 * (nu_y * nu_y));
                }
                setVectorElement(node_multiplicity, global_id, getVectorElement(node_multiplicity, global_id) + 1.0);
            }
        }
    }

    for (int i = 0; i < total_points; i++) {
        double mult = getVectorElement(node_multiplicity, i);
        if (mult > 0) {
            setVectorElement(d_dx, i, getVectorElement(d_dx, i) / mult);
            setVectorElement(d_dy, i, getVectorElement(d_dy, i) / mult);
            if (compute_second_order) {
                setVectorElement(d2_dx2, i, getVectorElement(d2_dx2, i) / mult);
                setVectorElement(d2_dy2, i, getVectorElement(d2_dy2, i) / mult);
            }
        }
    }
}
