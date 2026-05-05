#include "Linear-Algebra-C/linear-algebra.h"
#include "polylib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

// --- Global Variables ---
double time_int = 10.0;
int time_points = 5000;
double time_delta; 
int MESH_MODE = 1; // 0 = Manual Rectangle, 1 = Read from .dat


double x_0 = 0.0, x_f = 10.0;
double y_0 = 0.0, y_f = 6.0;
double delta_x, delta_y;
double Lx, Ly;
int n_x = 20, n_y = 20;            // Number of intervals per axis
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

        fscanf(fp, "%d %d", &n_elements, &total_points);
        rewind(fp);

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

                // Do not calculate if the node is in the edges (Neumann conditions)
                if (is_boundary(m) == 1) {
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

                if (is_boundary(d) == 1) {
                    setMatrixElement(u[2], 0, d, 0.0); 
                }
                else{
                    double Ku_d = getVectorElement(Ku, d);

                    // Now we use Ku_d in the central difference formula
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

double elasticity(double x, double y) {
    return 1.0;
}

double f(double x, double y, double t) {
    return 0.0;
}

void F_e(double xi, double nu, int e, double *output) {     // Assuming equidistant elements in x and y direction

    output[0] = xi*delta_x/2+getMatrixElement(center_element, e, 0);
    output[1] = nu*delta_y/2+getMatrixElement(center_element, e, 1);
}


double jacobian_calc(int e, int i, int j) {

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
    
    if (fscanf(fp, "%d %d", &num_elements, &num_nodes) != 2) {
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
        fscanf(fp, "%lf %lf", &x, &y);
        
        setMatrixElement(xy_points, i, 0, x);
        setMatrixElement(xy_points, i, 1, y);
    }

    // Read the Connectivity Matrix

    for (int e = 0; e < num_elements; e++) {
        for (int local_id = 0; local_id < n_nodes * n_nodes; local_id++) {
            int global_id;
            fscanf(fp, "%d", &global_id);
            
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
                fscanf(fp, "%d", &b_node);
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
// void initial_conditions(matrix* u, vector* u_vel) {
//     double cx = 5.0, cy = 5.0; // Center of the pulse
//     double sigma = 0.3;        // Width of the pulse
//     double amplitude = 1.0;

//     for (int i = 0; i < total_points; i++) {
//         double x = getMatrixElement(xy_points, i, 0);
//         double y = getMatrixElement(xy_points, i, 1);
        
//         // Gaussian pulse formula: A * exp(-(dist^2)/(2*sigma^2))
//         double dist_sq = (x - cx) * (x - cx) + (y - cy) * (y - cy);
//         double val = amplitude * exp(-dist_sq / (2 * sigma * sigma));
        
//         setMatrixElement(u, 0, i, val);    // Initial displacement (u_0)
//         // setMatrixElement(u, 0, i, val);    // Assume u_1 is same for zero velocity start
//         setVectorElement(u_vel, i, 0.0);   // Initial velocity is zero
//     }
// }
// ---------------------------------------With known analytical sol-------------------------------------------------
void initial_conditions(matrix* u, vector* u_vel) {
    for (int i = 0; i < total_points; i++) {
        double x = getMatrixElement(xy_points, i, 0);
        double y = getMatrixElement(xy_points, i, 1);
        double val = analytical_sol(x, y, 0.0);
        setMatrixElement(u, 0, i, val);
        setVectorElement(u_vel, i, 0.0); // Velocity is 0 at t=0 for cos(wt)
    }
}

double analytical_sol(double x, double y, double t) {
    double kx = M_PI / 1.00601787386982e+01;
    double ky = M_PI / 6.03174075151332e+00;
    double omega = sqrt(kx*kx + ky*ky); // for c=1
    return sin(kx * x) * sin(ky * y) * cos(omega * t);
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
        // Use %e for scientific notation (better for small wave amplitudes)
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
