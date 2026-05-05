#include "Linear-Algebra-C/linear-algebra.h"
#include "polylib.h"
#include <stdio.h>
#include <math.h>

// --- Global Variables ---
double time_int = 3.0;
int time_points = 15000;
double time_delta; 

double x_0 = 0.0, x_f = 10.0;
double L;
int n_elements = 20;
int n_nodes = 6;      // Nodes per element
int total_points;     // Will be calculated in main
vector* x_points = NULL;

// --- Function Prototypes ---
double density(double x);
double elasticity(double x);
double f(double x, double t);
double F_e(double xi, int e);

// Function to save results
void saveResults(matrix* u, vector* x_coords, int time_pts, int space_pts) {
    FILE *fp = fopen("simulation_results.csv", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // Write the x-coordinates as the first row (the header)
    for (int j = 0; j < space_pts; j++) {
        fprintf(fp, "%f%s", getVectorElement(x_coords, j), (j == space_pts - 1) ? "" : ",");
    }
    fprintf(fp, "\n");

    // Write each time step as a new row
    for (int i = 0; i < time_pts; i++) {
        for (int j = 0; j < space_pts; j++) {
            fprintf(fp, "%e%s", getMatrixElement(u, i, j), (j == space_pts - 1) ? "" : ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Results saved to simulation_results.csv\n");
}

int main() {
    // Initialize global variables that depend on other globals
    time_delta = time_int / time_points;
    total_points = n_elements * (n_nodes - 1) + 1;

    // Use the library to allocate memory

    matrix* u = zeroMatrix(time_points, total_points);
    vector* u_vel = zeroVector(total_points);
    x_points = nullVector(total_points);

    // --- Calculate Weights and Integration Points
    vector* gll_points = nullVector(n_nodes);
    vector* gll_weights = nullVector(n_nodes);

    // Call the Polylib function
    // Pass the internal 'data' pointers from our vector structs
    zwgll(gll_points->data, gll_weights->data, n_nodes);

    // Compute nodes in x coordinates
    for (int e = 0; e < n_elements; e++) {

        for (int i = 0; i < n_nodes; i++) {
            double xi = getVectorElement(gll_points, i);

            setVectorElement(x_points, e * (n_nodes - 1) + i, F_e(xi, e));
        }
    }
    
    // Initial condition
    L = x_f-x_0;
    for (int i = 0; i < total_points; i++) {
        double x_val = getVectorElement(x_points, i);
        setMatrixElement(u, 0, i, sin(M_PI * x_val / L));
    }
    

    // Mass and Stiffness must be square matrices of size total_points x total_points
    matrix* K = zeroMatrix(total_points, total_points);
    vector* M = zeroVector(total_points);
    vector* vector_f = zeroVector(total_points);

    // --- Printing Section ---
    // printf("--- Physical Grid Points (x_points) ---\n");
    // printVector(x_points, false);
    // printf("\n\n");

    // printf("--- Solution Matrix (u) ---\n");
    // printMatrix(u, false);
    // printf("\n");


    // Compute the derivatives of the Legendre Polynomials
    matrix* D = zeroMatrix(n_nodes, n_nodes);
    matrix* Dt = zeroMatrix(n_nodes, n_nodes); // Workspace for dgll
    Dglj(D->data, Dt->data, gll_points->data, n_nodes, 0.0, 0.0);


    printf("DEBUG: Starting M assembly\n");
    // Computation of M
    for (int e=0; e<n_elements; e++){
        double jacobian = (x_f-x_0)/(2*n_elements);      // For equidistant nodes
        double density_i;
        
        for(int i=0; i<n_nodes; i++){
            double xi_i = getVectorElement(gll_points, i);
            double w_i = getVectorElement(gll_weights, i);
            density_i = density(F_e(xi_i, e));
            
            // Calculation accounts for the possible overlapping

            double M_i = getVectorElement(M, e*(n_nodes-1)+i);
            double calc = density_i*w_i*jacobian;       
            setVectorElement(M, e*(n_nodes-1)+i, M_i+calc);
        }

    }

    printf("DEBUG: Starting K assembly\n");
    // Compute Stiffness Matrix K
    for(int e = 0; e < n_elements; e++) {
        double jacobian = (x_f - x_0) / (2.0 * n_elements);
        
        for(int p_row = 0; p_row < n_nodes; p_row++) {
            int global_row = e * (n_nodes - 1) + p_row; // Global index for row

            for(int p_col = 0; p_col < n_nodes; p_col++) {
                int global_col = e * (n_nodes - 1) + p_col; // Global index for col
                double sum = 0.0; 

                for(int i = 0; i < n_nodes; i++) {
                    double xi_i = getVectorElement(gll_points, i);
                    double w_i = getVectorElement(gll_weights, i);
                    double elasticity_i = elasticity(F_e(xi_i, e));

                    // D_ij is the derivative of basis j at node i
                    double d_phi_row = getMatrixElement(D, i, p_row);
                    double d_phi_col = getMatrixElement(D, i, p_col);

                    // Numerical integration: (E * dphi_i/dxi * dphi_j/dxi * w) / J
                    sum += (elasticity_i * d_phi_row * d_phi_col * w_i) / jacobian;
                }

                // Assemble into the GLOBAL matrix K
                double existing_K = getMatrixElement(K, global_row, global_col);
                setMatrixElement(K, global_row, global_col, existing_K + sum); 
            }
        }
    }
    
    // Time loop
    for(int t_i=0; t_i<time_points-1; t_i++){
        double t = t_i*time_delta;
        // printf("DEBUG: Entering Time Loop t_i=%d\n", t_i);
        fillVector(&vector_f, 0.0);
        // Comppute the vector f

        for(int e = 0; e<n_elements; e++){
            double jacobian = (x_f-x_0)/(2*n_elements);      // For equidistant nodes
            
            for(int i=0; i<n_nodes; i++){
                double xi_i = getVectorElement(gll_points, i);
                double w_i = getVectorElement(gll_weights, i);
                double f_eval_i = f(F_e(xi_i, e), t);
                
                // Calculation accounts for the possible overlapping

                double vector_f_i = getVectorElement(vector_f, e*(n_nodes-1)+i);
                double calc = f_eval_i*jacobian*w_i;       
                setVectorElement(vector_f, e*(n_nodes-1)+i, vector_f_i+calc);
            }
        }

        //  Case u_1
        // We first compute u_1 (first instant after 0) as it is done with a different approximation to \ddot u
        if(t_i == 0){
            for(int m=1; m<total_points-1; m++){
                double Ku_0_m = 0.0;
                for(int j = 0; j < total_points; j++) {
                    Ku_0_m += getMatrixElement(K, m, j) * getMatrixElement(u, 0, j);
                }

                double M_inv = 1.0/getVectorElement(M, m);
                double vector_f_i = getVectorElement(vector_f, m);
                double ai = (vector_f_i - Ku_0_m) * M_inv;     

                double sum = getMatrixElement(u, 0, m) + time_delta*getVectorElement(u_vel, m) + time_delta*time_delta/2*ai;
                setMatrixElement(u, 1, m, sum);
            }
        }

        // Other cases
        else{
            for (int d = 1; d < total_points - 1; d++) {
                // We manually calculate (K * u_prev) for the d-th node
                double Ku_d = 0.0;
                for (int j = 0; j < total_points; j++) {
                    // Dot product of d-th row of K and the (t_i-1) row of u
                    double K_val = getMatrixElement(K, d, j);
                    double u_val = getMatrixElement(u, t_i, j);
                    Ku_d += K_val * u_val;
                }

                // Now use Ku_d in your central difference formula
                double Mi_inv = 1.0 / getVectorElement(M, d);
                double Fi = getVectorElement(vector_f, d);
                double u_curr = getMatrixElement(u, t_i, d);
                double u_old  = getMatrixElement(u, t_i - 1, d);

                double accel = (Fi - Ku_d) * Mi_inv;
                double u_next = (time_delta * time_delta * accel) + (2.0 * u_curr) - u_old;

                setMatrixElement(u, t_i+1, d, u_next);
            }
        }
    }

    // printMatrix(u, true);



    saveResults(u, x_points, time_points, total_points);
    // Cleanup
    deleteVector(gll_points); 
    deleteVector(gll_weights);
    deleteVector(M);
    deleteVector(u_vel);
    deleteMatrix(K);
    deleteMatrix(u);
    deleteVector(x_points);
    deleteMatrix(D);
    deleteMatrix(Dt);
    deleteVector(vector_f);
    
    return 0;
}


double density(double x)
{
    return 1.0;
}

double elasticity(double x)
{
    return 1.0;
}

double f(double x, double t) {
    double L = 10.0;
    double spatial = sin(M_PI * x / L);
    double temporal = cos(t);
    double constant = pow(M_PI / L, 2) - 1.0;
    return constant * spatial * temporal;
}

double F_e(double xi, int e) {
    // Width of a single element
    double h_e = (x_f - x_0) / n_elements;
    
    // Physical start of element 'e'
    double x_start = x_0 + e * h_e;
    
    // Map xi (-1 to 1) to the physical element [x_start, x_start + h_e]
    // x = x_start + (xi + 1) * (h_e / 2)
    double x = x_start + (xi + 1.0) * (h_e / 2.0);
    
    return x;
}
