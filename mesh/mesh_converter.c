#include <stdio.h>
#include <stdlib.h>
#include "../Linear-Algebra-C/linear-algebra.h"
#include "../Linear-Algebra-C/polylib.h"
#include <math.h>

double l_pos(double x) {
    return (x+1)/2;
}

double l_neg(double x) {
    return -(x-1)/2;
}

void F_e(int e, double xi, double nu, double *output, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {

    // clockwise order starting from down-left corner

    double v1x = l_neg(xi)*l_neg(nu)*x1;
    double v1y = l_neg(xi)*l_neg(nu)*y1;

    double v2x = l_neg(xi)*l_pos(nu)*x2;
    double v2y = l_neg(xi)*l_pos(nu)*y2;

    double v3x = l_pos(xi)*l_pos(nu)*x3;
    double v3y = l_pos(xi)*l_pos(nu)*y3;

    double v4x = l_pos(xi)*l_neg(nu)*x4;
    double v4y = l_pos(xi)*l_neg(nu)*y4;

    output[0] = v1x + v2x + v3x + v4x;
    output[1] = v1y + v2y + v3y + v4y;
}

int equal(double x1, double y1, double x2, double y2) {

    double tol = 1e-6;
    double dist = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);

    if(dist > tol) {
        // not equal
        return 0;
    }
    else {
        return 1;
    }
}

// Function that checks if point is on boundary by checking if it lives on a line 
// that joins two boundary edges defined by Salome.
// Returns 1 if it lives, 0 if not.
int is_on_boundary(double x, double y, double x1, double y1, double x2, double y2) {
    double tol = 1e-6;

    // The check is done by computing the determinant between the vector defined by the vertices and the point and one vertex
    double dx = x2 - x1;
    double dy = y2 - y1;
    double length_sq = dx * dx + dy * dy;

    // Edge case: the boundary segment is actually just a single point
    if (length_sq < tol * tol) {
        if (fabs(x - x1) < tol && fabs(y - y1) < tol) return 1;
        return 0;
    }


    double cross_product = (x - x1) * dy - (y - y1) * dx;
    
    double distance = fabs(cross_product) / sqrt(length_sq);
    
    if (distance > tol) {
        return 0; // Point is too far from the line to be collinear
    }

    // Check if the point falls within the bounding box of the segment
    double min_x = (x1 < x2 ? x1 : x2) - tol;
    double max_x = (x1 > x2 ? x1 : x2) + tol;
    double min_y = (y1 < y2 ? y1 : y2) - tol;
    double max_y = (y1 > y2 ? y1 : y2) + tol;

    if (x >= min_x && x <= max_x && y >= min_y && y <= max_y) {
        return 1; // It's on the line and inside the segment bounds
    }

    return 0;
}

int main() {

    // Parameters
    int p = 6;

    // Open the file
    FILE *file = fopen("simple_rectangle.dat", "r");
    if (!file) {
        printf("Error: Could not open Mesh_3.dat. Make sure the file is in the same directory.\n");
        return 1;
    }

    int num_vertex, total_elements;
    
    // Read the first line (header)
    if (fscanf(file, "%d %d", &num_vertex, &total_elements) != 2) {
        printf("Error: Could not read the header.\n");
        fclose(file);
        return 1;
    }

    // 1. ALLOCATE AND READ COORDINATES MATRIX (num_vertex x 2)
    double **coords = (double **)malloc(num_vertex * sizeof(double *));
    for (int i = 0; i < num_vertex; i++) {
        coords[i] = (double *)malloc(2 * sizeof(double));
        
        int id;
        double x, y, z;
        fscanf(file, "%d %lf %lf %lf", &id, &x, &y, &z);
        
        // Store Y and Z only
        coords[id - 1][0] = y;
        coords[id - 1][1] = z;
    }

    // 2. ALLOCATE AND READ CONNECTIVITY MATRIX (4 x quad_count)
    // Conn is the connectivity matrix for the simple mesh
    // quad_count is the real number of elements (Salome creates elements with 2 vertices only, we eliminate them)

    int **conn = (int **)malloc(4 * sizeof(int *));
    for (int i = 0; i < 4; i++) {
        conn[i] = (int *)malloc(total_elements * sizeof(int));
    }

    int **edges = (int **)malloc(2 * sizeof(int *));
    for (int i = 0; i < 2; i++) {
        edges[i] = (int *)malloc(total_elements * sizeof(int));
    }


    int quad_count = 0;
    int edge_count = 0;

    for (int i = 0; i < total_elements; i++) {
        int id, type;
        fscanf(file, "%d %d", &id, &type);
        
        if (type == 102) {
            int n1, n2;
            fscanf(file, "%d %d", &n1, &n2);
            edges[0][edge_count] = n1 - 1;
            edges[1][edge_count] = n2 - 1;
            edge_count++;
        } 
        else if (type == 204) {
            // STORE QUADRILATERAL ELEMENTS
            int n1, n2, n3, n4;
            fscanf(file, "%d %d %d %d", &n1, &n2, &n3, &n4);
            conn[0][quad_count] = n1 - 1;
            conn[1][quad_count] = n2 - 1;
            conn[2][quad_count] = n3 - 1;
            conn[3][quad_count] = n4 - 1;
            quad_count++;
        }
    }
    printf("edge count = %d\n", edge_count);

    fclose(file);

    // Mask for the nodes that are on the boundary
    int *is_node_boundary = (int *)calloc(quad_count * p * p, sizeof(int));
    if (is_node_boundary == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    // initiate the GLL points
    double *z = (double *)malloc(p * sizeof(double));
    double *w = (double *)malloc(p * sizeof(double));
    if (z == NULL || w == NULL) {
        printf("Memory allocation failed for 1D arrays.\n");
        return 1;
    }

    zwgll(z, w, p);

    // We create the new connectivity matrix

    int **connectivity = (int **)malloc(p*p*sizeof(int *));

    if (connectivity == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    for (int i = 0; i < p*p; i++) {
        connectivity[i] = (int *)malloc(quad_count * sizeof(int));
    }

    double **aux_edges = (double **)malloc(quad_count*4*p*sizeof(double *));

    if (aux_edges == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    for (int i = 0; i < quad_count*4*p; i++) {
        aux_edges[i] = (double *)malloc(3 * sizeof(double));
    }

    double **xy_points = (double **)malloc(quad_count*p*p*sizeof(double *));

    if (xy_points == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    for (int i = 0; i < quad_count*p*p; i++) {
        xy_points[i] = (double *)malloc(2 * sizeof(double));
    }

    int aux_counter = 0;
    int node_counter = 0;
    for(int e = 0; e<quad_count; e++) {
        double v1[2], v2[2], v3[2], v4[2];
        int id1, id2, id3, id4;
        id1 = conn[0][e], id2 = conn[1][e], id3 = conn[2][e], id4 = conn[3][e];
        v1[0] = coords[id1][0], v1[1] = coords[id1][1];
        v2[0] = coords[id2][0], v2[1] = coords[id2][1];
        v3[0] = coords[id3][0], v3[1] = coords[id3][1];
        v4[0] = coords[id4][0], v4[1] = coords[id4][1];

        for(int i = 0; i<p; i++){
            double xi = z[i];
            for(int j = 0; j<p; j++){

                int local_id = i*p + j;
                int global_id; 
                int found = 0;
                double nu = z[j];
                double node[2];

                F_e(e, xi, nu, node, v1[0], v1[1], v2[0], v2[1], v3[0], v3[1], v4[0], v4[1]);

                if(i == 0 || j == 0 || i == p-1 || j == p-1) {

                    for(int s = 0; s<aux_counter; s++) {
                        // See if the point has been explored before
                        int eq = equal(node[0], node[1], aux_edges[s][0], aux_edges[s][1]);
                        if (eq == 1) {
                            // Node has appeared before
                            global_id = aux_edges[s][2];
                            found = 1;
                            connectivity[local_id][e] = global_id; 
                            break;
                        }
                    }
                    if (found == 0){
                        // Node not found

                        global_id = node_counter; 
                        connectivity[local_id][e] = global_id; 
                        xy_points[node_counter][0] = node[0], xy_points[node_counter][1] = node[1];
                        aux_edges[aux_counter][0] = node[0], aux_edges[aux_counter][1] = node[1], aux_edges[aux_counter][2] = global_id;
                        node_counter++;
                        aux_counter++;
                    }
                    
                    // Check to see if it is on the boundary
                    int boundary = 0;
                    for(int k = 0; k<edge_count; k++) {
                        int edge_id1 = edges[0][k];
                        int edge_id2 = edges[1][k];
                        boundary = boundary + is_on_boundary(node[0], node[1], coords[edge_id1][0], coords[edge_id1][1], 
                                                                        coords[edge_id2][0], coords[edge_id2][1]);
                        if (boundary == 1) {
                            is_node_boundary[global_id] = 1;
                        }
                    }
                }
                else {
                    global_id = node_counter; 
                    connectivity[local_id][e] = global_id; 
                    xy_points[node_counter][0] = node[0], xy_points[node_counter][1] = node[1];
                    node_counter++;
                }
            }
        }

    }

   // ==========================================
    // --- TEST/VERIFICATION PRINT ---
    // ==========================================
    printf("Total node count is %d and total nodes shared is %d \n", node_counter, aux_counter);

    int print_limit = 5; // Change to node_counter / quad_count to print everything

    // 1. Print the physical coordinates (xy_points)
    printf("\n--- High-Order Nodes (xy_points) ---\n");
    int max_nodes = (node_counter < print_limit) ? node_counter : print_limit;
    for (int i = 0; i < max_nodes; i++) {
        printf("Global Node %d: Y = %f, Z = %f\n", i, xy_points[i][0], xy_points[i][1]);
    }
    if (node_counter > print_limit) {
        printf("... (skipping the remaining %d nodes)\n", node_counter - print_limit);
    }

    // 2. Print the connectivity matrix
    printf("\n--- High-Order Connectivity Matrix ---\n");
    int max_quads = (quad_count < print_limit) ? quad_count : print_limit;
    for (int e = 0; e < max_quads; e++) {
        printf("Element %d mapping (Local ID -> Global ID):\n [", e);
        for (int local_id = 0; local_id < p * p; local_id++) {
            printf("%d", connectivity[local_id][e]);
            if (local_id < (p * p) - 1) {
                printf(", ");
            }
        }
        printf("]\n\n");
    }
    if (quad_count > print_limit) {
        printf("... (skipping the remaining %d elements)\n", quad_count - print_limit);
    }


    // ==========================================
    // --- OUTPUT TO .DAT FILE ---
    // ==========================================
    FILE *out_file = fopen("high_order_mesh.dat", "w");
    if (!out_file) {
        printf("Error: Could not create high_order_mesh.dat.\n");
    } else {
        // Line 1: Number of elements and total unique nodes
        fprintf(out_file, "%d %d\n", quad_count, node_counter);
        
        // Next lines: xy_points (coordinates)
        for (int i = 0; i < node_counter; i++) {
            fprintf(out_file, "%.15e %.15e\n", xy_points[i][0], xy_points[i][1]);
        }
        
// Final lines: Connectivity Matrix
        for (int e = 0; e < quad_count; e++) {
            for (int local_id = 0; local_id < p * p; local_id++) {
                fprintf(out_file, "%d", connectivity[local_id][e]);
                if (local_id < (p * p) - 1) fprintf(out_file, " ");
            }
            fprintf(out_file, "\n");
        }

        // --- Write Boundary Nodes ---
        int total_boundary_nodes = 0;
        for (int i = 0; i < node_counter; i++) {
            if (is_node_boundary[i] == 1) total_boundary_nodes++;
        }

        // Create a distinct header for the solver to search for
        fprintf(out_file, "BOUNDARY_NODES %d\n", total_boundary_nodes);
        
        // Print the IDs space-separated
        for (int i = 0; i < node_counter; i++) {
            if (is_node_boundary[i] == 1) {
                fprintf(out_file, "%d ", i);
            }
        }
        fprintf(out_file, "\n");
        printf("Successfully wrote high-order mesh to 'high_order_mesh.dat'\n");
    }

// ==========================================
    // --- 3. CLEANUP MEMORY ---
    // ==========================================
    
    // Original struct cleanup
    for (int i = 0; i < num_vertex; i++) {
        free(coords[i]);
    }
    free(coords);
    
    for (int i = 0; i < 4; i++) {
        free(conn[i]);
    }
    free(conn);

    // --- NEW: Free the edges matrix ---
    // Remember to only loop up to 2, since we fixed the allocation bug!
    for (int i = 0; i < 2; i++) {
        free(edges[i]);
    }
    free(edges);

    // 1D array cleanup
    free(z);
    free(w);

    // High-order matrices cleanup
    for (int i = 0; i < p * p; i++) {
        free(connectivity[i]);
    }
    free(connectivity);

    for (int i = 0; i < quad_count * 4 * p; i++) {
        free(aux_edges[i]);
    }
    free(aux_edges);

    for (int i = 0; i < quad_count * p * p; i++) {
        free(xy_points[i]);
    }
    free(xy_points);

    free(is_node_boundary);

    return 0;
}
