import numpy as np
import matplotlib.pyplot as plt
import os

# ==========================================
# --- TOGGLES ---
# ==========================================
SHOW_GLOBAL_IDS = True   # Set to False to hide red Global IDs
SHOW_LOCAL_IDS = True    # Set to False to hide green Local IDs

def visualize_mesh(filename):
    if not os.path.exists(filename):
        print(f"Error: Could not find {filename}")
        return

    with open(filename, 'r') as f:
        lines = f.readlines()
        
    # Parse header
    num_elements, num_nodes = map(int, lines[0].split())
    print(f"Reading {num_elements} elements and {num_nodes} nodes...")
    
    # Parse nodes (0-based indexing)
    nodes = np.zeros((num_nodes, 2)) 
    for i in range(num_nodes):
        x, y = map(float, lines[i + 1].split())
        nodes[i, 0] = x
        nodes[i, 1] = y
        
    # Parse connectivity
    elements = []
    start_line = num_nodes + 1
    for i in range(num_elements):
        element_nodes = list(map(int, lines[start_line + i].split()))
        elements.append(element_nodes)
        
    # --- PLOTTING ---
    plt.figure(figsize=(12, 10))
    
    p_squared = len(elements[0])
    p = int(np.sqrt(p_squared))
    
    # Dictionary to collect local IDs that share the same global node
    local_id_map = {i: [] for i in range(num_nodes)}
    
    for el_idx, el in enumerate(elements):
        grid = np.array(el).reshape((p, p))
        
        # Draw internal grid lines
        for i in range(p):
            line_x = [nodes[node_id, 0] for node_id in grid[i, :]]
            line_y = [nodes[node_id, 1] for node_id in grid[i, :]]
            plt.plot(line_x, line_y, color='gray', linewidth=0.5, zorder=1)
            
        for j in range(p):
            line_x = [nodes[node_id, 0] for node_id in grid[:, j]]
            line_y = [nodes[node_id, 1] for node_id in grid[:, j]]
            plt.plot(line_x, line_y, color='gray', linewidth=0.5, zorder=1)
            
        # Draw element boundaries
        boundary_i = [grid[0, :], grid[-1, :]]
        boundary_j = [grid[:, 0], grid[:, -1]]
        
        for b_nodes in boundary_i + boundary_j:
            b_x = [nodes[node_id, 0] for node_id in b_nodes]
            b_y = [nodes[node_id, 1] for node_id in b_nodes]
            plt.plot(b_x, b_y, color='black', linewidth=1.5, zorder=3)

        # Collect Local IDs instead of plotting them immediately
        if SHOW_LOCAL_IDS:
            for local_id, global_id in enumerate(el):
                # Format: "e0:L5" (Element 0, Local Node 5)
                local_id_map[global_id].append(f'e{el_idx}:L{local_id}')

    # Plot the GLL points
    plt.scatter(nodes[:, 0], nodes[:, 1], color='blue', s=1, zorder=4)

    # Plot the collected Local IDs
    if SHOW_LOCAL_IDS:
        for global_id, local_labels in local_id_map.items():
            if local_labels:
                # Join multiple labels with a newline so they stack cleanly
                combined_label = '\n'.join(local_labels)
                plt.text(nodes[global_id, 0], nodes[global_id, 1], 
                         f' {combined_label}', color='green', fontsize=6, 
                         ha='left', va='bottom', zorder=5)

    # Plot Global IDs
    if SHOW_GLOBAL_IDS:
        for i in range(num_nodes):
            # Offset slightly to the bottom-left
            plt.text(nodes[i, 0], nodes[i, 1], 
                     f'G{i} ', color='red', fontsize=6, 
                     ha='right', va='top', zorder=5)

    plt.title(f"High-Order SEM Mesh\n{num_elements} Elements, {num_nodes} Nodes (p={p})")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.axis('equal') 
    plt.savefig("plotted_mesh.png")

if __name__ == '__main__':
    visualize_mesh("high_order_mesh.dat")