# Spectral Elements For Wave Equation

Code that uses spectral elements to solve numerically the wave equation in 1-D and 2-D. This project implements the Spectral Element Method (SEM) using Gauss-Lobatto-Legendre (GLL) nodes to achieve high-order accuracy for seismic and acoustic wave simulations.

---

# `one_d_sem`: 1-D Solver

Solves the acoustic wave equation in a 1-D domain.

## Key Simulation Variables
The following parameters define the physics and resolution of the 1-D simulation:

### ⏳ Time Parameters
| Variable | Value | Description |
| :--- | :--- | :--- |
| `time_int` | `3.0` | Total duration of the simulation (seconds). |
| `time_points` | `15000` | Number of time steps. |
| `time_delta` | `3.0/15000` | Calculated time step size ($\Delta t$). |

### 📏 Space & Mesh Parameters
| Variable | Value | Description |
| :--- | :--- | :--- |
| `x_0`, `x_f` | `0.0`, `10.0` | Start and end points of the physical domain. |
| `n_elements` | `20` | Number of spectral elements. |
| `n_nodes` | `6` | Number of GLL nodes per element (polynomial order $N=5$). |
| `total_points` | `101` | Global unique nodes ($n_{elements} \times (n_{nodes} - 1) + 1$). |

## Output & Plotting
* **Output:** Generates `simulation_results.csv`. The first row contains spatial coordinates, and subsequent rows contain displacement $u(x)$ for every time step.
* **Plotting:** Visualize results and compute errors against analytical solutions using `plots.py`.

---

# `two_d_sem_working.c`: Updated 2-D Solver

A flexible 2-D wave equation solver that supports both manual grid generation and custom mesh imports.

## ⚙️ Core Parameters

### ⏳ Time Discretization
| Variable | Value | Description |
| :--- | :--- | :--- |
| `time_int` | `10.0` | Total simulation time. |
| `time_points` | `5000` | Number of time steps. |
| `time_delta` | `0.002` | Calculated time step size ($\Delta t$). |

### 🗺️ Mesh Modes
The program features a `MESH_MODE` toggle to determine the domain geometry:
* **Mode 0 (Manual):** Generates a standard rectangular grid based on user-defined boundaries (`x_0`, `x_f`, `y_0`, `y_f`) and intervals (`n_x`, `n_y`).
* **Mode 1 (File Import):** Reads high-order mesh information (coordinates, connectivity, and boundary nodes) from `high_order_mesh.dat`.

## 🧪 Physics & Boundary Conditions
* **Boundary Conditions:** Supports Dirichlet conditions (displacement fixed at $0.0$) using a lookup table of boundary node IDs loaded from the mesh file.
* **Properties:** Default configuration uses constant density ($\rho = 1.0$) and elasticity ($\mu = 1.0$).

---

# 🛠️ Mesh Utilities & Data

### `mesh_converter.c`
This utility transforms standard linear quadrilateral meshes into high-order SEM meshes.
* **Input:** Accepts a `.dat` file (e.g., `simple_rectangle.dat`) containing vertex coordinates and linear element connectivity.
* **Process:** Interpolates physical coordinates for GLL nodes within each element and computes unique global IDs for shared nodes to ensure $C^0$ continuity.
* **Output:** Generates `high_order_mesh.dat`, which includes a `BOUNDARY_NODES` section for the solver.

### `mesh_plotter.py`
A Python script used to visualize and verify the generated high-order meshes. It plots element boundaries, internal GLL grids, and both local and global node IDs.

---

# 🌊 Marmousi Test Case: `two_d_sem_marmousi.c`

A specialized version of the solver designed to run simulations on the complex **Marmousi velocity model**.

* **Data Source:** Reads `marmousi.xyz` to extract spatial velocity information.
* **Interpolation:** Uses bilinear interpolation to map velocities from the XYZ grid to the high-order GLL nodes.
* **Source Term:** Implements a **Ricker Wavelet** source located at the surface (center of the model) to simulate seismic acquisition.

---

# Requirements

*   **Polylib:** Routines for polynomial integration, differentiation, and interpolation. Used to calculate GLL points. [Nektar++ Polylib](https://www.nektar.info/2nd_edition/Polylib.html)
*   **Linear-Algebra-C:** A library for basic and advanced matrix/vector operations in C. [GitHub Repository](https://github.com/barrettotte/Linear-Algebra-C)
