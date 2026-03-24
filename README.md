# Spectral-Elements-For-Wave-Equation
Code that uses spectral elements to solve numerically the wave equation in 1-D and 2-D.

# `one_d_sem`: Solves the Wave equation in 1-D. 

## Key Simulation Variables

The following parameters in `main.c` define the physics and resolution of the simulation:

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

## Output

The program generates a file named `simulation_results.csv`.
* **First row:** Spatial $x$ coordinates.
* **Subsequent rows:** Displacement $u(x)$ for every time step.

## Plot
* Plotted using python with the file `plots.py`. Computes the error also between the analytical and the numerical solution.

## Physical Properties
The default configuration uses:
* **Density ($\rho$):** Constant 1.0.
* **Elasticity ($E$):** Constant 1.0.
* **Source Term ($f$):** Manufactured solution using a $sin(x)\cos(t)$ profile.

# `two_d_sem`: Solves the Wave equation in 2-D. 

## ⚙️ Core Parameters

### ⏳ Time Discretization
| Variable | Value | Description |
| :--- | :--- | :--- |
| `time_int` | `10.0` | Total simulation time. |
| `time_points` | `2000` | Number of time steps. |
| `time_delta` | `0.005` | Calculated time step size ($\Delta t$). |

### 🗺️ Spatial Domain & Mesh
| Variable | Value | Description |
| :--- | :--- | :--- |
| `x_0`, `x_f` | `0.0, 10.0` | Domain boundaries in X. |
| `y_0`, `y_f` | `0.0, 10.0` | Domain boundaries in Y. |
| `n_x`, `n_y` | `10, 10` | Intervals (elements) per axis (100 total elements). |
| `n_nodes` | `6` | GLL nodes per element per direction (Order $N=5$). |
| `total_points` | `2601` | Total unique global nodes. |

## 🧪 Physics & Initial Conditions

* **Initial State:** A Gaussian "droplet" pulse centered at $(5.0, 5.0)$ with $\sigma = 0.3$.
* **Boundary Conditions:** Dirichlet (fixed at $0.0$) implemented at the domain edges.
* **Properties:** Constant density ($\rho = 1.0$) and elasticity ($\mu = 1.0$).

## 🛠️ Data Structures

* **`connectivity`**: A 2D integer array mapping `[local_node_id][element_id]` to the `global_id`.
* **`local_K`**: An array of matrices storing the pre-computed Stiffness matrix for each element to optimize the time loop.
* **`u[3]`**: A matrix array storing three time generations: *Past*, *Present*, and *Future*.

## Plot
* Plotted using python with the file `plots2d.py`.
* To plot the error we use `error_plot_2d.py` if we use the initial conditions of which we have the analytical solution.


# Requirements
Polylib library: routines for polynomial integration, differentiation and interpolation. Used to calculate the GLL points.
https://www.nektar.info/2nd_edition/Polylib.html
Linear-Algebra-C: A linear algebra library for performing basic and advanced matrix/vector operations.
https://github.com/barrettotte/Linear-Algebra-C


