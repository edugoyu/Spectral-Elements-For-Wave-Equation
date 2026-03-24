import numpy as np
import matplotlib.pyplot as plt

# --- 1. Parameters (Must match your C code) ---
Lx, Ly = 10.0, 6.0
time_int = 3.0
time_points = 15000
save_freq = 50  # Frequency of saveStepToFile calls

dt = time_int / (time_points - 1)

def analytical_solution(x, y, t):
    """Analytical solution for a vibrating membrane with u=0 at boundaries."""
    kx = np.pi / Lx
    ky = np.pi / Ly
    omega = np.sqrt(kx**2 + ky**2)
    return np.sin(kx * x) * np.sin(ky * y) * np.cos(omega * t)

# --- 2. Load Simulation Data ---
# Row 0: X coords, Row 1: Y coords, Rows 2+: U at each saved timestep
data = np.genfromtxt('simulation_results_2d.csv', delimiter=',')

x_nodes = data[0, :]
y_nodes = data[1, :]
u_history = data[2:, :]

# --- 3. Compute Error Metrics ---
mean_errors = []
max_errors = []
times = []

print(f"Analyzing {len(u_history)} timesteps...")

for i in range(len(u_history)):
    # Calculate physical time for the current row
    # The first row of u_history corresponds to the first save (t_i = 0 or 1)
    # Adjust logic if your loop starts saving at t_i = 0
    t_curr = i * save_freq * dt
    times.append(t_curr)
    
    u_num = u_history[i, :]
    u_ana = analytical_solution(x_nodes, y_nodes, t_curr)
    
    # Calculate absolute error at each node
    abs_err = np.abs(u_num - u_ana)
    
    mean_errors.append(np.mean(abs_err))
    max_errors.append(np.max(abs_err))

# --- 4. 1D Plotting ---
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ... (Keep your existing data loading and error calculation logic here) ...

# --- 4. 1D Plotting ---
plt.figure(figsize=(10, 6))

# Plot both metrics
# Added "log-scale" to the labels to include it in the legend
plt.plot(times, mean_errors, label='Mean Absolute Error (log-scale)', color='blue', linewidth=1.5)
plt.plot(times, max_errors, label='Max Absolute Error (log-scale)', color='red', linestyle='--', linewidth=1.5)

# Formatting the Y-axis
plt.yscale('log')

# Custom formatter to show only the exponent (e.g., 10^-6 becomes -6)
def log_formatter(x, pos):
    return f"{int(np.log10(x))}"

plt.gca().yaxis.set_major_formatter(FuncFormatter(log_formatter))

plt.xlabel('Time (seconds)')
plt.ylabel('Absolute Error (Exponent Base 10)')
plt.title('Global Error over Time (SEM 2D Wave Equation)')

# Grid and Legend
plt.grid(True, which="major", ls="-", alpha=0.5)
plt.legend()

# Save and Show
plt.tight_layout()
plt.savefig('error_analysis_1d.png')
plt.show()

print(f"Final Max Error: {max_errors[-1]:.2e}")
print(f"Final Mean Error: {mean_errors[-1]:.2e}")