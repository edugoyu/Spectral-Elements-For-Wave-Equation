import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 1. Load the data
df = pd.read_csv('simulation_results.csv')
x = df.columns.astype(float).values
data = df.values

# 2. Parameters
L = 10.0           
time_int = 3.0     
num_steps = len(data)
time_delta = time_int / (num_steps - 1)

def true_solution(x, t):
    return np.sin(np.pi * x / L) * np.cos(t)

# 3. Setup the Figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [2, 1]})
plt.subplots_adjust(hspace=0.3)

# Top Plot: Wave Animation
ax1.set_xlim(0, L)
ax1.set_ylim(-1.2, 1.2)
ax1.grid(True, linestyle=':', alpha=0.6)
line_num, = ax1.plot([], [], 'o', color='teal', markersize=3, label='SEM (Numerical)')
line_true, = ax1.plot([], [], '--', color='orange', alpha=0.6, label='Analytical')
time_text = ax1.text(0.02, 0.90, '', transform=ax1.transAxes, fontweight='bold')
ax1.set_ylabel('Displacement (u)')
ax1.legend(loc='upper right')

# Bottom Plot: Logarithmic Error
ax2.set_yscale('log')
ax2.set_xlim(0, time_int)
ax2.set_ylim(1e-9, 1.0) 
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Error (Log Scale)')

# Two lines for the error plot
mean_error_line, = ax2.plot([], [], color='teal', lw=1.5, label='Mean Error (MAE)')
max_error_line, = ax2.plot([], [], color='crimson', lw=1.5, linestyle='-', label=r'Max Error ($L_\infty$)')
ax2.grid(True, which="both", linestyle=':', alpha=0.4)
ax2.legend(loc='lower right', fontsize='small')

# History containers
time_history = []
mean_history = []
max_history = []

def init():
    line_num.set_data([], [])
    line_true.set_data([], [])
    mean_error_line.set_data([], [])
    max_error_line.set_data([], [])
    time_text.set_text('')
    return line_num, line_true, mean_error_line, max_error_line, time_text

def update(frame):
    t_curr = frame * time_delta
    
    # 1. Update Waveforms
    line_num.set_data(x, data[frame])
    x_fine = np.linspace(0, L, 500)
    line_true.set_data(x_fine, true_solution(x_fine, t_curr))
    
    # 2. Calculate Errors
    analytical_at_nodes = true_solution(x, t_curr)
    errors = np.abs(data[frame] - analytical_at_nodes)
    
    mean_err = np.mean(errors)
    max_err = np.max(errors)
    
    # Floor at 1e-10 for the log plot
    plot_mean = max(mean_err, 1e-10)
    plot_max = max(max_err, 1e-10)
    
    time_history.append(t_curr)
    mean_history.append(plot_mean)
    max_history.append(plot_max)
    
    # 3. Update Lines
    mean_error_line.set_data(time_history, mean_history)
    max_error_line.set_data(time_history, max_history)
    
    time_text.set_text(f'Time: {t_curr:.3f}s | Max Err: {max_err:.2e}')
    return line_num, line_true, mean_error_line, max_error_line, time_text

skip_frames = 50 
ani = FuncAnimation(fig, update, frames=range(0, num_steps, skip_frames),
                    init_func=init, blit=False, interval=25, repeat=False)

plt.show()
fig.savefig('final_simulation_frame.png', dpi=300, bbox_inches='tight')