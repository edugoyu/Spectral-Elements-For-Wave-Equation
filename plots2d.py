import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import matplotlib.animation as animation

# --- 1. Load Data ---
data = np.genfromtxt('simulation_results_2d.csv', delimiter=',')
x = data[0, :]
y = data[1, :]
u_history = data[2:, :]

# --- 2. Interpolation Grid ---
xi = np.linspace(x.min(), x.max(), 100)
yi = np.linspace(y.min(), y.max(), 100)
X, Y = np.meshgrid(xi, yi)

# --- 3. Save Snapshots in One Figure ---
snapshot_indices = [0, len(u_history)//4, len(u_history)//2, len(u_history)-1]
fig_snap = plt.figure(figsize=(16, 12))
fig_snap.suptitle("Wave Propagation Snapshots", fontsize=20)

for i, idx in enumerate(snapshot_indices):
    ax_snap = fig_snap.add_subplot(2, 2, i+1, projection='3d')
    z_step = u_history[idx, :]
    Z = griddata((x, y), z_step, (X, Y), method='cubic')
    
    surf = ax_snap.plot_surface(X, Y, Z, cmap=cm.coolwarm, 
                               antialiased=False, rstride=2, cstride=2)
    
    ax_snap.set_title(f"Time Step: {idx}")
    ax_snap.set_zlim(-1.1, 1.1)
    # Adjusting view angle for better clarity
    ax_snap.view_init(elev=30, azim=45)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("simulation_snapshots.png", dpi=300)
print("Snapshots saved to simulation_snapshots.png")

# --- 4. Animation Logic (as before) ---
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

def animate(i):
    ax.clear()
    z_step = u_history[i, :]
    Z = griddata((x, y), z_step, (X, Y), method='cubic')
    
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, 
                           antialiased=False, rstride=2, cstride=2)
    
    ax.set_zlim(-1.1, 1.1) 
    ax.set_title(f"Wave Propagation - Frame {i}")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    return surf,

ani = animation.FuncAnimation(fig, animate, frames=range(0, len(u_history), 2), interval=50)

plt.show()