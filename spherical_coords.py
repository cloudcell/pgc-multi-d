import numpy as np

# --- Spherical <-> Cartesian Conversion Functions ---
def vec2deg(vec):
    """
    Convert a unit vector on an n-sphere to spherical coordinates in degrees.
    Returns an array of n-1 angles in degrees.
    For n=2: returns [phi], where phi = atan2(y, x)
    For n=3: returns [theta, phi], where theta = arccos(z), phi = atan2(y, x)
    For n>3: returns [theta_1, ..., theta_{n-2}, phi]
    """
    vec = np.asarray(vec, dtype=np.float64)
    n = vec.size
    if n < 2:
        raise ValueError("Vector must have at least 2 dimensions")
    norm = np.linalg.norm(vec)
    if not np.isclose(norm, 1.0):
        raise ValueError("Input vector must be normalized (unit length)")
    angles = []
    for i in range(n - 1):
        if i == n - 2:
            # Last angle: phi = atan2(x2, x1)
            phi = np.arctan2(vec[1], vec[0])
            angles.append(np.degrees(phi))
        else:
            denom = np.linalg.norm(vec[i:])
            if denom == 0:
                angle = 0.0
            else:
                angle = np.arccos(vec[i] / denom)
            angles.append(np.degrees(angle))
    return np.array(angles)

def deg2vec(angles):
    """
    Convert spherical coordinates in degrees to a unit vector on an n-sphere.
    Input: angles (n-1,)
    Output: vector (n,)
    For n=2: angles=[phi]
    For n=3: angles=[theta, phi]
    For n>3: angles=[theta_1, ..., theta_{n-2}, phi]
    """
    angles = np.asarray(angles, dtype=np.float64)
    n = angles.size + 1
    if n < 2:
        raise ValueError("Angles array must have at least 1 element (for 2D)")
    theta = np.radians(angles)
    vec = np.zeros(n)
    for i in range(n):
        prod = 1.0
        for j in range(i):
            if j == n - 2:
                prod *= np.sin(theta[-1])
            else:
                prod *= np.sin(theta[j])
        if i == n - 1:
            vec[i] = prod
        elif i == n - 2:
            vec[i] = prod * np.cos(theta[i - 1])
        else:
            vec[i] = prod * np.cos(theta[i])
    return vec

# --- Example Vectors and Usage ---
EXAMPLES = {
    "2D": np.array([np.cos(np.radians(30)), np.sin(np.radians(30))]),
    "3D": np.array([1/np.sqrt(3)]*3),
    "16D": np.ones(16) / np.sqrt(16),
}

import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, needed for 3D plots

def gui():
    root = tk.Tk()
    root.title("n-Sphere Spherical Coordinate Converter")

    # --- Frames ---
    input_frame = ttk.Frame(root, padding=10)
    input_frame.pack(side=tk.TOP, fill=tk.X)
    output_frame = ttk.Frame(root, padding=10)
    output_frame.pack(side=tk.TOP, fill=tk.X)
    chart_frame = ttk.Frame(root, padding=10)
    chart_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # --- Input ---
    dim_var = tk.StringVar(value="3")
    tk.Label(input_frame, text="Dimension (n):").grid(row=0, column=0)
    dim_entry = ttk.Entry(input_frame, textvariable=dim_var, width=4)
    dim_entry.grid(row=0, column=1)

    tk.Label(input_frame, text="Vector (comma-separated):").grid(row=1, column=0)
    vec_var = tk.StringVar(value="1,0,0")
    vec_entry = ttk.Entry(input_frame, textvariable=vec_var, width=40)
    vec_entry.grid(row=1, column=1, columnspan=3, sticky="ew")

    tk.Label(input_frame, text="Angles [deg] (comma-separated):").grid(row=2, column=0)
    ang_var = tk.StringVar(value="0,0")
    ang_entry = ttk.Entry(input_frame, textvariable=ang_var, width=40)
    ang_entry.grid(row=2, column=1, columnspan=3, sticky="ew")

    # --- Output ---
    out_var = tk.StringVar()
    tk.Label(output_frame, text="Output:").pack(side=tk.LEFT)
    out_label = tk.Label(output_frame, textvariable=out_var, fg="blue")
    out_label.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # --- Chart ---
    fig = plt.Figure(figsize=(6, 10))
    # We'll use 3 rows: main chart, vector bar, angles bar
    canvas = FigureCanvasTkAgg(fig, master=chart_frame)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def update_chart(vec, angles, n):
        fig.clf()
        # Main spherical chart (top)
        gs = fig.add_gridspec(3, 1, height_ratios=[2, 1, 1])
        if n == 2:
            ax_main = fig.add_subplot(gs[0, 0])
            ax_main.set_aspect('equal')
            ax_main.set_xlim(-1.1, 1.1)
            ax_main.set_ylim(-1.1, 1.1)
            circle = plt.Circle((0, 0), 1, color='gray', fill=False)
            ax_main.add_patch(circle)
            ax_main.arrow(0, 0, vec[0], vec[1], head_width=0.05, head_length=0.1, color='red', length_includes_head=True)
            ax_main.plot([0, vec[0]], [0, vec[1]], 'r-', label="Vector")
            ax_main.set_title(f"2D: φ={angles[0]:.1f}°")
            ax_main.legend()
        elif n == 3:
            ax_main = fig.add_subplot(gs[0, 0], projection='3d')
            u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
            x = np.cos(u) * np.sin(v)
            y = np.sin(u) * np.sin(v)
            z = np.cos(v)
            ax_main.plot_surface(x, y, z, color='lightgray', alpha=0.2, linewidth=0)
            ax_main.quiver(0, 0, 0, vec[0], vec[1], vec[2], color='r', linewidth=2, arrow_length_ratio=0.1)
            ax_main.set_xlim([-1, 1])
            ax_main.set_ylim([-1, 1])
            ax_main.set_zlim([-1, 1])
            ax_main.set_title(f"3D: θ={angles[0]:.1f}°, φ={angles[1]:.1f}°")
        else:
            ax_main = fig.add_subplot(gs[0, 0])
            ax_main.text(0.5, 0.5, f"No chart for n={n}", ha='center', va='center', fontsize=16)
            ax_main.axis('off')

        # Bar chart for vector (middle)
        ax_vec = fig.add_subplot(gs[1, 0])
        ax_vec.bar(range(1, len(vec)+1), vec, color='skyblue')
        ax_vec.set_title("Vector Components")
        ax_vec.set_xlabel("Index")
        ax_vec.set_ylabel("Value")
        ax_vec.set_xticks(range(1, len(vec)+1))
        ax_vec.set_ylim(-1.1, 1.1)

        # Bar chart for angles (bottom)
        ax_ang = fig.add_subplot(gs[2, 0])
        if len(angles) > 0:
            ax_ang.bar(range(1, len(angles)+1), angles, color='orange')
            ax_ang.set_xticks(range(1, len(angles)+1))
        ax_ang.set_title("Spherical Angles (deg)")
        ax_ang.set_xlabel("Angle Index")
        ax_ang.set_ylabel("Degrees")
        ax_ang.set_ylim(-180, 180)
        fig.tight_layout()
        canvas.draw()

    def convert_vec2deg():
        try:
            n = int(dim_var.get())
            vec = np.fromstring(vec_var.get(), sep=',')
            if vec.size != n:
                raise ValueError(f"Vector must have {n} elements.")
            norm = np.linalg.norm(vec)
            if not np.isclose(norm, 1.0):
                vec = vec / norm
            angles = vec2deg(vec)
            ang_var.set(','.join(f"{a:.6f}" for a in angles))
            out_var.set(f"Angles [deg]: {angles}")
            update_chart(vec, angles, n)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def convert_deg2vec():
        try:
            n = int(dim_var.get())
            angles = np.fromstring(ang_var.get(), sep=',')
            if angles.size != n - 1:
                raise ValueError(f"Angles must have {n-1} elements.")
            vec = deg2vec(angles)
            norm = np.linalg.norm(vec)
            if not np.isclose(norm, 1.0):
                vec = vec / norm
            vec_var.set(','.join(f"{v:.6f}" for v in vec))
            out_var.set(f"Vector: {vec}")
            update_chart(vec, angles, n)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    ttk.Button(input_frame, text="Vec → Deg", command=convert_vec2deg).grid(row=3, column=1, pady=5)
    ttk.Button(input_frame, text="Deg → Vec", command=convert_deg2vec).grid(row=3, column=2, pady=5)

    # Initial chart
    try:
        n = int(dim_var.get())
        vec = np.fromstring(vec_var.get(), sep=',')
        if vec.size == n:
            norm = np.linalg.norm(vec)
            if not np.isclose(norm, 1.0):
                vec = vec / norm
            angles = vec2deg(vec)
        else:
            angles = np.zeros(n-1)
        update_chart(vec, angles, n)
    except:
        pass

    root.mainloop()

if __name__ == "__main__":
    gui()
