import numpy as np
import tkinter as tk
from tkinter import ttk

# --- Spherical <-> Cartesian Conversion Functions ---
def vec2deg(vec):
    vec = np.asarray(vec, dtype=np.float64)
    n = vec.size
    if n < 2:
        raise ValueError("Vector must have at least 2 dimensions")
    if np.allclose(vec, 0):
        return np.zeros(n-1)
    norm = np.linalg.norm(vec)
    vec_unit = vec / norm
    angles = []
    for i in range(n - 1):
        if i == 0:
            # Polar angle θ from z-axis
            theta = np.arccos(vec_unit[-1])
            angles.append(np.degrees(theta))
        elif i < n - 2:
            # For higher dimensions: generalized θ
            denom = np.linalg.norm(vec_unit[-(i+1):])
            if denom == 0:
                angle = 0.0
            else:
                angle = np.arccos(vec_unit[-(i+2)] / denom)
            angles.append(np.degrees(angle))
        else:
            # Last angle is azimuth φ in xy-plane
            phi = np.arctan2(vec_unit[1], vec_unit[0])
            angles.append(np.degrees(phi))
    return np.array(angles)

def deg2vec(angles):
    angles = np.asarray(angles, dtype=np.float64)
    n = angles.size + 1
    theta = np.radians(angles)
    vec = np.zeros(n)
    for i in range(n):
        prod = 1.0
        for j in range(i):
            prod *= np.sin(theta[j])
        if i < n - 1:
            vec[i] = prod * np.cos(theta[i])
        else:
            vec[i] = prod
    return vec

# --- SphericalSliderApp ---
class SphericalSliderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("n-Sphere Spherical Coordinate Sliders")
        self.max_dim = 16
        self.min_dim = 2
        self.cur_dim = 3
        self.suppress_callback = False

        # Top frame for dimension slider
        top_frame = ttk.Frame(root, padding=5)
        top_frame.pack(side=tk.TOP, fill=tk.X)
        tk.Label(top_frame, text="Dimension (n):").pack(side=tk.LEFT)
        self.dim_slider = tk.Scale(top_frame, from_=self.min_dim, to=self.max_dim, orient=tk.HORIZONTAL, command=self.on_dim_change)
        self.dim_slider.set(self.cur_dim)
        self.dim_slider.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # Main frame for sliders
        main_frame = ttk.Frame(root, padding=5)
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Vector sliders (left)
        self.vec_frame = ttk.LabelFrame(main_frame, text="Vector Components", padding=5)
        self.vec_frame.pack(side=tk.LEFT, fill=tk.Y)
        # Angle sliders (right)
        self.ang_frame = ttk.LabelFrame(main_frame, text="Spherical Angles [deg]", padding=5)
        self.ang_frame.pack(side=tk.RIGHT, fill=tk.Y)

        # Info label (bottom)
        self.info_var = tk.StringVar()
        tk.Label(root, textvariable=self.info_var, fg="blue").pack(side=tk.BOTTOM, fill=tk.X)

        # Initialize sliders
        self.vec_sliders = []
        self.ang_sliders = []
        self.init_sliders(self.cur_dim)

    def init_sliders(self, n):
        for s in self.vec_sliders + self.ang_sliders:
            s.destroy()
        self.vec_sliders.clear()
        self.ang_sliders.clear()

        self.vec_vals = np.zeros(n)
        self.vec_vals[0] = 1.0  # x1 locked
        self.ang_vals = vec2deg(self.vec_vals)

        for i in range(n):
            s = tk.Scale(self.vec_frame, from_=-1.0, to=1.0, resolution=0.01, orient=tk.VERTICAL, length=200,
                         label=f'x{i+1}', command=lambda val, idx=i: self.on_vec_slider(idx, val))
            s.set(self.vec_vals[i])
            s.pack(side=tk.LEFT, padx=2)
            self.vec_sliders.append(s)

        for i in range(n-1):
            # Polar angles (θ): [0,180], Azimuthal angle (φ): [-180,180]
            if i < n - 2:
                from_, to_ = (0, 180)
                label = f'θ{i+1}'
            else:
                from_, to_ = (-180, 180)
                label = 'φ'
            s = tk.Scale(self.ang_frame, from_=from_, to=to_, resolution=0.1, orient=tk.VERTICAL, length=200,
                         label=label, command=lambda val, idx=i: self.on_ang_slider(idx, val))
            s.set(self.ang_vals[i])
            s.pack(side=tk.LEFT, padx=2)
            self.ang_sliders.append(s)

        self.update_info()

    def on_dim_change(self, val):
        self.cur_dim = int(float(val))
        self.init_sliders(self.cur_dim)

    def on_vec_slider(self, idx, val):
        if self.suppress_callback:
            return
        self.suppress_callback = True
        self.vec_vals[idx] = float(val)
        # Do NOT update angles—vector sliders are raw values
        for i, s in enumerate(self.vec_sliders):
            s.configure(command=lambda val, idx=i: self.on_vec_slider(idx, val))
        self.suppress_callback = False
        self.update_info()

    def on_ang_slider(self, idx, val):
        if self.suppress_callback:
            return
        self.suppress_callback = True
        self.ang_vals[idx] = float(val)
        # Recompute all vector components from angles
        self.vec_vals = deg2vec(self.ang_vals)
        for i, s in enumerate(self.vec_sliders):
            s.configure(command=None)
            s.set(self.vec_vals[i])
            s.configure(command=lambda val, idx=i: self.on_vec_slider(idx, val))
        self.suppress_callback = False
        self.update_info()

    def update_info(self):
        self.info_var.set(f"Vector: {np.round(self.vec_vals, 4)}  |  Angles [deg]: {np.round(self.ang_vals, 2)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = SphericalSliderApp(root)
    root.mainloop()
