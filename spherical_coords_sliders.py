import numpy as np
import tkinter as tk
from tkinter import ttk

# --- Spherical <-> Cartesian Conversion Functions ---
def vec2deg(vec):
    vec = np.asarray(vec, dtype=np.float64)
    n = vec.size
    if n < 2:
        raise ValueError("Vector must have at least 2 dimensions")
    norm = np.linalg.norm(vec)
    if not np.isclose(norm, 1.0):
        vec = vec / norm
    angles = []
    for i in range(n - 1):
        if i == n - 2:
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

# --- Tkinter GUI with Sliders ---
class SphericalSliderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("n-Sphere Spherical Coordinate Sliders")
        self.max_dim = 16
        self.min_dim = 2
        self.cur_dim = 3
        self.suppress_callback = False  # Prevent feedback loops

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
        # Remove old sliders
        for s in self.vec_sliders:
            s.destroy()
        for s in self.ang_sliders:
            s.destroy()
        self.vec_sliders.clear()
        self.ang_sliders.clear()
        # Default: unit vector along x
        self.vec_vals = np.zeros(n)
        self.vec_vals[0] = 1.0
        self.ang_vals = vec2deg(self.vec_vals)
        # Vector sliders
        for i in range(n):
            s = tk.Scale(self.vec_frame, from_=-1.0, to=1.0, resolution=0.01, orient=tk.VERTICAL, length=200,
                         label=f'x{i+1}', command=lambda val, idx=i: self.on_vec_slider(idx, val))
            s.set(self.vec_vals[i])
            s.pack(side=tk.LEFT, padx=2)
            self.vec_sliders.append(s)
        # Angle sliders
        for i in range(n-1):
            # For last angle (phi), use -180 to 180; others use 0 to 180
            if i == n-2:
                from_, to_ = -180, 180
            else:
                from_, to_ = 0, 180
            s = tk.Scale(self.ang_frame, from_=from_, to=to_, resolution=0.1, orient=tk.VERTICAL, length=200,
                         label=f'θ{i+1}' if i < n-2 else 'φ', command=lambda val, idx=i: self.on_ang_slider(idx, val))
            s.set(self.ang_vals[i])
            s.pack(side=tk.LEFT, padx=2)
            self.ang_sliders.append(s)
        self.update_info()

    def on_dim_change(self, val):
        n = int(float(val))
        self.cur_dim = n
        self.init_sliders(n)

    def on_vec_slider(self, idx, val):
        if self.suppress_callback:
            return
        self.suppress_callback = True
        # Update vector value
        self.vec_vals[idx] = float(val)
        norm = np.linalg.norm(self.vec_vals)
        if norm > 1.0:
            self.vec_vals = self.vec_vals / norm
        # Update all angle sliders
        self.ang_vals = vec2deg(self.vec_vals)
        for i, s in enumerate(self.ang_sliders):
            s.configure(command=None)
            s.set(self.ang_vals[i])
            s.configure(command=lambda val, idx=i: self.on_ang_slider(idx, val))
        self.suppress_callback = False
        self.update_info()

    def on_ang_slider(self, idx, val):
        if self.suppress_callback:
            return
        self.suppress_callback = True
        # Update only the angle value that was changed
        self.ang_vals[idx] = float(val)
        # Update vector
        self.vec_vals = deg2vec(self.ang_vals)
        # Normalize (should already be unit, but for safety)
        norm = np.linalg.norm(self.vec_vals)
        if not np.isclose(norm, 1.0):
            self.vec_vals = self.vec_vals / norm
        # Update all vector sliders
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
