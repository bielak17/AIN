# Python script: 3D viewer with button to switch between animated functions
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.widgets import Button

files = {
    "rosenbrock": "rosenbrock_2d.csv",
    "salomon": "salomon_2d.csv",
    "whitley": "whitley_2d.csv"
}


def load_csv_or_demo(kind, fname):
    if os.path.exists(fname):
        try:
            df = pd.read_csv(fname)
            if not set(["x","y","z"]).issubset(df.columns):
                raise ValueError(f"{fname} doesn't contain x,y,z columns")
            return df
        except Exception as e:
            print(f"Failed to load {fname}: {e} -> creating demo data instead.")
            return make_demo_data(fname, kind)
    else:
        print(f"{fname} not found - creating demo data file.")
        return make_demo_data(fname, kind)

datasets = {}
for k, fname in files.items():
    datasets[k] = load_csv_or_demo(k, fname)

def try_grid(df):
    xs = np.sort(df['x'].unique())
    ys = np.sort(df['y'].unique())
    if len(xs) * len(ys) == len(df):
        try:
            X = xs
            Y = ys
            Z = df.pivot(index='y', columns='x', values='z').values
            return np.meshgrid(xs, ys), Z
        except Exception:
            return None
    return None

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(bottom=0.15)

keys = list(datasets.keys())
idx = 0
current_plot = {"type": None, "obj": None}

def plot_dataset(key):
    global current_plot
    df = datasets[key]
    ax.clear()
    ax.set_box_aspect([1,1,0.6])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(f"{key} - {files[key]}")
    grid_result = try_grid(df)
    if grid_result is not None:
        (Xg, Yg), Zg = grid_result
        surf = ax.plot_surface(Xg, Yg, Zg, rstride=1, cstride=1, linewidth=0, antialiased=True)
        current_plot = {"type":"surface", "obj": surf, "X": Xg, "Y": Yg, "Z": Zg}
    else:
        sc = ax.scatter(df['x'], df['y'], df['z'], s=6)
        current_plot = {"type":"scatter", "obj": sc, "data": df}

plot_dataset(keys[idx])

angle = 0
pulse_phase = 0.0


ax_button = plt.axes([0.4, 0.03, 0.2, 0.06])
btn = Button(ax_button, 'Next function')

def on_next(event):
    global idx
    idx = (idx + 1) % len(keys)
    plot_dataset(keys[idx])
    fig.canvas.draw_idle()

btn.on_clicked(on_next)

print("Interactive viewer ready.")
print(" - Button 'Next function' switches between datasets: " + ", ".join(keys))
print(" - The plot auto-rotates. Close the window to end the program if running as a script.")
plt.show()
