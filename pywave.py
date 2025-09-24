import numpy as np
import os

nx, ny, nz = 30, 30, 30   # keep it smaller for speed
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
z = np.linspace(-1, 1, nz)
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

c = 1.0
times = np.linspace(0, 1.0, 20)  # 20 frames

for idx, t in enumerate(times):
    r = np.sqrt(X**2 + Y**2 + Z**2)
    wavefield = np.sin(10 * (r - c*t)) / (r + 1e-6)

    output_dir = "timestep"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filename = f"timestep/wave3d_{idx:03d}.txt"
    with open(filename, "w") as f:
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    f.write(f"{x[i]} {y[j]} {z[k]} {wavefield[i,j,k]}\n")
