import numpy as np
import struct
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#read the bunch of binary data
def read_frame(filename):
    with open(filename, 'rb') as f:
        nx = struct.unpack('Q', f.read(8))[0]
        ny = struct.unpack('Q', f.read(8))[0]
        t  = struct.unpack('d', f.read(8))[0]
        data = np.frombuffer(f.read(), dtype=np.float64).reshape(nx, ny, 4)
    return t, data

#primitive values from state vector
def compute_primitives(data, gamma=1.4):
    rho = data[:,:,0]
    u   = data[:,:,1] / rho
    v   = data[:,:,2] / rho
    E   = data[:,:,3]
    p   = (gamma - 1.0) * (E - 0.5 * rho * (u**2 + v**2))
    return rho, u, v, p

#load frames
outdir = input("Output directory: ").strip()
gamma  = float(input("Gamma: ").strip())

files = sorted(glob.glob(os.path.join(outdir, "frame_*.bin")))
if not files:
    raise FileNotFoundError(f"No frames found in '{outdir}'")

print(f"Found {len(files)} frames, loading...")
frames = [read_frame(f) for f in files]
times  = [fr[0] for fr in frames]
data   = [fr[1] for fr in frames]

#global color limits
def lims(arr_list):
    return min(a.min() for a in arr_list), max(a.max() for a in arr_list)

prims    = [compute_primitives(d, gamma) for d in data]
rho_lim  = lims([p[0] for p in prims])
u_lim    = lims([p[1] for p in prims])
v_lim    = lims([p[2] for p in prims])
p_lim    = lims([p[3] for p in prims])

#figure
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
fig.subplots_adjust(hspace=0.4, wspace=0.35)

titles = [r"Density  $\rho$", r"$x$-velocity  $u$", r"$y$-velocity  $v$", r"Pressure  $p$"]
clims  = [rho_lim, u_lim, v_lim, p_lim]
cmaps  = ["inferno", "RdBu_r", "RdBu_r", "viridis"]

rho0, u0, v0, p0 = prims[0]
fields0 = [rho0, u0, v0, p0]

meshes = []
for ax, field, title, clim, cmap in zip(axes.flat, fields0, titles, clims, cmaps):
    mesh = ax.pcolormesh(field.T, cmap=cmap, vmin=clim[0], vmax=clim[1], shading='auto')
    fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    meshes.append(mesh)

suptitle = fig.suptitle(f"t = {times[0]:.4f}", fontsize=13, fontweight='bold')

#update
def update(n):
    rho, u, v, p = prims[n]
    for mesh, field in zip(meshes, [rho, u, v, p]):
        mesh.set_array(field.T.ravel())
    suptitle.set_text(f"t = {times[n]:.4f}")
    return meshes + [suptitle]

ani = animation.FuncAnimation(fig, update, frames=len(data), interval=50, blit=False)
plt.show()

#save
if input("Save as mp4? (y/n): ").strip().lower() == 'y':
    out = os.path.join(outdir, "animation.mp4")
    ani.save(out, writer='ffmpeg', fps=200, dpi=300)
    print(f"Saved to {out}")

#check at t=0.2
rho_final = prims[-1][0]  #shape: (nx, ny)
u_final = prims[-1][1]
p_final = prims[-1][3]
j_mid = rho_final.shape[1] // 2
rho_slice = rho_final[:, j_mid]
u_slice = u_final[:, j_mid]
p_slice = p_final[:, j_mid]
x = np.linspace(0, 1, len(rho_slice))

plt.figure(figsize=(8, 4))
plt.plot(x, rho_slice, 'b-', lw=1.5, label=r'$\rho(x,y=0.5)$')
plt.plot(x, u_slice, 'r--', lw=1.5, label=r'$u(x,y=0.5)$')
plt.plot(x, p_slice, 'g-.', lw=1.5, label=r'$p(x,y=0.5)$')
plt.xlabel("x")
plt.ylabel(f"Primitive vals HLLC t={times[-1]:.4f}")
plt.title("shock tube snippet primitive vals at t~0.2")
plt.legend()
plt.grid(ls=":",alpha=0.5)
plt.tight_layout()
plt.savefig("./Sod_shock_tube/sod_shock_tube_1d.png", bbox_inches="tight", dpi=300)
plt.show()
