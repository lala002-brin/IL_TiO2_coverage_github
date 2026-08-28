from ase.io import read
from ase.geometry import get_distances
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il0  = read("input/BMP_FTFSI_trans.xyz")

Lx, Ly, Lz = slab.cell.lengths()
z_surface = slab.positions[:, 2].max()

rng = np.random.default_rng(12345)
angles = np.arange(0, 360, 30)
ntrial = 5000

base = il0.copy()
base.translate(-base.get_center_of_mass())

def make_il(xcenter, ycenter, angle):
    il = base.copy()
    il.rotate(angle, "z", center="COM")
    com = il.get_center_of_mass()
    il.translate([xcenter - com[0], ycenter - com[1], 0.0])
    zmin = il.positions[:, 2].min()
    il.translate([0.0, 0.0, z_surface + 3.0 - zmin])
    return il

def minimum_distance_pbc(a, b):
    diff = a.positions[:, None, :] - b.positions[None, :, :]
    diff[:, :, 0] -= np.rint(diff[:, :, 0] / Lx) * Lx
    diff[:, :, 1] -= np.rint(diff[:, :, 1] / Ly) * Ly
    dist = np.linalg.norm(diff, axis=2)
    index = np.unravel_index(np.argmin(dist), dist.shape)
    return dist[index], index

best_score = -1.0
best_data = None

print("Searching", ntrial, "configurations ...")

for trial in range(ntrial):
    ils = []
    parameters = []

    for n in range(3):
        x = rng.uniform(0.0, Lx)
        y = rng.uniform(0.0, Ly)
        angle = int(rng.choice(angles))
        il = make_il(x, y, angle)
        ils.append(il)
        parameters.append((x, y, angle))

    d12, _ = minimum_distance_pbc(ils[0], ils[1])
    d13, _ = minimum_distance_pbc(ils[0], ils[2])
    d23, _ = minimum_distance_pbc(ils[1], ils[2])

    score = min(d12, d13, d23)

    if score > best_score:
        best_score = score
        best_data = {
            "parameters": parameters,
            "distances": (d12, d13, d23),
            "ils": [x.copy() for x in ils],
        }

print("\n=== BEST CONFIGURATION ===")
for n, (x, y, angle) in enumerate(best_data["parameters"], start=1):
    print(f"IL{n}: x = {x:.6f} y = {y:.6f} angle = {angle}")

d12, d13, d23 = best_data["distances"]
print("\n=== IL - IL DISTANCES ===")
print("IL1 - IL2 :", d12, "Angstrom")
print("IL1 - IL3 :", d13, "Angstrom")
print("IL2 - IL3 :", d23, "Angstrom")
print("Worst minimum distance:", best_score, "Angstrom")

print("\n=== TiO2 - IL ===")
for n, il in enumerate(best_data["ils"], start=1):
    _, distances = get_distances(slab.positions, il.positions, cell=slab.cell, pbc=slab.pbc)
    idx = np.unravel_index(np.argmin(distances), distances.shape)
    print(
        f"TiO2 - IL{n}:",
        distances[idx], "Angstrom",
        slab[idx[0]].symbol, "-", il[idx[1]].symbol
    )

print("\nBelum ada file CIF yang dibuat.")
