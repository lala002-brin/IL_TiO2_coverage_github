from ase.io import read
from ase.geometry import get_distances
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il0  = read("input/BMP_FTFSI_trans.xyz")

Lx, Ly, Lz = slab.cell.lengths()
z_surface = slab.positions[:, 2].max()

def place_first_layer(xcenter, ycenter, angle):
    il = il0.copy()
    il.translate(-il.get_center_of_mass())
    il.rotate(angle, "z", center="COM")
    com = il.get_center_of_mass()
    il.translate([xcenter - com[0], ycenter - com[1], 0.0])
    zmin = il.positions[:, 2].min()
    il.translate([0.0, 0.0, z_surface + 3.0 - zmin])
    return il

il1 = place_first_layer(Lx * 0.25, Ly * 0.50, 90)
il2 = place_first_layer(Lx * 0.75, Ly * 0.50, 270)

il3 = il0.copy()
il3.translate(-il3.get_center_of_mass())
il3.rotate(0, "z", center="COM")
com = il3.get_center_of_mass()
il3.translate([Lx / 2 - com[0], Ly / 2 - com[1], 0.0])

z_first_top = max(il1.positions[:, 2].max(), il2.positions[:, 2].max())
z3_min = il3.positions[:, 2].min()
il3.translate([0.0, 0.0, z_first_top + 3.0 - z3_min])

system = slab + il1 + il2 + il3
system.set_cell(slab.cell)
system.set_pbc([True, True, True])

def minimum_distance(a, b):
    _, distances = get_distances(a.positions, b.positions, cell=system.cell, pbc=system.pbc)
    idx = np.unravel_index(np.argmin(distances), distances.shape)
    return distances[idx], idx

for label, a, b in [
    ("IL1 - IL2", il1, il2),
    ("IL1 - IL3", il1, il3),
    ("IL2 - IL3", il2, il3),
]:
    d, p = minimum_distance(a, b)
    print(label, ":", d, "Angstrom", a[p[0]].symbol, "-", b[p[1]].symbol)

print("\n=== TiO2 - IL ===")
for n, il in enumerate([il1, il2, il3], start=1):
    d, p = minimum_distance(slab, il)
    print(f"TiO2 - IL{n}:", d, "Angstrom", slab[p[0]].symbol, "-", il[p[1]].symbol)

print("\n=== Z POSITION ===")
print("Surface top :", z_surface)
print("IL1 range   :", il1.positions[:, 2].min(), "-", il1.positions[:, 2].max())
print("IL2 range   :", il2.positions[:, 2].min(), "-", il2.positions[:, 2].max())
print("IL3 range   :", il3.positions[:, 2].min(), "-", il3.positions[:, 2].max())
