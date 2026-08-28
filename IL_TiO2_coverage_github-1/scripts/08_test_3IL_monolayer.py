from ase.io import read
from ase.geometry import get_distances
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il0  = read("input/BMP_FTFSI_trans.xyz")

Lx, Ly, Lz = slab.cell.lengths()
z_surface = slab.positions[:, 2].max()

def place_il(xcenter, ycenter, angle_z):
    il = il0.copy()
    il.translate(-il.get_center_of_mass())
    il.rotate(angle_z, "z", center="COM")
    com = il.get_center_of_mass()
    il.translate([xcenter - com[0], ycenter - com[1], 0.0])
    zmin = il.positions[:, 2].min()
    il.translate([0.0, 0.0, z_surface + 3.0 - zmin])
    return il

il1 = place_il(Lx / 6, Ly / 2, 90)
il2 = place_il(Lx / 2, Ly / 2, 270)
il3 = place_il(5 * Lx / 6, Ly / 2, 90)

system = slab + il1 + il2 + il3
system.set_cell(slab.cell)
system.set_pbc([True, True, True])

def minimum_distance(a, b):
    _, distances = get_distances(a.positions, b.positions, cell=system.cell, pbc=system.pbc)
    idx = np.unravel_index(np.argmin(distances), distances.shape)
    return distances[idx], idx

for label, a, b in [
    ("IL1 - IL2", il1, il2),
    ("IL2 - IL3", il2, il3),
    ("IL1 - IL3", il1, il3),
]:
    d, idx = minimum_distance(a, b)
    print(label, ":", d, "Angstrom", a[idx[0]].symbol, "-", b[idx[1]].symbol)

for n, il in enumerate([il1, il2, il3], start=1):
    d, idx = minimum_distance(slab, il)
    print(f"TiO2 - IL{n}:", d, "Angstrom", slab[idx[0]].symbol, "-", il[idx[1]].symbol)
