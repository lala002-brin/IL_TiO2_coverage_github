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

il1 = place_il(Lx * 0.25, Ly * 0.50, 90)
il2 = place_il(Lx * 0.75, Ly * 0.50, 270)

system = slab + il1 + il2
system.set_cell(slab.cell)
system.set_pbc([True, True, True])

def mindist(a, b):
    _, d = get_distances(a.positions, b.positions, cell=system.cell, pbc=system.pbc)
    idx = np.unravel_index(np.argmin(d), d.shape)
    return d[idx], idx

d12, p12 = mindist(il1, il2)
d1, p1 = mindist(slab, il1)
d2, p2 = mindist(slab, il2)

print("Total atoms:", len(system))
print("IL1 - IL2 :", d12, "Angstrom", il1[p12[0]].symbol, "-", il2[p12[1]].symbol)
print("TiO2 - IL1:", d1, "Angstrom", slab[p1[0]].symbol, "-", il1[p1[1]].symbol)
print("TiO2 - IL2:", d2, "Angstrom", slab[p2[0]].symbol, "-", il2[p2[1]].symbol)
