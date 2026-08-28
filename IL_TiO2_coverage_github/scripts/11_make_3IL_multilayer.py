from ase.io import read, write
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

top_vacuum = 15.0
z_top = system.positions[:, 2].max()
cell = system.cell.array.copy()

if cell[2, 2] - z_top < top_vacuum:
    cell[2, 2] = z_top + top_vacuum
    system.set_cell(cell, scale_atoms=False)

def minimum_distance(a, b):
    _, distances = get_distances(a.positions, b.positions, cell=system.cell, pbc=system.pbc)
    idx = np.unravel_index(np.argmin(distances), distances.shape)
    return distances[idx], idx

d12, _ = minimum_distance(il1, il2)
d13, _ = minimum_distance(il1, il3)
d23, _ = minimum_distance(il2, il3)

print("Total atoms :", len(system))
print("Formula     :", system.get_chemical_formula())
print("Cell        :", system.cell.lengths())
print("IL1 - IL2  :", d12, "Angstrom")
print("IL1 - IL3  :", d13, "Angstrom")
print("IL2 - IL3  :", d23, "Angstrom")
print("Top vacuum :", system.cell.lengths()[2] - system.positions[:, 2].max(), "Angstrom")

write("outputs/TiO2_3BMP_FTFSI_multilayer.cif", system)
write("outputs/TiO2_3BMP_FTFSI_multilayer.xyz", system)
