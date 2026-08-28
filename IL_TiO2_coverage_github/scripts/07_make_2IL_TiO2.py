from ase.io import read, write
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

top_vacuum = 15.0
z_top = system.positions[:, 2].max()
cell = system.cell.array.copy()

if cell[2, 2] - z_top < top_vacuum:
    cell[2, 2] = z_top + top_vacuum
    system.set_cell(cell, scale_atoms=False)

_, d12 = get_distances(il1.positions, il2.positions, cell=system.cell, pbc=system.pbc)
_, d1  = get_distances(slab.positions, il1.positions, cell=system.cell, pbc=system.pbc)
_, d2  = get_distances(slab.positions, il2.positions, cell=system.cell, pbc=system.pbc)

print("Total atoms :", len(system))
print("Formula     :", system.get_chemical_formula())
print("Cell        :", system.cell.lengths())
print("IL1 - IL2  :", d12.min(), "Angstrom")
print("TiO2 - IL1 :", d1.min(), "Angstrom")
print("TiO2 - IL2 :", d2.min(), "Angstrom")
print("Top vacuum :", system.cell.lengths()[2] - system.positions[:, 2].max(), "Angstrom")

write("outputs/TiO2_2BMP_FTFSI.cif", system)
write("outputs/TiO2_2BMP_FTFSI.xyz", system)
