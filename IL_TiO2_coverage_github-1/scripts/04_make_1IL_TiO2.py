from ase.io import read, write
from ase.geometry import get_distances
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il   = read("input/BMP_FTFSI_trans.xyz")

il.translate(-il.get_center_of_mass())
Lx, Ly, Lz = slab.cell.lengths()
com = il.get_center_of_mass()
il.translate([Lx / 2 - com[0], Ly / 2 - com[1], 0.0])

z_surface = slab.positions[:, 2].max()
z_il_min  = il.positions[:, 2].min()
il.translate([0.0, 0.0, z_surface + 3.0 - z_il_min])

system = slab + il
system.set_cell(slab.cell)
system.set_pbc([True, True, True])

top_vacuum = 15.0
z_top = system.positions[:, 2].max()
cell = system.cell.array.copy()

if cell[2, 2] - z_top < top_vacuum:
    cell[2, 2] = z_top + top_vacuum
    system.set_cell(cell, scale_atoms=False)

_, distances = get_distances(
    slab.positions, il.positions, cell=system.cell, pbc=system.pbc
)
dmin = distances.min()

print("Total atoms       :", len(system))
print("Formula           :", system.get_chemical_formula())
print("Cell              :", system.cell.lengths())
print("Minimum IL-TiO2   :", dmin, "Angstrom")
print("Top vacuum        :", system.cell.lengths()[2] - system.positions[:, 2].max())

write("outputs/TiO2_1BMP_FTFSI.cif", system)
write("outputs/TiO2_1BMP_FTFSI.xyz", system)

print("Generated:")
print("outputs/TiO2_1BMP_FTFSI.cif")
print("outputs/TiO2_1BMP_FTFSI.xyz")
