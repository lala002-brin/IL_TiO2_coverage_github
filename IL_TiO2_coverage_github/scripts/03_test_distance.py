from ase.io import read
from ase.geometry import get_distances
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il   = read("input/BMP_FTFSI_trans.xyz")

il.translate(-il.get_center_of_mass())
Lx, Ly, Lz = slab.cell.lengths()
com = il.get_center_of_mass()
il.translate([Lx / 2 - com[0], Ly / 2 - com[1], 0.0])

z_surface = slab.positions[:, 2].max()
z_il_min = il.positions[:, 2].min()
il.translate([0.0, 0.0, z_surface + 3.0 - z_il_min])

system = slab + il
system.set_cell(slab.cell)
system.set_pbc([True, True, True])

_, distances = get_distances(
    slab.positions, il.positions, cell=system.cell, pbc=system.pbc
)

dmin = distances.min()
index = np.unravel_index(np.argmin(distances), distances.shape)

print("Minimum distance :", dmin, "Angstrom")
print("Slab atom index  :", index[0])
print("IL atom index    :", index[1])
print("Slab atom        :", slab[index[0]].symbol)
print("IL atom          :", il[index[1]].symbol)
