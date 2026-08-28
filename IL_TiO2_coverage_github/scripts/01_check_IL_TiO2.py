from ase.io import read
import numpy as np

slab = read("input/TiO2_110_3x5.cif")
il   = read("input/BMP_FTFSI_trans.xyz")

print("=== TiO2 ===")
print("Atoms   :", len(slab))
print("Formula :", slab.get_chemical_formula())
print("Cell    :", slab.cell.lengths())

print("\n=== Ionic Liquid ===")
print("Atoms   :", len(il))
print("Formula :", il.get_chemical_formula())

size = np.ptp(il.positions, axis=0)
print("Size X  :", size[0], "Angstrom")
print("Size Y  :", size[1], "Angstrom")
print("Size Z  :", size[2], "Angstrom")
