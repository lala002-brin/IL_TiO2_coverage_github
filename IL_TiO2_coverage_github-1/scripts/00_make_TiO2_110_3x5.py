from ase.io import read, write
from ase.build import surface

# ============================================================
# INPUT BULK TiO2
# ============================================================

bulk_file = "input/TiO2.cif"

bulk = read(bulk_file)

print("=== BULK TiO2 ===")
print("Atoms   :", len(bulk))
print("Formula :", bulk.get_chemical_formula())
print("Cell    :", bulk.cell.lengths())

# ============================================================
# BUILD RUTILE TiO2(110) SLAB
# ============================================================
#
# Miller index : (1 1 0)
# Layers       : 4
# Vacuum       : 15 Angstrom
#
# periodic=True keeps periodic boundary conditions in the slab.
# ============================================================

slab = surface(
    bulk,
    indices=(1, 1, 0),
    layers=4,
    vacuum=15.0,
    periodic=True
)

# Ensure the slab is centered along z with 15 Å vacuum
slab.center(
    vacuum=15.0,
    axis=2
)

print("\n=== TiO2(110) 1x1 SLAB ===")
print("Atoms   :", len(slab))
print("Formula :", slab.get_chemical_formula())
print("Cell    :", slab.cell.lengths())
print(
    "z range :",
    slab.positions[:, 2].min(),
    "-",
    slab.positions[:, 2].max()
)

# Save the base slab
write(
    "outputs/TiO2_110_1x1.cif",
    slab
)

write(
    "outputs/TiO2_110_1x1.xyz",
    slab
)

# ============================================================
# MAKE 3 x 5 SURFACE SUPERCELL
# ============================================================

slab_3x5 = slab.repeat((3, 5, 1))

print("\n=== TiO2(110) 3x5 SLAB ===")
print("Atoms   :", len(slab_3x5))
print("Formula :", slab_3x5.get_chemical_formula())
print("Cell    :", slab_3x5.cell.lengths())
print(
    "z range :",
    slab_3x5.positions[:, 2].min(),
    "-",
    slab_3x5.positions[:, 2].max()
)

# This is the slab used by the ionic-liquid scripts
write(
    "input/TiO2_110_3x5.cif",
    slab_3x5
)

write(
    "outputs/TiO2_110_3x5.xyz",
    slab_3x5
)

print("\nGenerated:")
print("outputs/TiO2_110_1x1.cif")
print("outputs/TiO2_110_1x1.xyz")
print("input/TiO2_110_3x5.cif")
print("outputs/TiO2_110_3x5.xyz")
