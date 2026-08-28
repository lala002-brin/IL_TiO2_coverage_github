from ase.io import read

print("STEP 1 - mulai")

slab = read("input/TiO2_110_3x5.cif")
print("STEP 2 - slab terbaca:", len(slab))

il = read("input/BMP_FTFSI_trans.xyz")
print("STEP 3 - IL terbaca:", len(il))

il.translate(-il.get_center_of_mass())
print("STEP 4 - IL dicenter")

Lx, Ly, Lz = slab.cell.lengths()
com = il.get_center_of_mass()

il.translate([Lx / 2 - com[0], Ly / 2 - com[1], 0.0])
print("STEP 5 - posisi XY selesai")

z_surface = slab.positions[:, 2].max()
z_il_min = il.positions[:, 2].min()

il.translate([0.0, 0.0, z_surface + 3.0 - z_il_min])
print("STEP 6 - posisi Z selesai")

system = slab + il
print("STEP 7 - sistem digabung")
print("Total atoms:", len(system))

system.set_cell(slab.cell)
system.set_pbc([True, True, True])
print("STEP 8 - selesai")
