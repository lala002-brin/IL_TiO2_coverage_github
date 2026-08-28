from ase.io import read
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

def minimum_periodic_self_distance(il):
    positions = il.positions
    translations = [
        [ Lx,  0, 0], [-Lx,  0, 0],
        [0,  Ly, 0],  [0, -Ly, 0],
        [ Lx,  Ly, 0], [ Lx, -Ly, 0],
        [-Lx,  Ly, 0], [-Lx, -Ly, 0],
    ]

    dmin = 999.0
    best_translation = None
    best_pair = None

    for T in translations:
        image = positions + np.array(T)
        diff = positions[:, None, :] - image[None, :, :]
        distances = np.linalg.norm(diff, axis=2)
        index = np.unravel_index(np.argmin(distances), distances.shape)
        value = distances[index]

        if value < dmin:
            dmin = value
            best_translation = T
            best_pair = index

    return dmin, best_translation, best_pair

for name, il in [("IL1", il1), ("IL2", il2)]:
    d, t, p = minimum_periodic_self_distance(il)
    print(f"=== {name} periodic image ===")
    print("Minimum distance :", d, "Angstrom")
    print("Translation      :", t)
    print("Atoms            :", il[p[0]].symbol, "-", il[p[1]].symbol)
