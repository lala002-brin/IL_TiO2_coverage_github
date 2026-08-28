from ase.io import read, write
from ase.constraints import FixAtoms
from pathlib import Path
import numpy as np

structure_file = "outputs/TiO2_1BMP_FTFSI.cif"
pseudo_dir = Path("/mgpfs/home/lala002/pseudo")

pseudos = {
    "Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "O" : "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "C" : "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H" : "H.pbe-kjpaw_psl.1.0.0.UPF",
    "N" : "N.pbe-n-kjpaw_psl.1.0.0.UPF",
    "F" : "F.pbe-n-kjpaw_psl.1.0.0.UPF",
    "S" : "S.pbe-n-kjpaw_psl.1.0.0.UPF",
}

print("=== CHECK PSEUDOPOTENTIALS ===")
missing = []

for element, filename in pseudos.items():
    path = pseudo_dir / filename
    if path.exists():
        print(f"{element:2s} : OK      {filename}")
    else:
        print(f"{element:2s} : MISSING {filename}")
        missing.append((element, filename))

if missing:
    print("\nERROR: pseudopotential belum lengkap.")
    print("Edit mapping 'pseudos' sesuai file yang tersedia di pseudo_dir.")
    raise SystemExit(1)

atoms = read(structure_file)

print("\n=== STRUCTURE ===")
print("Atoms   :", len(atoms))
print("Formula :", atoms.get_chemical_formula())
print("Cell    :", atoms.cell.lengths())

symbols = np.array(atoms.get_chemical_symbols())
z = atoms.positions[:, 2]

ti_mask = symbols == "Ti"
z_ti_min = z[ti_mask].min()
z_ti_max = z[ti_mask].max()
z_fix_cutoff = 0.5 * (z_ti_min + z_ti_max)

fixed_indices = [
    i for i, atom in enumerate(atoms)
    if atom.symbol in ["Ti", "O"] and atom.position[2] <= z_fix_cutoff
]

atoms.set_constraint(FixAtoms(indices=fixed_indices))

print("\n=== CONSTRAINT ===")
print("Ti z min       :", z_ti_min)
print("Ti z max       :", z_ti_max)
print("Fixed cutoff   :", z_fix_cutoff)
print("Fixed atoms    :", len(fixed_indices))
print("Mobile atoms   :", len(atoms) - len(fixed_indices))

input_data = {
    "control": {
        "calculation": "relax",
        "restart_mode": "from_scratch",
        "prefix": "TiO2_1BMP_FTFSI",
        "pseudo_dir": str(pseudo_dir),
        "outdir": "./tmp",
        "tstress": True,
        "tprnfor": True,
        "etot_conv_thr": 1.0e-4,
        "forc_conv_thr": 1.0e-3,
    },
    "system": {
        "ibrav": 0,
        "ecutwfc": 60,
        "ecutrho": 480,
        "nspin": 1,
        "tot_charge": 0.0,
        "occupations": "smearing",
        "smearing": "mv",
        "degauss": 0.01,
        "vdw_corr": "grimme-d3",
        "nosym": True,
    },
    "electrons": {
        "conv_thr": 1.0e-6,
        "electron_maxstep": 200,
        "mixing_beta": 0.2,
        "diagonalization": "david",
    },
    "ions": {
        "ion_dynamics": "bfgs",
    },
}

output_file = "qe/TiO2_1BMP_FTFSI_relax.in"

write(
    output_file,
    atoms,
    format="espresso-in",
    pseudopotentials=pseudos,
    input_data=input_data,
    kpts=None,
)

print("\nGenerated:", output_file)
print("Gamma point")
print("ecutwfc = 60 Ry")
print("ecutrho = 480 Ry")
print("PBE-D3")
