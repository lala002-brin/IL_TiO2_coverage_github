# [BMP][FTFSI] / TiO2(110) Coverage Workflow with ASE + Quantum ESPRESSO

Repository ini merangkum workflow yang dibuat langkah demi langkah untuk membangun model
ionic liquid **[BMP][FTFSI]** pada slab **rutile TiO2(110) 3×5**, mengecek overlap,
menguji coverage 1–3 ion pair, membangun model multilayer, dan menyiapkan input awal
Quantum ESPRESSO.

## Sistem yang digunakan

Slab awal:

- Formula: `O240Ti120`
- Jumlah atom: `360`
- Cell sekitar: `19.742161 × 14.846015 × 41.875436 Å`

Ion pair:

- `[BMP][FTFSI]`
- Formula: `C10H20F4N2O4S2`
- Jumlah atom: `42`
- Ukuran awal sekitar: `12.114 × 5.279 × 6.585 Å`

Struktur yang dihasilkan:

| Model | Deskripsi | Jumlah atom |
|---|---|---:|
| N=1 | 1 ion pair pada permukaan | 402 |
| N=2 | 2 ion pair pada first layer | 444 |
| N=3 | 2 first-layer + 1 second-layer | 486 |

> Catatan: model `N=3` bukan tiga ion pair dalam satu monolayer. Percobaan pencarian
> konfigurasi lateral menunjukkan kontak terlalu pendek, sehingga ion pair ketiga
> ditempatkan sebagai awal second layer.

## Struktur folder

```text
IL_TiO2_coverage_github/
├── input/
│   ├── BMP_FTFSI_trans.xyz
│   └── TiO2_110_3x5.cif        # copy file slab Anda ke sini
├── scripts/
│   ├── 01_check_IL_TiO2.py
│   ├── 02_test_add_1IL.py
│   ├── 03_test_distance.py
│   ├── 04_make_1IL_TiO2.py
│   ├── 05_test_2IL_TiO2.py
│   ├── 06_check_periodic_IL.py
│   ├── 07_make_2IL_TiO2.py
│   ├── 08_test_3IL_monolayer.py
│   ├── 09_search_3IL_staggered.py
│   ├── 10_test_3IL_multilayer.py
│   └── 11_make_3IL_multilayer.py
├── qe/
│   └── make_qe_N1.py
├── outputs/
├── requirements.txt
└── .gitignore
```

## Persiapan

```bash
git clone <URL-REPOSITORY-ANDA>
cd IL_TiO2_coverage_github
```

Copy slab TiO2 Anda ke folder `input`:

```bash
cp /path/TiO2_110_3x5.cif input/
```

Install dependency jika diperlukan:

```bash
python3 -m pip install -r requirements.txt
```

Semua script dijalankan dari root repository.

## 1. Cek input

```bash
python3 scripts/01_check_IL_TiO2.py
```

Expected:

```text
TiO2 atoms : 360
IL atoms   : 42
```

## 2. Tes penempatan 1 ion pair

```bash
python3 scripts/02_test_add_1IL.py
python3 scripts/03_test_distance.py
```

Pada workflow awal didapat minimum contact TiO2–IL sekitar:

```text
3.16 Å
```

## 3. Buat model N=1

```bash
python3 scripts/04_make_1IL_TiO2.py
```

Output:

```text
outputs/TiO2_1BMP_FTFSI.cif
outputs/TiO2_1BMP_FTFSI.xyz
```

Target:

- total 402 atom
- top vacuum 15 Å

## 4. Tes dan buat model N=2

```bash
python3 scripts/05_test_2IL_TiO2.py
python3 scripts/06_check_periodic_IL.py
python3 scripts/07_make_2IL_TiO2.py
```

Hasil konfigurasi yang diperoleh pada workflow awal:

- IL1–IL2 ≈ `4.91 Å`
- TiO2–IL1 ≈ `3.12 Å`
- TiO2–IL2 ≈ `3.12 Å`
- closest periodic image ≈ `3.32 Å`

Output:

```text
outputs/TiO2_2BMP_FTFSI.cif
outputs/TiO2_2BMP_FTFSI.xyz
```

## 5. Tes 3 ion pair dalam satu layer

```bash
python3 scripts/08_test_3IL_monolayer.py
python3 scripts/09_search_3IL_staggered.py
```

Percobaan satu baris memberi kontak IL–IL sekitar `1.87–2.35 Å`.
Pencarian 5000 konfigurasi staggered juga masih menemukan worst minimum distance sekitar
`1.77 Å`.

Karena terlalu rapat, konfigurasi ini **tidak dipakai** untuk input QE.

## 6. Model N=3 multilayer

```bash
python3 scripts/10_test_3IL_multilayer.py
python3 scripts/11_make_3IL_multilayer.py
```

Konfigurasi:

- IL1 dan IL2 = first layer
- IL3 = second layer

Hasil workflow awal:

- IL1–IL2 ≈ `4.91 Å`
- IL1–IL3 ≈ `3.56 Å`
- IL2–IL3 ≈ `5.30 Å`
- top vacuum = `15 Å`
- total = `486 atom`
- cell Z ≈ `61.046 Å`

Output:

```text
outputs/TiO2_3BMP_FTFSI_multilayer.cif
outputs/TiO2_3BMP_FTFSI_multilayer.xyz
```

## 7. Quantum ESPRESSO — N=1

Edit pseudopotential names pada:

```text
qe/make_qe_N1.py
```

Lalu:

```bash
python3 qe/make_qe_N1.py
```

Generator menggunakan baseline awal:

- `calculation = 'relax'`
- PBE
- Grimme D3
- Gamma point
- `ecutwfc = 60 Ry`
- `ecutrho = 480 Ry`
- `nspin = 1`
- neutral system
- sebagian bawah TiO2 fixed

Cutoff, k-point, smearing, slab constraints, dan detail pseudopotential harus diuji/ditetapkan
sesuai convergence study sebelum calculation produksi.

## Pseudopotential

Default path pada generator QE:

```text
/mgpfs/home/lala002/pseudo
```

Nama pseudopotential Ti/O/C/H/N mengikuti library yang sebelumnya digunakan. Nama F dan S
harus dipastikan tersedia di mesin Anda. Script akan berhenti dengan pesan `MISSING` jika file
yang disebut tidak ditemukan.

## Catatan ilmiah

1. Selalu tambahkan **ion pair netral** `[BMP]+[FTFSI]-` untuk menjaga total charge cell = 0.
2. `N=1` dan `N=2` merepresentasikan perubahan coverage pada luas permukaan yang sama.
3. `N=3` pada repository ini merepresentasikan onset pembentukan multilayer/film, bukan
   monolayer coverage 3.
4. Struktur yang dibuat merupakan **initial configurations**, bukan minimum-energy structures.
   Geometry relaxation tetap diperlukan.
5. Untuk final production study, lakukan convergence test terhadap cutoff, k-point, slab
   thickness, vacuum, dan metodologi elektronik yang dipakai.

## Lisensi

Tambahkan lisensi sesuai kebutuhan proyek Anda sebelum repository dipublikasikan.
