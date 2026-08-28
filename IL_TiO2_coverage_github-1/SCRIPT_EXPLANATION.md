# Penjelasan Masing-Masing Script

Dokumen ini menjelaskan fungsi setiap script di repository **IL_TiO2_coverage_github** yang digunakan untuk membangun model ionic liquid **[BMP][FTFSI]** pada permukaan **rutile TiO2(110) 3×5**, mengecek overlap, membandingkan beberapa jumlah ion pair, membangun multilayer, dan menyiapkan input awal Quantum ESPRESSO.

---


## 0. `scripts/00_make_TiO2_110_3x5.py`

### Fungsi

Script ini membuat **surface rutile TiO2(110)** dari struktur bulk `TiO2.cif`, lalu memperbesar surface menjadi supercell **3×5**.

Input:

```text
input/TiO2.cif
```

Bulk dibaca dengan:

```python
bulk = read("input/TiO2.cif")
```

Surface dibangun menggunakan ASE:

```python
slab = surface(
    bulk,
    indices=(1, 1, 0),
    layers=4,
    vacuum=15.0,
    periodic=True
)
```

Parameter yang digunakan:

```text
Miller index = (1 1 0)
layers       = 4
vacuum       = 15 Å
```

Setelah surface dibuat, slab dipusatkan pada arah z:

```python
slab.center(
    vacuum=15.0,
    axis=2
)
```

Tujuannya adalah memberikan ruang vacuum di atas dan di bawah slab sehingga periodic image sepanjang z tidak berinteraksi terlalu kuat.

### Membuat supercell 3×5

Surface dasar kemudian diperbesar:

```python
slab_3x5 = slab.repeat((3, 5, 1))
```

Artinya:

```text
3 kali sepanjang arah cell X
5 kali sepanjang arah cell Y
1 kali sepanjang arah Z
```

Pada bulk rutile yang digunakan dalam workflow ini, surface dasar menghasilkan kira-kira:

```text
TiO2(110) 1×1
Formula : O16Ti8
Atoms   : 24
Cell    : 6.5807 × 2.9692 × 41.8754 Å
```

Setelah `repeat((3,5,1))`:

```text
TiO2(110) 3×5
Formula : O240Ti120
Atoms   : 360
Cell    : 19.7422 × 14.8460 × 41.8754 Å
```

### Output

```text
outputs/TiO2_110_1x1.cif
outputs/TiO2_110_1x1.xyz
input/TiO2_110_3x5.cif
outputs/TiO2_110_3x5.xyz
```

File:

```text
input/TiO2_110_3x5.cif
```

kemudian menjadi input untuk semua script penempatan ionic liquid.

### Tujuan

Urutan workflow lengkap sekarang menjadi:

```text
bulk TiO2.cif
      ↓
buat TiO2(110)
      ↓
4-layer slab + 15 Å vacuum
      ↓
repeat 3×5
      ↓
TiO2_110_3x5.cif
      ↓
tambahkan [BMP][FTFSI]
```

---

## 1. `scripts/01_check_IL_TiO2.py`

### Fungsi

Script ini digunakan untuk **memeriksa struktur awal** sebelum TiO2 dan ionic liquid digabungkan.

Input:

```text
input/TiO2_110_3x5.cif
input/BMP_FTFSI_trans.xyz
```

Script mengecek:

- jumlah atom,
- formula kimia,
- ukuran cell TiO2,
- ukuran molekul ionic liquid pada arah X, Y, dan Z.

Bagian utama:

```python
slab = read("input/TiO2_110_3x5.cif")
il   = read("input/BMP_FTFSI_trans.xyz")
```

ASE membaca kedua struktur.

Ukuran ionic liquid dihitung dengan:

```python
size = np.ptp(il.positions, axis=0)
```

Secara matematis:

\[
\Delta x = x_{\max} - x_{\min}
\]

dan sama untuk Y dan Z.

### Hasil pada sistem ini

```text
TiO2:
360 atom
O240Ti120

[BMP][FTFSI]:
42 atom
C10H20F4N2O4S2

Ukuran IL:
X ≈ 12.11 Å
Y ≈ 5.28 Å
Z ≈ 6.59 Å
```

### Tujuan

Memastikan input terbaca dengan benar dan memperkirakan apakah ukuran permukaan cukup untuk menempatkan ionic liquid.

---

## 2. `scripts/02_test_add_1IL.py`

### Fungsi

Script diagnostik sederhana untuk memastikan bahwa:

```text
TiO2 + 1 [BMP][FTFSI]
```

dapat digabungkan tanpa masalah.

Urutan proses:

```text
read TiO2
    ↓
read ionic liquid
    ↓
center ionic liquid
    ↓
center pada X-Y slab
    ↓
letakkan 3 Å di atas surface
    ↓
gabungkan
```

Ionic liquid dipusatkan dengan:

```python
il.translate(-il.get_center_of_mass())
```

Kemudian dipindahkan ke tengah cell:

```python
il.translate([
    Lx / 2 - com[0],
    Ly / 2 - com[1],
    0.0
])
```

Posisi vertikal ditentukan dari atom paling atas TiO2 dan atom paling bawah IL:

```python
z_surface = slab.positions[:, 2].max()
z_il_min = il.positions[:, 2].min()
```

Kemudian IL diletakkan sekitar 3 Å dari permukaan:

```python
il.translate([
    0.0,
    0.0,
    z_surface + 3.0 - z_il_min
])
```

### Jumlah atom

\[
360 + 42 = 402
\]

### Tujuan

Script ini belum menghitung jarak minimum dan belum menulis file. Fungsinya hanya memastikan proses placement dasar berhasil.

---

## 3. `scripts/03_test_distance.py`

### Fungsi

Menghitung **jarak minimum antara TiO2 dan ionic liquid** setelah placement awal.

Digunakan:

```python
from ase.geometry import get_distances
```

Kemudian:

```python
_, distances = get_distances(
    slab.positions,
    il.positions,
    cell=system.cell,
    pbc=system.pbc
)
```

Jarak minimum:

```python
dmin = distances.min()
```

### Hasil pada konfigurasi N=1

```text
Minimum TiO2–IL ≈ 3.16 Å
Closest pair = O – F
```

### Tujuan

Melakukan **overlap checking** sebelum struktur disimpan atau digunakan dalam perhitungan DFT.

Jarak sangat pendek seperti:

```text
0.8 Å
1.0 Å
1.3 Å
```

biasanya menandakan atom saling menembus dan placement harus diperbaiki.

---

# Model N = 1

## 4. `scripts/04_make_1IL_TiO2.py`

### Fungsi

Membuat struktur final awal:

```text
TiO2 + 1 [BMP][FTFSI]
```

Jumlah atom:

\[
360 + 42 = 402
\]

Selain placement, script juga memastikan adanya vacuum yang cukup di atas adsorbat.

Atom tertinggi dicari dengan:

```python
z_top = system.positions[:, 2].max()
```

Jika vacuum terlalu kecil:

```python
if cell[2, 2] - z_top < top_vacuum:
    cell[2, 2] = z_top + top_vacuum
```

Cell diperpanjang tanpa mengubah posisi atom:

```python
system.set_cell(cell, scale_atoms=False)
```

### Output

```text
outputs/TiO2_1BMP_FTFSI.cif
outputs/TiO2_1BMP_FTFSI.xyz
```

### Hasil pada workflow ini

```text
Total atom = 402
Minimum TiO2–IL ≈ 3.16 Å
Top vacuum = 15 Å
```

### Tujuan

Menghasilkan model N=1 yang siap diperiksa secara visual dan selanjutnya digunakan untuk generator input QE.

---

# Model N = 2

## 5. `scripts/05_test_2IL_TiO2.py`

### Fungsi

Menguji penempatan dua ion pair pada permukaan yang sama.

Posisi yang dipakai:

```python
il1 = place_il(Lx * 0.25, Ly * 0.50, 90)
il2 = place_il(Lx * 0.75, Ly * 0.50, 270)
```

Artinya:

```text
IL1 → sekitar 1/4 panjang cell X
IL2 → sekitar 3/4 panjang cell X
```

Rotasi:

```text
IL1 = 90°
IL2 = 270°
```

Tujuannya agar bentuk molekul yang memanjang dapat tersusun lebih baik.

Script menghitung:

```text
IL1–IL2
TiO2–IL1
TiO2–IL2
```

### Hasil pada workflow ini

```text
IL1–IL2 ≈ 4.91 Å
TiO2–IL1 ≈ 3.12 Å
TiO2–IL2 ≈ 3.12 Å
```

### Tujuan

Memastikan dua ion pair tidak saling overlap dan tetap berada pada jarak yang masuk akal terhadap permukaan.

---

## 6. `scripts/06_check_periodic_IL.py`

### Fungsi

Memeriksa jarak ionic liquid terhadap **periodic image dirinya sendiri**.

Dalam simulasi periodik, satu molekul memiliki salinan pada cell tetangga:

```text
| cell -1 | cell 0 | cell +1 |
```

Script membuat translasi:

```python
translations = [
    [ Lx,  0, 0],
    [-Lx,  0, 0],
    [0,  Ly, 0],
    [0, -Ly, 0],
    ...
]
```

Kemudian dihitung jarak atom pada IL terhadap copy periodiknya.

### Hasil pada workflow ini

```text
Minimum periodic distance ≈ 3.32 Å
Closest pair = F – H
```

### Tujuan

Mencegah struktur yang tampak aman di dalam satu cell tetapi sebenarnya bertabrakan dengan periodic image.

---

## 7. `scripts/07_make_2IL_TiO2.py`

### Fungsi

Membuat model final awal:

```text
TiO2 + 2 [BMP][FTFSI]
```

Jumlah atom:

\[
360 + 2(42) = 444
\]

Script:

- memakai konfigurasi N=2 yang telah diuji,
- menghitung jarak IL–IL,
- menghitung jarak IL–TiO2,
- memperbaiki top vacuum menjadi 15 Å,
- menulis CIF dan XYZ.

### Output

```text
outputs/TiO2_2BMP_FTFSI.cif
outputs/TiO2_2BMP_FTFSI.xyz
```

### Formula

```text
C20H40F8N4O248S4Ti120
```

### Interpretasi

N=2 merepresentasikan **coverage yang lebih tinggi** dibanding N=1 pada luas permukaan TiO2 yang sama.

---

# Pengujian N = 3 Monolayer

## 8. `scripts/08_test_3IL_monolayer.py`

### Fungsi

Menguji apakah tiga ion pair dapat ditempatkan semuanya pada first layer.

Susunan awal:

```text
IL1       IL2       IL3

=======================
        TiO2
=======================
```

Posisi lateral:

```python
Lx / 6
Lx / 2
5 * Lx / 6
```

### Hasil pada workflow ini

```text
IL1–IL2 ≈ 1.87 Å
IL2–IL3 ≈ 2.35 Å
IL1–IL3 ≈ 2.33 Å
```

### Interpretasi

Jarak antar-IL terlalu pendek. Konfigurasi ini dianggap terlalu padat untuk digunakan sebagai struktur awal QE.

### Tujuan

Script ini adalah **uji kelayakan**, bukan generator struktur produksi.

---

## 9. `scripts/09_search_3IL_staggered.py`

### Fungsi

Mencari konfigurasi tiga IL yang lebih baik secara otomatis menggunakan pencarian posisi dan rotasi.

Random seed:

```python
rng = np.random.default_rng(12345)
```

Seed membuat hasil dapat direproduksi.

Jumlah percobaan:

```text
5000 konfigurasi
```

Rotasi yang dicoba:

```text
0°
30°
60°
...
330°
```

Untuk setiap trial dihitung:

```text
d12 = IL1–IL2
d13 = IL1–IL3
d23 = IL2–IL3
```

Skor:

```python
score = min(d12, d13, d23)
```

Artinya, algoritma memilih konfigurasi yang memaksimalkan jarak terburuk:

\[
S = \max [\min(d_{12}, d_{13}, d_{23})]
\]

### Hasil pada workflow ini

```text
Worst minimum distance ≈ 1.77 Å
```

### Interpretasi

Walaupun 5000 konfigurasi dicoba, tiga ion pair masih terlalu padat untuk monolayer pada cell TiO2 3×5 ini.

### Tujuan

Memberikan dasar bahwa masalah bukan hanya posisi awal, tetapi keterbatasan ruang lateral pada cell yang digunakan.

---

# Model N = 3 Multilayer

## 10. `scripts/10_test_3IL_multilayer.py`

### Fungsi

Menguji konfigurasi:

```text
             IL3
        second layer

      IL1       IL2
        first layer

========================
        TiO2
========================
```

IL1 dan IL2 memakai konfigurasi N=2 yang sebelumnya berhasil.

Posisi IL3 ditentukan berdasarkan atom tertinggi first layer:

```python
z_first_top = max(
    il1.positions[:, 2].max(),
    il2.positions[:, 2].max()
)
```

Kemudian bagian bawah IL3 ditempatkan sekitar 3 Å di atas first layer.

### Hasil pada workflow ini

```text
IL1–IL2 ≈ 4.91 Å
IL1–IL3 ≈ 3.56 Å
IL2–IL3 ≈ 5.30 Å
```

Jarak TiO2–IL3:

```text
≈ 10.89 Å
```

Nilai tersebut wajar karena IL3 berada pada second layer dan tidak langsung mengadsorpsi ke TiO2.

### Tujuan

Menguji apakah jumlah tiga ion pair lebih masuk akal sebagai awal pembentukan multilayer daripada dipaksakan dalam satu monolayer.

---

## 11. `scripts/11_make_3IL_multilayer.py`

### Fungsi

Membuat model N=3 final awal:

```text
TiO2
+
2 IL first layer
+
1 IL second layer
```

Jumlah atom:

\[
360 + 3(42) = 486
\]

Karena IL3 berada cukup tinggi, cell Z asli tidak lagi cukup.

Pada workflow ini:

```text
Old cell Z ≈ 41.875 Å
Highest z ≈ 46.046 Å
```

Cell diperpanjang menjadi:

\[
46.046 + 15 \approx 61.046\ \text{Å}
\]

### Output

```text
outputs/TiO2_3BMP_FTFSI_multilayer.cif
outputs/TiO2_3BMP_FTFSI_multilayer.xyz
```

### Hasil

```text
Total atom = 486
Cell Z ≈ 61.046 Å
Top vacuum = 15 Å
```

### Interpretasi

Model ini harus dibaca sebagai:

```text
N=1 → first-layer coverage rendah
N=2 → first-layer coverage lebih tinggi
N=3 → onset multilayer / film formation
```

N=3 **bukan** tiga ion pair dalam satu monolayer.

---

# Quantum ESPRESSO

## 12. `qe/make_qe_N1.py`

### Fungsi

Mengubah:

```text
outputs/TiO2_1BMP_FTFSI.cif
```

menjadi input awal Quantum ESPRESSO untuk geometry relaxation.

---

### A. Pseudopotential

Pseudopotential dipetakan melalui:

```python
pseudos = {
    "Ti": "...UPF",
    "O" : "...UPF",
    "C" : "...UPF",
    "H" : "...UPF",
    "N" : "...UPF",
    "F" : "...UPF",
    "S" : "...UPF",
}
```

Script memeriksa keberadaan masing-masing file.

Jika salah satu tidak ditemukan:

```text
F : MISSING
```

atau:

```text
S : MISSING
```

script berhenti.

### Tujuan

Mencegah pembuatan input QE dengan nama pseudopotential yang tidak tersedia.

---

### B. Constraint pada slab

Bagian bawah TiO2 dibuat fixed.

Script mencari range koordinat z atom Ti:

```python
z_ti_min
z_ti_max
```

Kemudian cutoff awal:

```python
z_fix_cutoff = 0.5 * (z_ti_min + z_ti_max)
```

Atom Ti/O di bawah cutoff dimasukkan ke:

```python
FixAtoms
```

Ionic liquid tetap mobile.

### Tujuan

Membiarkan bagian atas surface dan adsorbat relax, sambil menjaga bagian bawah slab mendekati struktur bulk.

---

### C. Jenis calculation

Digunakan:

```text
calculation = 'relax'
```

Artinya:

```text
atom bergerak
cell tetap
```

Cell tidak dioptimasi seperti pada `vc-relax`.

---

### D. Cutoff

Baseline awal:

```text
ecutwfc = 60 Ry
ecutrho = 480 Ry
```

Nilai ini adalah parameter awal. Untuk perhitungan final tetap perlu convergence test.

---

### E. Dispersion

Digunakan:

```python
"vdw_corr": "grimme-d3"
```

Dispersion penting untuk interaksi ionic liquid, ion pair, dan permukaan.

---

### F. Spin dan charge

Digunakan:

```text
nspin = 1
tot_charge = 0
```

Satu complete ion pair terdiri dari:

\[
[BMP]^+ + [FTFSI]^-
\]

sehingga total charge cell:

\[
q = 0
\]

---

### G. K-point

Generator memakai:

```python
kpts=None
```

yang menghasilkan Gamma-point-only untuk large supercell.

Gamma point cocok sebagai initial test, tetapi tetap perlu convergence check untuk calculation final.

---

### H. Output

```text
qe/TiO2_1BMP_FTFSI_relax.in
```

Input ini kemudian dapat dijalankan menggunakan `pw.x`.

---

# `run_structure_workflow.sh`

### Fungsi

Shell script untuk menjalankan hampir seluruh workflow struktur secara berurutan.

Contoh:

```bash
python3 scripts/01_check_IL_TiO2.py
python3 scripts/03_test_distance.py
python3 scripts/04_make_1IL_TiO2.py
...
python3 scripts/11_make_3IL_multilayer.py
```

Dijalankan dengan:

```bash
./run_structure_workflow.sh
```

### Catatan

Untuk debugging dan belajar, disarankan tetap menjalankan script satu per satu. Script ini terutama berguna setelah seluruh workflow sudah tervalidasi.

---

# Ringkasan Workflow

```text
bulk TiO2.cif
        │
        ▼
00_make_TiO2_110_3x5.py
        │
        ▼
TiO2_110_3x5.cif
        +
BMP_FTFSI_trans.xyz
        │
        ▼
01_check_IL_TiO2.py
        │
        ▼
placement + distance checking
        │
        ├───────────────┐
        ▼               ▼
      N=1             N=2
   402 atoms         444 atoms
        │               │
        ▼               ▼
  make_1IL         make_2IL
                        │
                        ▼
                 test 3 IL monolayer
                        │
                    terlalu padat
                        │
                        ▼
                 search staggered
                        │
                    masih padat
                        │
                        ▼
                3 IL multilayer
                   486 atoms
                        │
                        ▼
                 QE relaxation
```

Prinsip utama repository ini adalah:

```text
BUILD
  ↓
CHECK
  ↓
ACCEPT / REJECT
  ↓
GENERATE
  ↓
QUANTUM ESPRESSO
```

Struktur awal tidak langsung digunakan untuk DFT tanpa pengecekan overlap, periodic image, jarak antarmolekul, dan kecukupan vacuum terlebih dahulu.
