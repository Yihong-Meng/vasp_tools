#!/usr/bin/env python3
"""
Band structure plotting from VASPKIT-style BAND.dat.
If BAND.dat doesn't exist, generate it from EIGENVAL first.

Usage:
  python3 plot_band.py <work_dir> <struct_name>
Example:
  python3 plot_band.py ./1/1_1 1_1
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ======== 1. Args ========
if len(sys.argv) >= 3:
    work_dir = sys.argv[1]
    struct_name = sys.argv[2]
elif len(sys.argv) == 2:
    work_dir = sys.argv[1]
    struct_name = "band"
else:
    work_dir = os.getcwd()
    struct_name = "band"
band_dat = os.path.join(work_dir, "BAND.dat")
eigenval = os.path.join(work_dir, "EIGENVAL")
outcar   = os.path.join(work_dir, "OUTCAR")
poscar   = os.path.join(work_dir, "POSCAR")
kpoints  = os.path.join(work_dir, "KPOINTS")

print(f"[Plot] {struct_name}")

# ======== 2. If BAND.dat doesn't exist, generate from EIGENVAL ========
if not os.path.exists(band_dat):
    print("  BAND.dat not found, generating from EIGENVAL ...")

    if not os.path.exists(eigenval):
        print(f"  ERROR: {eigenval} not found either!")
        sys.exit(1)

    # Read EIGENVAL
    with open(eigenval) as f:
        lines = f.readlines()

    # Find NKPTS and NBANDS by scanning forward.
    # Header format varies, but NKPTS/NBANDS line has 3 numbers:
    #   NELECT  NKPTS  NBANDS   (e.g. "336  7  224")
    # or 2 numbers: "NKPTS  NBANDS"
    nkpts, nbands = 0, 0
    header_line_idx = -1
    for i, line in enumerate(lines):
        parts = line.strip().split()
        if len(parts) == 3:
            try:
                nk = int(parts[1])
                nb = int(parts[2])
                if 1 <= nk <= 10000 and 1 <= nb <= 10000:
                    nkpts, nbands = nk, nb
                    header_line_idx = i
                    break
            except:
                continue
        elif len(parts) == 2:
            try:
                nk = int(parts[0])
                nb = int(parts[1])
                if 1 <= nk <= 10000 and 1 <= nb <= 10000:
                    nkpts, nbands = nk, nb
                    header_line_idx = i
                    break
            except:
                continue

    if nkpts == 0 or nbands == 0:
        print("  ERROR: Could not parse NKPTS/NBANDS from EIGENVAL header")
        print("  First 10 lines:")
        for k in range(min(10, len(lines))):
            print(f"    {k}: {lines[k].rstrip()}")
        sys.exit(1)

    print(f"  EIGENVAL: NKPTS={nkpts}, NBANDS={nbands}")

    # Read lattice from POSCAR for reciprocal space distance
    lattice = np.zeros((3, 3))
    with open(poscar) as f:
        plines = f.readlines()
    scale = float(plines[1].strip())
    for i in range(3):
        lattice[i] = [float(x) for x in plines[2+i].split()[:3]]
    lattice *= scale

    volume = np.dot(lattice[0], np.cross(lattice[1], lattice[2]))
    recip = np.zeros((3, 3))
    recip[0] = 2*np.pi * np.cross(lattice[1], lattice[2]) / volume
    recip[1] = 2*np.pi * np.cross(lattice[2], lattice[0]) / volume
    recip[2] = 2*np.pi * np.cross(lattice[0], lattice[1]) / volume

    def kcart(k):
        return np.dot(k, recip)

    # Parse k-points and bands (data starts after blank line following header)
    kpts_frac = []
    e_up_all = []
    e_down_all = []

    i = header_line_idx + 1
    # Skip blank lines until first k-point
    while i < len(lines) and not lines[i].strip():
        i += 1

    while i < len(lines) and len(kpts_frac) < nkpts:
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        parts = line.split()
        if len(parts) == 4:
            try:
                kx, ky, kz = float(parts[0]), float(parts[1]), float(parts[2])
            except:
                continue
            kpts_frac.append([kx, ky, kz])

            # Read bands until next k-point or EOF
            up_band, down_band = [], []
            while i < len(lines) and len(up_band) < nbands:
                bline = lines[i].strip()
                i += 1
                if not bline:
                    continue
                bp = bline.split()
                # Check if next k-point (4 numbers)
                if len(bp) == 4:
                    try:
                        _ = [float(x) for x in bp]
                        i -= 1  # un-consume, outer loop will handle
                        break
                    except:
                        pass
                if len(bp) >= 3:
                    try:
                        up_band.append(float(bp[1]))
                        down_band.append(float(bp[2]))
                    except:
                        continue

            e_up_all.append(up_band)
            e_down_all.append(down_band)

    nparsed_k = len(kpts_frac)
    nparsed_b = len(e_up_all[0]) if e_up_all else 0
    print(f"  Parsed {nparsed_k} k-points, {nparsed_b} bands")

    if nparsed_k < 2:
        print("  ERROR: Too few k-points!")
        sys.exit(1)

    # Truncate bands to min common count
    min_b = min(len(e) for e in e_up_all)
    e_up_arr = np.array([e[:min_b] for e in e_up_all])
    e_down_arr = np.array([e[:min_b] for e in e_down_all])
    kpts_arr = np.array(kpts_frac)

    # Calculate k-distances in 1/A
    kcart_arr = np.array([kcart(k) for k in kpts_arr])
    kdist = np.zeros(len(kcart_arr))
    for i in range(1, len(kcart_arr)):
        kdist[i] = kdist[i-1] + np.linalg.norm(kcart_arr[i] - kcart_arr[i-1])

    # Write BAND.dat (VASPKIT format, all bands forward)
    nk_actual = nparsed_k
    nb_actual = min_b
    with open(band_dat, 'w') as f:
        f.write("#K-Path(1/A)         Spin-Up(eV)   Spin-down(eV)\n")
        f.write(f"# NKPTS & NBANDS:   {nk_actual} {nb_actual}\n")
        for b in range(nb_actual):
            f.write(f"\n# Band-Index    {b+1}\n")
            for k_idx in range(nk_actual):
                f.write(f"   {kdist[k_idx]:.8f}    {e_up_arr[k_idx, b]:.6f}    {e_down_arr[k_idx, b]:.6f}\n")

    print(f"  Generated BAND.dat")

# ======== 3. Read BAND.dat ========
with open(band_dat) as f:
    lines = f.readlines()

# Parse NKPTS and NBANDS
nkpts, nbands = 0, 0
for line in lines[:5]:
    if "NKPTS" in line or "NBANDS" in line:
        parts = line.split(':')[1].strip().split()
        nkpts = int(parts[0])
        nbands = int(parts[1])
        break

print(f"  BAND.dat: NKPTS={nkpts}, NBANDS={nbands}")

# Read Fermi energy
efermi = 0.0
fermi_file = os.path.join(work_dir, "FERMI_ENERGY.in")
if os.path.exists(fermi_file):
    with open(fermi_file) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith('#') and not s.startswith('!'):
                try:
                    efermi = float(s.split()[0])
                except:
                    pass
                break
elif os.path.exists(outcar):
    with open(outcar) as f:
        for line in f:
            if "E-fermi" in line:
                try:
                    efermi = float(line.strip().split()[2])
                except:
                    pass
                break

print(f"  E_F = {efermi:.4f} eV (shifted to 0 in plot)")

# Parse band data
all_kdist = []
all_up = []
all_down = []

i = 0
while i < len(lines) and ("#K-Path" in lines[i] or "NKPTS" in lines[i] or not lines[i].strip()):
    i += 1

band_count = 0
while i < len(lines) and band_count < nbands:
    line = lines[i].strip()
    i += 1
    if "Band-Index" not in line:
        continue

    band_count += 1
    kdist_local = []
    up = []
    down = []
    read = 0
    while read < nkpts and i < len(lines):
        dl = lines[i].strip()
        i += 1
        if not dl or dl.startswith('#'):
            continue
        parts = dl.split()
        if len(parts) >= 3:
            try:
                kdist_local.append(float(parts[0]))
                up.append(float(parts[1]))
                down.append(float(parts[2]))
                read += 1
            except:
                continue

    all_kdist.append(np.array(kdist_local))
    all_up.append(np.array(up))
    all_down.append(np.array(down))

print(f"  Parsed {band_count} bands")

if band_count == 0:
    print("ERROR: No band data")
    sys.exit(1)

# ======== 5. Parse KPOINTS for high-symmetry point labels ========
hs_positions = []
hs_labels = []
if os.path.exists(kpoints) and os.path.exists(poscar):
    try:
        with open(poscar) as f:
            plines = f.readlines()
        scale = float(plines[1].strip())
        lattice = np.zeros((3, 3))
        for i in range(3):
            lattice[i] = [float(x) for x in plines[2+i].split()[:3]]
        lattice *= scale

        volume = np.dot(lattice[0], np.cross(lattice[1], lattice[2]))
        recip = np.zeros((3, 3))
        recip[0] = 2*np.pi * np.cross(lattice[1], lattice[2]) / volume
        recip[1] = 2*np.pi * np.cross(lattice[2], lattice[0]) / volume
        recip[2] = 2*np.pi * np.cross(lattice[0], lattice[1]) / volume

        def kcart(k):
            return np.dot(k, recip)

        with open(kpoints) as f:
            klines = f.readlines()

        # Collect all segment endpoints in order
        path_k = []
        path_lbl = []
        i = 4
        while i < len(klines):
            line = klines[i].strip()
            i += 1
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            kvec = [float(x) for x in parts[:3]]
            label = line.split('!')[1].strip() if '!' in line else ''
            path_k.append(kvec)
            path_lbl.append(label)

        # Remove consecutive duplicates (e.g. G→M→K→G keeps all 4)
        uniq_k = [path_k[0]]
        uniq_lbl = [path_lbl[0]]
        for k, l in zip(path_k[1:], path_lbl[1:]):
            if any(abs(k[j] - uniq_k[-1][j]) > 1e-6 for j in range(3)):
                uniq_k.append(k)
                uniq_lbl.append(l)

        # Cumulative k-distances along the path
        hs_positions = [0.0]
        for i in range(1, len(uniq_k)):
            d = np.linalg.norm(kcart(uniq_k[i]) - kcart(uniq_k[i-1]))
            hs_positions.append(hs_positions[-1] + d)
        hs_labels = uniq_lbl

        print(f"  K-path: {', '.join(f'{l}={p:.3f}' for l, p in zip(hs_labels, hs_positions))}")
    except Exception as e:
        print(f"  Warning: KPOINTS parse failed ({e}), using default labels")

# ======== 6. Plot ========
fig, ax = plt.subplots(figsize=(8, 6))

for b in range(band_count):
    kd = all_kdist[b]
    ax.plot(kd, all_up[b] - efermi, 'r-', lw=0.8, alpha=0.8)
    ax.plot(kd, all_down[b] - efermi, 'b--', lw=0.8, alpha=0.8)

ax.axhline(y=0, color='gray', ls='--', lw=0.5, alpha=0.7)

if hs_positions:
    for pos in hs_positions:
        ax.axvline(x=pos, color='gray', ls='-', lw=0.5, alpha=0.5)
    ax.set_xticks(hs_positions)
    ax.set_xticklabels(hs_labels)
else:
    k_min, k_max = 0, max(kd[-1] for kd in all_kdist)
    for pos in [k_min, k_max]:
        ax.axvline(x=pos, color='gray', ls='-', lw=0.5, alpha=0.5)
    ax.set_xticks([k_min, k_max])
    ax.set_xticklabels(['M', 'K'])
ax.set_ylabel('E - E_F (eV)')
ax.set_title(f'Band Structure - {struct_name}')

ax.set_ylim(-3, 3)

ax.legend(
    [Line2D([0],[0], color='red', lw=1.5),
     Line2D([0],[0], color='blue', lw=1.5, ls='--')],
    ['Spin up', 'Spin down'],
    loc='upper right')

ax.grid(True, alpha=0.3)
plt.tight_layout()

pdf = os.path.join(work_dir, f'{struct_name}.pdf')
jpg = os.path.join(work_dir, f'{struct_name}.jpg')
plt.savefig(pdf, dpi=150)
plt.savefig(jpg, dpi=200, format='jpg', bbox_inches='tight')
print(f'  Saved: {pdf}')
print(f'  Saved: {jpg}')
plt.close()
