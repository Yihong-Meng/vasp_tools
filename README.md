# vasp_tools
VASP post-processing and workflow automation toolkit: DFT/Wannier band structures, Fatband projections, optical properties (Kubo/JDOS), NEB analysis, and interface sliding generation.

## 📌 Disclaimer & Introduction

This repository contains a collection of scripts I have written and accumulated during my daily **VASP (Vienna Ab initio Simulation Package)** calculations. These are **personal usage records and notes** rather than polished commercial software.

They were born out of the need to automate repetitive tasks, visualize complex data (bands, optics, NEB), and manage high-throughput sliding calculations for 2D heterostructures.

> ⚠️ **Legal Notice**: This repository **only** contains user-generated analysis scripts and job submission templates. It **does NOT** include the VASP source code, pseudopotentials (POTCAR), or any other proprietary VASP files. Users must have their own VASP license to utilize the generated input files.

---

## 📂 Script Categories & Usage

### 1. Electronic Band Structure (DFT & Wannier90)
*Use case: Extracting DFT bands, comparing with Wannier interpolation, and plotting Fatbands.*

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `plot_dft_wannier.py` / `plot_comparar.py` | Overlay DFT (solid) and Wannier (dashed) bands for spin-up/down channels. | Modify `FERMI_ENERGY` and `SPIN_CHANNEL` in the script, then run `python3 plot_comparar.py`. |
| `plot_Spinsplit_band.py` / `plot_band_.py` | Quick spin-split band plotting from `BAND.dat` or `EIGENVAL`. | If `BAND.dat` is missing, the script auto-generates it from `EIGENVAL`. |
| `fatband_color.py` / `fatband_bubble_improved.py` | **Projected (Fatband) visualization**. First uses yellow→red colormaps; second uses bubble sizes to represent orbital weights (publication-ready). | Prepare `KLABELS` and `PBAND_*.dat` files, then configure paths and run. |
| `engwin.py` / `analyze_eigenval.py` | Extract specific band energies from `EIGENVAL` and count bands around the Fermi level. | `python3 engwin.py EIGENVAL e 10` to extract band #10. |

### 2. Optical Properties
*Use case: Analyzing Kubo conductivity and JDOS from Wannier90 outputs.*

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `plot_optical.py` | Plots the imaginary dielectric function (ε₂), optical conductivity (σ), and JDOS. **Feature**: Adds a rainbow background for the visible light region (1.65–3.26 eV). | Place the script in the folder with `wannier90-kubo_*.dat` and run `python3 plot_optical.py`. |

### 3. NEB Transition State
*Use case: Visualizing Climbing Image NEB results.*

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `plot_neb_energy.py` | Plots relative energy (in meV) along the normalized reaction coordinate (0→1). **Auto-detects the Transition State (TS)** and calculates activation energy. | Paste your NEB energy data into the `data` array within the script and run it. |

### 4. High-Throughput Sliding Workflow (Interface Sliding)
*Use case: Scanning the interlayer sliding potential energy surface for 2D heterostructures (e.g., Mn₂BrI/BN). This is my most frequently used workflow.*

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `gen_slip.py` / `gen_slip_middle_bn.py` | **Structure Generation**: Automatically detects the stacking axis (Z) and generates 8×8=64 slid `POSCAR` files along a1/a2. | `python3 gen_slip.py -p POSCAR -d 7` |
| `setup_slip.py` | **Environment Setup**: Creates independent calculation folders (`y_x`) and symlinks `POTCAR`/`KPOINTS`. | `python3 setup_slip.py` |
| `submit_slip.sh` | **Batch Submission**: SLURM controller with **max concurrent job limit (default 4)**. Supports auto-resubmission of failed jobs. | `nohup bash submit_slip.sh > submit.log 2>&1 &` |
| `check_results_v2.py` | **Result Diagnosis**: Batch checks `log` files for convergence flags and identifies VASP errors like `ZBRENT`. | `python3 check_results_v2.py` |
| `change_IBRION.py` / `potim_Iibrion.py` | **Auto-Restart**: Scans for unconverged folders and switches `IBRION=2→1` while lowering `POTIM=0.2` to stabilize optimization. | Run `python3 change_IBRION.py --dry-run` to preview, then remove `--dry-run` to execute. |
| `extract_slip_energy.py` | **Energy Extraction**: Extracts final energies from all `OUTCAR` files and outputs a relative energy matrix (in meV). | `python3 extract_slip_energy.py -o energy.csv` |

### 5. Phonon Workflow
*Use case: Batch processing for Phonopy finite-displacement method.*

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `setup_phono.sh` | Batch creates `disp-*` directories and symlinks VASP input files. | Run `bash setup_phono.sh` after preparing `POSCAR-*` files. |
| `submit_phonon.sh` | Batch submits VASP jobs in `disp-*` directories with a job limit. | `nohup bash submit_phonon.sh > phonon.log 2>&1 &` |
| `plot_phonon_GKMG.py` | Reads Phonopy's `band.yaml` and plots high-quality phonon dispersion with a secondary `cm⁻¹` y-axis. | `python3 plot_phonon_GKMG.py` (requires `pyyaml`). |

### 6. Utilities

| Script | Description | Basic Usage |
| :--- | :--- | :--- |
| `sort_poscar.py` | **POSCAR Atom Sorter**: Sorts atoms by x/y/z coordinates (within species or globally). | `python3 sort_poscar.py POSCAR -a z -o POSCAR_sorted` |
| `vasp_ab.sh` | A template SLURM submission script with Intel oneAPI support, suitable for hexagonal lattice parameter scanning. | Modify the `a` loop range and submit via `sbatch vasp_ab.sh`. |

---

## 🛠️ Dependencies

- **System**: Linux cluster with SLURM workload manager.
- **Python**: 3.0+
- **Python Libraries**:
  ```bash
  pip install numpy matplotlib pyyaml
