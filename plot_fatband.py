#!/usr/bin/env python3
"""
投影能带 Fatband 气泡图 - 自动批量处理自旋通道版
支持按 SPIN_CHANNEL 自动搜索所有 PBAND_*_UP/DW.dat 文件，逐个绘图。
"""

import numpy as np
import os
import glob

# ============================================================================
# ################################ 用户配置区 ################################
# ============================================================================

KLABELS_FILE = "KLABELS"

# ---------- 自旋通道与文件自动查找 ----------
SPIN_CHANNEL = "UP"           # 可选 "UP" 或 "DW"，脚本自动查找 PBAND_*_UP.dat 或 PBAND_*_DW.dat
# 如果仍想手动指定文件列表，请取消注释下面一行并注释掉 SPIN_CHANNEL
# PBAND_FILES = ["PBAND_N_UP.dat", "PBAND_Mn_UP.dat"]

# ---------- 轨道过滤设置 ----------
AUTO_D_ORBITALS = False          # 改为 False，使用手动过滤
ORBITAL_FILTER = ["s"]   # 手动指定轨道列名
OUTPUT_SUFFIX = "s"              # 输出后缀

# ---------- 其他参数 ----------
YLIM = (-4, 6)
SAVE_FIG = True
OUTPUT_PREFIX = "fatband_"
FIG_FORMAT = "png"
DPI = 400

# 气泡视觉参数（优化连续感）
BUBBLE_COLOR = "#1f4e79"
BUBBLE_ALPHA = 0.35
MIN_SIZE = 0.5
MAX_SIZE = 80
EDGECOLOR = "none"
WEIGHT_THRESHOLD = 1e-6

DRAW_EFERMI_LINE = True
EFERMI = 0.0

# ############################################################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "mathtext.default": "regular",
    "axes.linewidth": 0.8,
    "xtick.major.size": 4,
    "xtick.major.width": 0.8,
    "ytick.major.size": 4,
    "ytick.major.width": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

# ========================== KLABELS 解析 ==========================
def read_klabels(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    start_idx = 1
    if len(lines) > 1:
        first_parts = lines[1].split()
        if len(first_parts) == 1 and first_parts[0].isdigit():
            start_idx = 2
    labels, coords = [], []
    for line in lines[start_idx:]:
        parts = line.split()
        if len(parts) >= 2:
            try:
                coords.append(float(parts[1]))
                labels.append(parts[0])
            except ValueError:
                continue
    return labels, coords

def format_label(raw):
    mapping = {
        "GAMMA": r"$\Gamma$", "GAMMA_": r"$\Gamma'$",
        "M": "M", "M_": r"$\mathrm{M}'$",
        "K": "K", "K_": r"$\mathrm{K}'$",
        "L": "L", "L_": r"$\mathrm{L}'$",
        "X": "X", "X_": r"$\mathrm{X}'$",
    }
    return mapping.get(raw, raw)

# ========================== PBAND 读取 ==========================
def read_pband(filename):
    with open(filename, 'r') as f:
        header = f.readline().strip()
    col_names = header.replace('#', '').split()
    idx_energy = -1
    try:
        idx_energy = col_names.index('Energy')
    except ValueError:
        for i, n in enumerate(col_names):
            if 'energy' in n.lower():
                idx_energy = i
                break
    if idx_energy == -1:
        raise RuntimeError(f"在 {filename} 中找不到 Energy 列！")
    excluded_lower = {'k-path', 'energy', 'tot', 'band', 'index'}
    orbital_names = [n for n in col_names if n.lower() not in excluded_lower]
    data_blocks, current = [], []
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('# Band-Index'):
                if current:
                    data_blocks.append(np.array(current))
                    current = []
                continue
            parts = line.split()
            if len(parts) >= len(col_names):
                current.append([float(x) for x in parts[:len(col_names)]])
        if current:
            data_blocks.append(np.array(current))
    nkpts = data_blocks[0].shape[0]
    kpts = data_blocks[0][:, 0]
    energies = np.array([block[:, idx_energy] for block in data_blocks])
    weights = {}
    for orb in orbital_names:
        idx = col_names.index(orb)
        weights[orb] = np.array([block[:, idx] for block in data_blocks])
    return orbital_names, kpts, energies, weights

# ========================== 绘图函数 ==========================
def plot_fatband_bubble(ax, kpts, energies, weight, ylim,
                        high_sym_coords, high_sym_labels, title,
                        global_max_weight):
    nbands = energies.shape[0]
    if global_max_weight > 0:
        size = np.where(weight > WEIGHT_THRESHOLD,
                       MIN_SIZE + (weight / global_max_weight) * (MAX_SIZE - MIN_SIZE),
                       0)
    else:
        size = np.full_like(weight, MIN_SIZE)
    for ib in range(nbands):
        e = energies[ib, :]
        w = weight[ib, :]
        s = size[ib, :]
        mask = (e >= ylim[0]) & (e <= ylim[1]) & (s > 0)
        if not np.any(mask):
            continue
        ax.scatter(kpts[mask], e[mask],
                   s=s[mask],
                   c=BUBBLE_COLOR,
                   alpha=BUBBLE_ALPHA,
                   edgecolors=EDGECOLOR,
                   linewidth=0,
                   rasterized=True)
    ax.set_xlim(kpts[0], kpts[-1])
    ax.set_ylim(ylim)
    for coord in high_sym_coords:
        ax.axvline(x=coord, color='gray', linestyle='--', linewidth=0.6, alpha=0.7)
    if DRAW_EFERMI_LINE and ylim[0] <= EFERMI <= ylim[1]:
        ax.axhline(y=EFERMI, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_xticks(high_sym_coords)
    ax.set_xticklabels([format_label(lab) for lab in high_sym_labels])
    ax.set_ylabel("Energy (eV)")
    if title:
        ax.set_title(title, fontsize=10)

def add_bubble_legend(fig, max_weight):
    if max_weight <= 0:
        return
    weights_ref = np.linspace(0, max_weight, 3)
    sizes_ref = MIN_SIZE + (weights_ref / max_weight) * (MAX_SIZE - MIN_SIZE)
    fmt = '{:.2f}' if max_weight < 10 else '{:.0f}'
    handles, labels = [], []
    for w, s in zip(weights_ref, sizes_ref):
        h = plt.scatter([], [], s=s, c=BUBBLE_COLOR, alpha=BUBBLE_ALPHA, edgecolors='none')
        handles.append(h)
        labels.append(fmt.format(w))
    fig.legend(handles, labels,
               title="Weight",
               loc='upper center',
               bbox_to_anchor=(0.5, 1.12),
               frameon=False,
               scatterpoints=1,
               fontsize=8,
               title_fontsize=9,
               labelspacing=1.5,
               borderpad=1.0,
               ncol=3)

# ========================== 主程序 ==========================
def main():
    # ---------- 确定要处理的文件列表 ----------
    if 'PBAND_FILES' in globals() and PBAND_FILES:
        pband_files = PBAND_FILES
    else:
        pattern = f"PBAND_*_{SPIN_CHANNEL}.dat"
        pband_files = sorted(glob.glob(pattern))
        if not pband_files:
            print(f"错误：未找到任何匹配 {pattern} 的文件！")
            return
        print(f"根据自旋通道 {SPIN_CHANNEL} 自动找到以下文件：")
        for f in pband_files:
            print(f"  {f}")
        print("")

    if not os.path.isfile(KLABELS_FILE):
        print(f"错误：找不到 {KLABELS_FILE}")
        return
    sym_labels, sym_coords = read_klabels(KLABELS_FILE)
    print(f"高对称点: {list(zip(sym_labels, sym_coords))}")

    for pband_file in pband_files:
        if not os.path.isfile(pband_file):
            print(f"警告：{pband_file} 不存在，跳过。")
            continue
        print(f"\n正在处理 {pband_file} ...")
        orbital_names_full, kpts, energies, weights = read_pband(pband_file)
        print(f"  读取到的全部轨道: {orbital_names_full}")

        # ---------- 应用过滤 ----------
        suffix = ""
        if AUTO_D_ORBITALS:
            d_orbitals = []
            for orb in orbital_names_full:
                if 'd' in orb.lower() or orb.lower() in ['x2-y2', 'dx2-y2', 'x2y2']:
                    d_orbitals.append(orb)
            if not d_orbitals:
                print("  警告：自动识别未找到任何 d 轨道，将绘制全部轨道。")
                orbital_names = orbital_names_full
                weights = {orb: weights[orb] for orb in orbital_names}
                suffix = "all"
            else:
                orbital_names = d_orbitals
                weights = {orb: weights[orb] for orb in orbital_names}
                print(f"  自动识别到 d 轨道: {orbital_names}")
                suffix = "d"
        else:
            if ORBITAL_FILTER:
                filtered = []
                for orb in orbital_names_full:
                    for f in ORBITAL_FILTER:
                        if f in orb:
                            filtered.append(orb)
                            break
                filtered_ordered = []
                for orb in orbital_names_full:
                    if orb in filtered and orb not in filtered_ordered:
                        filtered_ordered.append(orb)
                if not filtered_ordered:
                    print(f"  警告：过滤器 {ORBITAL_FILTER} 未匹配任何轨道，将绘制全部轨道！")
                    orbital_names = orbital_names_full
                    weights = {orb: weights[orb] for orb in orbital_names}
                    suffix = "all"
                else:
                    orbital_names = filtered_ordered
                    weights = {orb: weights[orb] for orb in orbital_names}
                    print(f"  手动过滤后保留轨道: {orbital_names}")
                    if any('p' in f for f in ORBITAL_FILTER):
                        suffix = "p"
                    elif any('s' in f for f in ORBITAL_FILTER):
                        suffix = "s"
                    else:
                        suffix = "custom"
            else:
                orbital_names = orbital_names_full
                print("  未使用过滤器，绘制全部轨道。")
                suffix = "all"

        if OUTPUT_SUFFIX:
            suffix = OUTPUT_SUFFIX

        n_orbs = len(orbital_names)
        if n_orbs == 0:
            print("  错误：没有可绘制的轨道，跳过该文件。")
            continue

        # ---------- 全局最大权重 ----------
        global_max = 0.0
        for orb in orbital_names:
            max_w = np.max(weights[orb])
            if max_w > global_max:
                global_max = max_w
        print(f"  全局最大权重: {global_max:.4f} (所有轨道统一以此为基准)")

        # ---------- 创建子图 ----------
        cols = int(np.ceil(np.sqrt(n_orbs)))
        rows = int(np.ceil(n_orbs / cols))

        if n_orbs == 1:
            figsize = (6, 5)
        else:
            figsize = (3.2 * cols + 0.8, 2.8 * rows)

        fig, axes = plt.subplots(rows, cols,
                                 figsize=figsize,
                                 sharex=True, sharey=True)
        axes = np.atleast_1d(axes).flatten()

        for idx, orb in enumerate(orbital_names):
            plot_fatband_bubble(axes[idx], kpts, energies, weights[orb],
                                ylim=YLIM,
                                high_sym_coords=sym_coords,
                                high_sym_labels=sym_labels,
                                title=orb,
                                global_max_weight=global_max)

        for idx in range(n_orbs, len(axes)):
            axes[idx].set_visible(False)

        add_bubble_legend(fig, global_max)

        orb_list_str = ', '.join(orbital_names)
        fig.suptitle(f"Projected Band: {os.path.basename(pband_file)} - {orb_list_str}",
                     fontsize=12, y=0.98)

        if n_orbs == 1:
            fig.subplots_adjust(top=0.80)
        else:
            fig.subplots_adjust(top=0.85)

        if SAVE_FIG:
            base = os.path.splitext(os.path.basename(pband_file))[0]
            suffix_str = f"_{suffix}" if suffix else ""
            out_name = f"{OUTPUT_PREFIX}_{base}{suffix_str}.{FIG_FORMAT}"
            plt.savefig(out_name, format=FIG_FORMAT, dpi=DPI)
            print(f"  图片已保存至 {out_name}")
            plt.close(fig)
        else:
            plt.show()

if __name__ == "__main__":
    main()
