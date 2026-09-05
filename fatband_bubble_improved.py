#!/usr/bin/env python3
"""
投影能带 fatband 图（气泡图版，期刊风格，全局统一缩放）- 改进版
- 横坐标与高对称点标签自动从 KLABELS 文件读取
- 纵坐标锁定在 [-1, 2] eV
- 气泡大小使用全局统一缩放（所有轨道可比）
- 自动解析轨道列名，多子图展示

改进点：
1. excluded 集合改为大小写不敏感，避免 Tot/TOT 列被误识别为轨道
2. 零权重点不显示，避免形成假能带线
3. 增加费米面 (E=0) 参考线
4. suptitle 位置调整，避免被裁切
5. 图例使用 fig.legend，避免压住子图
6. 保存后关闭 figure，避免内存累积
7. format_label 默认分支更安全
8. axes 包装使用 np.atleast_1d
9. read_pband 数据行长度校验更宽松
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ========================== 用户配置 ==========================
KLABELS_FILE = "KLABELS"
PBAND_FILES = [
    "PBAND_Mn_UP.dat",
    # "PBAND_Mn_DW.dat",
]
YLIM = (-1, 2)
SAVE_FIG = True
OUTPUT_PREFIX = "fatband_bubble"
FIG_FORMAT = "pdf"
DPI = 300

# 气泡图参数
BUBBLE_COLOR = "#1f4e79"
BUBBLE_ALPHA = 0.6
MIN_SIZE = 0.5          # 最小气泡大小 (pt^2)
MAX_SIZE = 40           # 最大气泡大小 (pt^2)
EDGECOLOR = "none"

# 零权重阈值：低于此值的权重不显示气泡
WEIGHT_THRESHOLD = 1e-4

# 是否绘制费米面参考线
DRAW_EFERMI_LINE = True
EFERMI = 0.0

# ========================== 期刊风格设置 ==========================
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

# ========================== 1. KLABELS 解析 ==========================
def read_klabels(filename):
    labels, coords = [], []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            label = parts[0]
            try:
                coord = float(parts[1])
            except ValueError:
                continue
            labels.append(label)
            coords.append(coord)
    return labels, coords

def format_label(raw):
    """格式化高对称点标签。未知标签原样返回，避免误解释下划线。"""
    mapping = {
        "GAMMA": r"$\Gamma$",
        "GAMMA_": r"$\Gamma'$",
        "M": "M",
        "M_": r"$\mathrm{M}'$",
        "K": "K",
        "K_": r"$\mathrm{K}'$",
        "L": "L",
        "L_": r"$\mathrm{L}'$",
        "X": "X",
        "X_": r"$\mathrm{X}'$",
    }
    return mapping.get(raw, raw)

# ========================== 2. 投影能带文件读取 ==========================
def read_pband(filename):
    with open(filename, 'r') as f:
        header = f.readline().strip()
    col_names = header.replace('#', '').split()

    try:
        idx_energy = col_names.index('Energy')
    except ValueError:
        idx_energy = next(i for i, n in enumerate(col_names) if 'energy' in n.lower())

    # 改进 1: 大小写不敏感排除 K-Path / Energy / Tot 等非轨道列
    excluded_lower = {'k-path', 'energy', 'tot'}
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
            # 改进 9: 长度校验更宽松，多字段时取前 len(col_names) 个
            if len(parts) >= len(col_names):
                current.append([float(x) for x in parts[:len(col_names)]])
            elif len(parts) == len(col_names):
                current.append([float(x) for x in parts])
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

# ========================== 3. 气泡图绘制（全局缩放） ==========================
def plot_fatband_bubble(ax, kpts, energies, weight, ylim,
                        high_sym_coords, high_sym_labels, title,
                        global_max_weight):
    nbands = energies.shape[0]
    if global_max_weight > 0:
        # 改进 2: 零权重不显示
        size = np.where(weight > WEIGHT_THRESHOLD,
                       MIN_SIZE + (weight / global_max_weight) * (MAX_SIZE - MIN_SIZE),
                       0)
    else:
        size = np.full_like(weight, MIN_SIZE)

    for ib in range(nbands):
        e = energies[ib, :]
        w = weight[ib, :]
        s = size[ib, :]
        # 同时按能量范围和权重阈值过滤
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

    # 高对称点垂直线
    for coord in high_sym_coords:
        ax.axvline(x=coord, color='gray', linestyle='--', linewidth=0.6, alpha=0.7)

    # 改进 3: 费米面参考线
    if DRAW_EFERMI_LINE and ylim[0] <= EFERMI <= ylim[1]:
        ax.axhline(y=EFERMI, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.set_xticks(high_sym_coords)
    ax.set_xticklabels([format_label(lab) for lab in high_sym_labels])
    ax.set_ylabel("Energy (eV)")
    if title:
        ax.set_title(title, fontsize=10)

# ========================== 4. 气泡大小图例 ==========================
def add_bubble_legend(fig, max_weight):
    """改进 5: 使用 fig.legend，自动避免压住子图"""
    if max_weight <= 0:
        return
    weights_ref = np.linspace(0, max_weight, 3)
    sizes_ref = MIN_SIZE + (weights_ref / max_weight) * (MAX_SIZE - MIN_SIZE)

    handles = []
    labels = []
    for w, s in zip(weights_ref, sizes_ref):
        h = plt.scatter([], [], s=s, c=BUBBLE_COLOR, alpha=BUBBLE_ALPHA,
                        edgecolors='none')
        handles.append(h)
        labels.append(f"{w:.2f}")

    fig.legend(handles, labels,
               title="Weight",
               loc='center left',
               bbox_to_anchor=(0.99, 0.5),
               frameon=False,
               scatterpoints=1,
               fontsize=8,
               title_fontsize=9,
               labelspacing=1.5,
               borderpad=1.0)

# ========================== 5. 主程序 ==========================
def main():
    if not os.path.isfile(KLABELS_FILE):
        print(f"错误：找不到 {KLABELS_FILE}")
        return
    sym_labels, sym_coords = read_klabels(KLABELS_FILE)
    print(f"高对称点: {list(zip(sym_labels, sym_coords))}")

    for pband_file in PBAND_FILES:
        if not os.path.isfile(pband_file):
            print(f"警告：{pband_file} 不存在，跳过。")
            continue
        print(f"正在处理 {pband_file} ...")
        orbital_names, kpts, energies, weights = read_pband(pband_file)
        n_orbs = len(orbital_names)
        print(f"  轨道: {orbital_names}")

        # ---------- 计算全局最大权重 ----------
        global_max = 0.0
        for orb in orbital_names:
            max_w = np.max(weights[orb])
            if max_w > global_max:
                global_max = max_w

        # ---------- 创建子图 ----------
        cols = int(np.ceil(np.sqrt(n_orbs)))
        rows = int(np.ceil(n_orbs / cols))
        # 改进 5: 右侧预留空间给图例
        fig, axes = plt.subplots(rows, cols,
                                 figsize=(3.2*cols + 0.8, 2.8*rows),
                                 sharex=True, sharey=True)
        # 改进 8: 统一 axes 包装
        axes = np.atleast_1d(axes).flatten()

        # ---------- 绘制所有轨道 ----------
        for idx, orb in enumerate(orbital_names):
            plot_fatband_bubble(axes[idx], kpts, energies, weights[orb],
                                ylim=YLIM,
                                high_sym_coords=sym_coords,
                                high_sym_labels=sym_labels,
                                title=orb,
                                global_max_weight=global_max)

        # 隐藏多余子图
        for idx in range(n_orbs, len(axes)):
            axes[idx].set_visible(False)

        # 添加图例
        add_bubble_legend(fig, global_max)

        # 改进 4: suptitle 位置调整
        fig.suptitle(f"Projected Band Structure: {pband_file}",
                     fontsize=12, y=0.98)
        fig.subplots_adjust(top=0.92, right=0.92)

        if SAVE_FIG:
            base = os.path.splitext(os.path.basename(pband_file))[0]
            out_name = f"{OUTPUT_PREFIX}_{base}.{FIG_FORMAT}"
            plt.savefig(out_name, format=FIG_FORMAT, dpi=DPI)
            print(f"  图片已保存至 {out_name}")
            # 改进 6: 保存后关闭
            plt.close(fig)
        else:
            plt.show()

if __name__ == "__main__":
    main()
