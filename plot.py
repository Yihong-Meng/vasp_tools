#!/usr/bin/env python3
"""
投影能带 fatband 图（气泡图版，期刊风格，全局统一缩放）- 改进版 v2
- 横坐标与高对称点标签自动从 KLABELS 文件读取
- 纵坐标锁定在 [-1, 2] eV
- 气泡大小使用全局统一缩放（所有轨道可比）
- 自动解析轨道列名，多子图展示
- 支持自旋极化 UP/DW 同图不同色展示（蓝/红）

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
10. [新增] 支持自旋极化 UP/DW 同图不同色（蓝/红），全局最大权重跨自旋计算
11. [新增] 极小权重更不显眼：size 用 power 非线性映射 + 极小权重单独低 alpha
12. [新增] Origin 风格：底层细能带骨架线 + 气泡叠加，分布更连贯；取消硬阈值过滤避免断裂
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import re

# ========================== 用户配置 ==========================
KLABELS_FILE = "KLABELS"

# 自旋极化文件配对：(UP文件, DW文件)
# DW 文件可为 None（仅画 UP）；UP 文件可为 None（仅画 DW）
PBAND_PAIRS = [
    ("PBAND_Mn_UP.dat", "PBAND_Mn_DW.dat"),
    # ("PBAND_Fe_UP.dat", "PBAND_Fe_DW.dat"),
    # ("PBAND_O_UP.dat", None),   # 仅 UP
]

YLIM = (-1, 2)
SAVE_FIG = True
OUTPUT_PREFIX = "fatband_bubble"
FIG_FORMAT = "pdf"
DPI = 300

# 气泡图参数 - 自旋极化颜色
BUBBLE_COLOR_UP = "#1f4e79"   # 蓝（UP / majority spin）
BUBBLE_COLOR_DW = "#c0392b"   # 红（DW / minority spin）
BUBBLE_ALPHA = 0.6            # 正常权重透明度
BUBBLE_ALPHA_SMALL = 0.2      # 极小权重透明度（更不显眼）
MIN_SIZE = 0.1                # 最小气泡大小 (pt^2)，降低以让小权重更不显眼
MAX_SIZE = 40                 # 最大气泡大小 (pt^2)
EDGECOLOR = "none"

# 零权重阈值：低于此值的权重不显示气泡
WEIGHT_THRESHOLD = 1e-4

# 极小权重阈值（占全局最大权重的比例）：低于此值用更小 size + 更低 alpha
SMALL_WEIGHT_FRACTION = 0.05  # 5%

# size 非线性映射指数：>1 让小权重更小、大权重更突出（推荐 1.5~2.5）
SIZE_POWER = 2.0

# 底层能带骨架线（Origin 风格：让气泡分布连贯）
DRAW_SKELETON_LINE = True
SKELETON_COLOR = "lightgray"
SKELETON_LINEWIDTH = 0.5
SKELETON_ALPHA = 0.6

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
        if current:
            data_blocks.append(np.array(current))

    kpts = data_blocks[0][:, 0]
    energies = np.array([block[:, idx_energy] for block in data_blocks])

    weights = {}
    for orb in orbital_names:
        idx = col_names.index(orb)
        weights[orb] = np.array([block[:, idx] for block in data_blocks])

    return orbital_names, kpts, energies, weights

# ========================== 3. 气泡图绘制 ==========================
def plot_fatband_single_spin(ax, kpts, energies, weight, color,
                              ylim, global_max_weight):
    """画单个自旋方向的气泡（Origin 风格：底层骨架线 + 气泡叠加）
    改进：
    - 底层画细能带线作为骨架，保证视觉连贯
    - 取消硬阈值过滤，极小权重也画（用极小 size + 极低 alpha），避免断裂
    - size 用 power 非线性映射，小权重更小
    - 极小权重单独低 alpha
    """
    nbands = energies.shape[0]
    if global_max_weight > 0:
        w_norm = weight / global_max_weight
        # 非线性映射：power > 1 让小权重更小
        size = MIN_SIZE + (w_norm ** SIZE_POWER) * (MAX_SIZE - MIN_SIZE)
    else:
        size = np.full_like(weight, MIN_SIZE)
        w_norm = np.zeros_like(weight)

    # 极小权重阈值（绝对值）
    small_threshold = SMALL_WEIGHT_FRACTION * global_max_weight

    for ib in range(nbands):
        e = energies[ib, :]
        s = size[ib, :]
        w = weight[ib, :]

        # 能量范围 mask（不再用 s>0 过滤，保证连贯）
        in_range = (e >= ylim[0]) & (e <= ylim[1])
        if not np.any(in_range):
            continue

        # ---- 1. 底层骨架线（Origin 风格）----
        if DRAW_SKELETON_LINE:
            ax.plot(kpts[in_range], e[in_range],
                    color=SKELETON_COLOR, linewidth=SKELETON_LINEWIDTH,
                    alpha=SKELETON_ALPHA, zorder=1, rasterized=True)

        # ---- 2. 气泡叠加 ----
        # 分两组：极小权重 vs 正常权重
        small_mask = in_range & (w < small_threshold)
        normal_mask = in_range & (w >= small_threshold)

        # 先画极小权重（低 alpha，更不显眼）
        if np.any(small_mask):
            ax.scatter(kpts[small_mask], e[small_mask],
                       s=s[small_mask], c=color, alpha=BUBBLE_ALPHA_SMALL,
                       edgecolors=EDGECOLOR, linewidth=0,
                       zorder=2, rasterized=True)
        # 再画正常权重（正常 alpha）
        if np.any(normal_mask):
            ax.scatter(kpts[normal_mask], e[normal_mask],
                       s=s[normal_mask], c=color, alpha=BUBBLE_ALPHA,
                       edgecolors=EDGECOLOR, linewidth=0,
                       zorder=3, rasterized=True)

def plot_fatband_bubble(ax, spin_data, ylim,
                        high_sym_coords, high_sym_labels, title,
                        global_max_weight):
    """画一个子图，可包含 UP 和 DW
    spin_data: list of dict, 每个含 kpts, energies, weight, color, label
    """
    # 画各自旋方向
    for sd in spin_data:
        plot_fatband_single_spin(ax, sd['kpts'], sd['energies'],
                                  sd['weight'], sd['color'],
                                  ylim, global_max_weight)

    # 坐标轴范围（用第一个自旋的 kpts）
    kpts0 = spin_data[0]['kpts']
    ax.set_xlim(kpts0[0], kpts0[-1])
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

    # 改进 10: 多自旋时加颜色图例
    if len(spin_data) > 1:
        legend_elements = [
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=sd['color'], markeredgecolor='none',
                   markersize=7, label=sd['label'], alpha=BUBBLE_ALPHA)
            for sd in spin_data
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
                  handletextpad=0.3, borderpad=0.3)

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
        h = plt.scatter([], [], s=s, c=BUBBLE_COLOR_UP, alpha=BUBBLE_ALPHA,
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

# ========================== 5. 辅助函数 ==========================
def extract_element_name(filename):
    """从文件名提取元素标识，如 PBAND_Mn_UP.dat -> Mn"""
    base = os.path.basename(filename)
    m = re.search(r'PBAND[_-]([A-Za-z]+)[_-](UP|DW)', base)
    if m:
        return m.group(1)
    # 兜底：去掉 PBAND_ 前缀和扩展名
    m = re.search(r'PBAND[_-]([A-Za-z]+)', base)
    if m:
        return m.group(1)
    return os.path.splitext(base)[0]

# ========================== 6. 主程序 ==========================
def main():
    if not os.path.isfile(KLABELS_FILE):
        print(f"错误：找不到 {KLABELS_FILE}")
        return
    sym_labels, sym_coords = read_klabels(KLABELS_FILE)
    print(f"高对称点: {list(zip(sym_labels, sym_coords))}")

    for up_file, dw_file in PBAND_PAIRS:
        # ---------- 读取 UP ----------
        up_data = None
        if up_file:
            if os.path.isfile(up_file):
                print(f"正在处理 {up_file} ...")
                up_data = read_pband(up_file)
            else:
                print(f"警告：{up_file} 不存在")

        # ---------- 读取 DW ----------
        dw_data = None
        if dw_file:
            if os.path.isfile(dw_file):
                print(f"正在处理 {dw_file} ...")
                dw_data = read_pband(dw_file)
            else:
                print(f"警告：{dw_file} 不存在")

        if up_data is None and dw_data is None:
            print("  UP 和 DW 都不存在，跳过。")
            continue

        # ---------- 轨道名（UP 优先，否则 DW） ----------
        orbital_names = up_data[0] if up_data is not None else dw_data[0]
        n_orbs = len(orbital_names)
        print(f"  轨道: {orbital_names}")

        # ---------- 全局最大权重（跨 UP/DW） ----------
        global_max = 0.0
        for orb in orbital_names:
            if up_data is not None:
                global_max = max(global_max, np.max(up_data[3][orb]))
            if dw_data is not None:
                global_max = max(global_max, np.max(dw_data[3][orb]))
        print(f"  全局最大权重: {global_max:.4f}")

        # ---------- 创建子图 ----------
        cols = int(np.ceil(np.sqrt(n_orbs)))
        rows = int(np.ceil(n_orbs / cols))
        # 改进 5: 右侧预留空间给图例
        fig, axes = plt.subplots(rows, cols,
                                 figsize=(3.2*cols + 0.8, 2.8*rows),
                                 sharex=True, sharey=True)
        # 改进 8: 统一 axes 包装
        axes = np.atleast_1d(axes).flatten()

        # ---------- 绘制每个轨道 ----------
        for idx, orb in enumerate(orbital_names):
            spin_data = []
            if up_data is not None:
                spin_data.append({
                    'kpts': up_data[1],
                    'energies': up_data[2],
                    'weight': up_data[3][orb],
                    'color': BUBBLE_COLOR_UP,
                    'label': 'UP'
                })
            if dw_data is not None:
                spin_data.append({
                    'kpts': dw_data[1],
                    'energies': dw_data[2],
                    'weight': dw_data[3][orb],
                    'color': BUBBLE_COLOR_DW,
                    'label': 'DW'
                })
            plot_fatband_bubble(axes[idx], spin_data, ylim=YLIM,
                                high_sym_coords=sym_coords,
                                high_sym_labels=sym_labels,
                                title=orb,
                                global_max_weight=global_max)

        # 隐藏多余子图
        for idx in range(n_orbs, len(axes)):
            axes[idx].set_visible(False)

        # 添加权重图例
        add_bubble_legend(fig, global_max)

        # ---------- 标题与保存 ----------
        ref_file = up_file or dw_file
        element = extract_element_name(ref_file)
        if up_data is not None and dw_data is not None:
            spin_label = "UP + DW"
        elif up_data is not None:
            spin_label = "UP"
        else:
            spin_label = "DW"

        # 改进 4: suptitle 位置调整
        fig.suptitle(f"Projected Band Structure: {element} ({spin_label})",
                     fontsize=12, y=0.98)
        fig.subplots_adjust(top=0.92, right=0.92)

        if SAVE_FIG:
            out_name = f"{OUTPUT_PREFIX}_{element}.{FIG_FORMAT}"
            plt.savefig(out_name, format=FIG_FORMAT, dpi=DPI)
            print(f"  图片已保存至 {out_name}")
            # 改进 6: 保存后关闭
            plt.close(fig)
        else:
            plt.show()

if __name__ == "__main__":
    main()
