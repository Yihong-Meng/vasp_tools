#!/usr/bin/env python3
"""
投影能带 fatband 图（颜色版，带完整能带背景）
- 从 BAND.dat 读取完整能带，用浅色绘制作为背景
- 从 PBAND_*.dat 读取投影能带，用深色绘制叠加在背景上
- 权重高：深色/显眼颜色（深红色）
- 权重低：浅色/淡颜色（黄色）
- 右侧显示颜色直方图（colorbar）

使用方法：
python fatband_color.py

需要的文件：
1. BAND.dat 文件（完整能带数据）
2. PBAND_*.dat 文件（投影能带数据）
3. KLABELS 文件（高对称点标签）
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端，避免显示错误
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import os

# ========================== 用户配置 ==========================
# 输入文件
KLABELS_FILE = "KLABELS"
BAND_FILE = "BAND.dat"             # 完整能带文件
PBAND_FILES = [
    "PBAND_Mn_UP.dat",
    # "PBAND_Mn_DW.dat",
    # "PBAND_Br_UP.dat",
    # "PBAND_I_UP.dat",
]

# 输出设置
YLIM = (-8, 8)                    # Y 轴范围 (eV)
SAVE_FIG = True
OUTPUT_PREFIX = "fatband_color"
FIG_FORMAT = "png"
DPI = 300

# 颜色设置（参考图片：黄色→橙色→红色→深红色）
# 自定义颜色渐变，从浅黄到深红
CUSTOM_COLORS = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", 
                 "#fc4e2a", "#e31a1c", "#bd0026", "#800026"]
COLORMAP = LinearSegmentedColormap.from_list("fatband_custom", CUSTOM_COLORS, N=256)

# 颜色范围
WEIGHT_MIN_COLOR = 0.0            # 最小权重对应的颜色强度 (0-1)
WEIGHT_MAX_COLOR = 1.0            # 最大权重对应的颜色强度 (0-1)

# 绘图参数
BAND_LINEWIDTH = 1.2               # 能带线宽
WEIGHT_THRESHOLD = 1e-4            # 权重低于此值不绘制
EFERMI = 0.0                       # 费米能级位置
DRAW_EFERMI_LINE = True            # 是否绘制费米能级线

# 完整能带背景设置
BACKGROUND_COLOR = "#d9d9d9"       # 背景能带颜色（浅灰色）
BACKGROUND_LINEWIDTH = 0.8         # 背景能带线宽
BACKGROUND_ALPHA = 0.6             # 背景能带透明度

# 坐标轴设置
Y_TICK_STEP = 2                    # Y 轴刻度步长 (eV)

# 颜色条设置
COLORBAR_LABEL = "Orbital Weight"
COLORBAR_FRACTION = 0.046          # 颜色条宽度比例
COLORBAR_PAD = 0.04                # 颜色条与图的间距

# 子图布局（自动计算）
FIGURE_WIDTH = 12.0                # 图总宽度
FIGURE_HEIGHT = 8.0                # 图总高度

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
    "backend": "Agg",  # 使用非交互式后端
})

# ========================== 1. KLABELS 解析 ==========================
def read_klabels(filename):
    """读取高对称点标签文件"""
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
    """格式化高对称点标签"""
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

# ========================== 2. 完整能带文件读取 ==========================
def read_band(filename, spin='up'):
    """从 BAND.dat 读取完整能带数据"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    nkpts = nbands = 0
    for line in lines[:10]:
        if 'NKPTS' in line and 'NBANDS' in line:
            parts = line.split(':')[1].strip().split()
            nkpts = int(parts[0])
            nbands = int(parts[1])
            break
    
    if nkpts == 0:
        raise ValueError(f"无法解析 {filename} 中的 NKPTS/NBANDS")
    
    k_list, e_list = [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if '# Band-Index' not in line:
            continue
        k_band, e_band = [], []
        while len(k_band) < nkpts and i < len(lines):
            data_line = lines[i].strip()
            i += 1
            if not data_line or data_line.startswith('#'):
                continue
            parts = data_line.split()
            if len(parts) < 3:
                continue
            try:
                k_band.append(float(parts[0]))
                if spin.lower() == 'up':
                    e_band.append(float(parts[1]))
                else:
                    e_band.append(float(parts[2]))
            except ValueError:
                continue
        if len(k_band) == nkpts:
            k_arr = np.array(k_band)
            e_arr = np.array(e_band)
            sort_idx = np.argsort(k_arr)
            k_list.append(k_arr[sort_idx])
            e_list.append(e_arr[sort_idx])
    
    return k_list, e_list

# ========================== 3. 投影能带文件读取 ==========================
def read_pband(filename):
    """读取投影能带数据文件"""
    with open(filename, 'r') as f:
        header = f.readline().strip()
    col_names = header.replace('#', '').split()

    try:
        idx_energy = col_names.index('Energy')
    except ValueError:
        idx_energy = next(i for i, n in enumerate(col_names) if 'energy' in n.lower())

    # 排除非轨道列
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

# ========================== 4. 颜色映射函数 ==========================
def create_colormap():
    """创建颜色映射"""
    return COLORMAP

def weight_to_color(weight, max_weight, cmap):
    """将权重转换为颜色"""
    if max_weight <= 0:
        return cmap(0.5)  # 返回中间颜色
    # 归一化权重到 [0, 1]
    norm_weight = np.clip(weight / max_weight, 0, 1)
    return cmap(norm_weight)

# ========================== 5. 绘制完整能带背景 ==========================
def plot_band_background(ax, kpts_list, energies_list, ylim):
    """
    绘制完整能带作为背景
    
    参数:
        ax: matplotlib axes 对象
        kpts_list: k 点坐标列表（每个元素是一个 k 点数组）
        energies_list: 能量列表（每个元素是一个能量数组）
        ylim: Y 轴范围
    """
    for kpts, energies in zip(kpts_list, energies_list):
        # 处理一维数组的情况（只有一个 k 点路径）
        if energies.ndim == 1:
            # 一维数组：只有一个能带在多个 k 点上
            e = energies
            mask = (e >= ylim[0]) & (e <= ylim[1])
            if np.any(mask):
                ax.plot(kpts[mask], e[mask], 
                       color=BACKGROUND_COLOR, 
                       linewidth=BACKGROUND_LINEWIDTH, 
                       alpha=BACKGROUND_ALPHA)
        else:
            # 二维数组：多个能带在多个 k 点上
            nbands = energies.shape[0]
            for ib in range(nbands):
                e = energies[ib, :]
                # 过滤能量范围
                mask = (e >= ylim[0]) & (e <= ylim[1])
                if np.any(mask):
                    ax.plot(kpts[mask], e[mask], 
                           color=BACKGROUND_COLOR, 
                           linewidth=BACKGROUND_LINEWIDTH, 
                           alpha=BACKGROUND_ALPHA)

# ========================== 6. 颜色 Fatband 绘制 ==========================
def plot_fatband_color(ax, kpts, energies, weight, ylim,
                       high_sym_coords, high_sym_labels, title,
                       global_max_weight, cmap):
    """
    绘制颜色表示的 fatband
    
    参数:
        ax: matplotlib axes 对象
        kpts: k 点坐标数组
        energies: 能量数组 (nbands × nkpts)
        weight: 权重数组 (nbands × nkpts)
        ylim: Y 轴范围
        high_sym_coords: 高对称点坐标
        high_sym_labels: 高对称点标签
        title: 子图标题
        global_max_weight: 全局最大权重（用于归一化）
        cmap: 颜色映射
    """
    nbands = energies.shape[0]
    
    # 绘制每条能带
    for ib in range(nbands):
        e = energies[ib, :]
        w = weight[ib, :]
        
        # 过滤能量范围
        mask = (e >= ylim[0]) & (e <= ylim[1])
        if not np.any(mask):
            continue
        
        k_masked = kpts[mask]
        e_masked = e[mask]
        w_masked = w[mask]
        
        # 对能带进行分段绘制，每段使用对应权重的颜色
        for i in range(len(k_masked) - 1):
            if w_masked[i] < WEIGHT_THRESHOLD and w_masked[i+1] < WEIGHT_THRESHOLD:
                continue  # 两点权重都太低，跳过
            
            # 计算线段的平均权重
            avg_weight = (w_masked[i] + w_masked[i+1]) / 2
            color = weight_to_color(avg_weight, global_max_weight, cmap)
            
            # 设置线段的 alpha 值基于权重
            alpha = 0.3 + 0.7 * (avg_weight / global_max_weight if global_max_weight > 0 else 0)
            alpha = np.clip(alpha, 0.3, 1.0)
            
            ax.plot([k_masked[i], k_masked[i+1]], 
                   [e_masked[i], e_masked[i+1]],
                   color=color, linewidth=BAND_LINEWIDTH, alpha=alpha)

# ========================== 7. 坐标轴设置 ==========================
def setup_axes(ax, ylim, high_sym_coords, high_sym_labels, show_xlabel=True):
    """
    设置坐标轴格式
    
    参数:
        ax: matplotlib axes 对象
        ylim: Y 轴范围
        high_sym_coords: 高对称点坐标
        high_sym_labels: 高对称点标签
        show_xlabel: 是否显示 x 轴标签
    """
    ax.set_xlim(high_sym_coords[0], high_sym_coords[-1])
    ax.set_ylim(ylim)
    
    # 设置 Y 轴刻度（以 2 为步长）
    y_ticks = np.arange(ylim[0], ylim[1] + 1, Y_TICK_STEP)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{int(y)}" for y in y_ticks])
    
    # 高对称点垂直线
    for coord in high_sym_coords:
        ax.axvline(x=coord, color='gray', linestyle='--', linewidth=0.6, alpha=0.7)
    
    # 费米面参考线
    if DRAW_EFERMI_LINE and ylim[0] <= EFERMI <= ylim[1]:
        ax.axhline(y=EFERMI, color='black', linestyle='--', linewidth=0.8, alpha=0.8)
    
    # 设置 x 轴刻度
    ax.set_xticks(high_sym_coords)
    ax.set_xticklabels([format_label(lab) for lab in high_sym_labels])
    
    # 设置标签
    ax.set_ylabel("E - E$_F$ (eV)")
    if show_xlabel:
        ax.set_xlabel("k-path")

# ========================== 8. 颜色条图例 ==========================
def add_colorbar(fig, ax, cmap, max_weight, label=COLORBAR_LABEL):
    """添加颜色条图例"""
    # 创建一个 scalar mappable
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_weight))
    sm.set_array([])
    
    # 添加颜色条
    cbar = fig.colorbar(sm, ax=ax, fraction=COLORBAR_FRACTION, pad=COLORBAR_PAD)
    cbar.set_label(label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    # 设置刻度
    cbar.set_ticks([0, max_weight * 0.5, max_weight])
    cbar.set_ticklabels([f'{0:.2f}', f'{max_weight * 0.5:.2f}', f'{max_weight:.2f}'])
    
    return cbar

# ========================== 9. 主程序 ==========================
def main():
    # 检查输入文件
    if not os.path.isfile(KLABELS_FILE):
        print(f"错误：找不到 {KLABELS_FILE}")
        return
    
    if not os.path.isfile(BAND_FILE):
        print(f"警告：找不到 {BAND_FILE}，将只绘制投影能带，不显示完整能带背景。")
        band_data = None
    else:
        print(f"读取完整能带文件: {BAND_FILE}")
        band_kpts, band_energies = read_band(BAND_FILE)
        band_data = (band_kpts, band_energies)
        print(f"  读取到 {len(band_kpts)} 个 k 点路径，每个路径包含 {band_energies[0].shape[0]} 条能带")
    
    # 读取高对称点
    sym_labels, sym_coords = read_klabels(KLABELS_FILE)
    print(f"高对称点: {list(zip(sym_labels, sym_coords))}")
    
    # 创建颜色映射
    cmap = create_colormap()
    print(f"使用颜色映射: 黄-橙-红渐变")
    
    # 处理每个 PBAND 文件
    for pband_file in PBAND_FILES:
        if not os.path.isfile(pband_file):
            print(f"警告：{pband_file} 不存在，跳过。")
            continue
        
        print(f"正在处理 {pband_file} ...")
        orbital_names, kpts, energies, weights = read_pband(pband_file)
        n_orbs = len(orbital_names)
        print(f"  轨道: {orbital_names}")
        
        # 计算全局最大权重
        global_max = 0.0
        for orb in orbital_names:
            max_w = np.max(weights[orb])
            if max_w > global_max:
                global_max = max_w
        
        print(f"  全局最大权重: {global_max:.4f}")
        
        # 计算子图布局（均匀排列）
        cols = int(np.ceil(np.sqrt(n_orbs)))
        rows = int(np.ceil(n_orbs / cols))
        
        # 创建子图
        fig, axes = plt.subplots(rows, cols, 
                                 figsize=(FIGURE_WIDTH, FIGURE_HEIGHT),
                                 sharex=True,
                                 sharey=True,
                                 gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        
        # 确保 axes 是数组
        axes = np.atleast_1d(axes).flatten()
        
        # 绘制所有轨道
        for idx, orb in enumerate(orbital_names):
            # 绘制完整能带背景
            if band_data is not None:
                plot_band_background(axes[idx], band_data[0], band_data[1], YLIM)
            
            # 绘制投影能带
            plot_fatband_color(axes[idx], kpts, energies, weights[orb],
                              ylim=YLIM,
                              high_sym_coords=sym_coords,
                              high_sym_labels=sym_labels,
                              title=orb,
                              global_max_weight=global_max,
                              cmap=cmap)
            
            # 设置坐标轴
            show_xlabel = (idx >= cols * (rows - 1))  # 只在最后一行显示 x 轴标签
            setup_axes(axes[idx], YLIM, sym_coords, sym_labels, show_xlabel)
        
        # 隐藏多余子图
        for idx in range(n_orbs, len(axes)):
            axes[idx].set_visible(False)
        
        # 为最后一个有内容的子图添加颜色条
        last_visible_idx = n_orbs - 1
        add_colorbar(fig, axes[last_visible_idx], cmap, global_max)
        
        # 设置总标题
        fig.suptitle(f"Projected Band Structure: {pband_file}",
                     fontsize=12, y=1.02)
        
        # 调整布局
        fig.tight_layout()
        
        # 保存图片
        if SAVE_FIG:
            base = os.path.splitext(os.path.basename(pband_file))[0]
            out_name = f"{OUTPUT_PREFIX}_{base}.{FIG_FORMAT}"
            plt.savefig(out_name, format=FIG_FORMAT, dpi=DPI, 
                       bbox_inches='tight', pad_inches=0.1)
            print(f"  图片已保存至 {out_name}")
            plt.close(fig)
        else:
            plt.show()

# ========================== 10. 单独绘制每个轨道 ==========================
def plot_each_orbital_separately():
    """
    为每个轨道单独生成一张图（可选功能）
    """
    # 检查输入文件
    if not os.path.isfile(KLABELS_FILE):
        print(f"错误：找不到 {KLABELS_FILE}")
        return
    
    if not os.path.isfile(BAND_FILE):
        print(f"警告：找不到 {BAND_FILE}，将只绘制投影能带，不显示完整能带背景。")
        band_data = None
    else:
        band_kpts, band_energies = read_band(BAND_FILE)
        band_data = (band_kpts, band_energies)
    
    sym_labels, sym_coords = read_klabels(KLABELS_FILE)
    
    # 处理每个 PBAND 文件
    for pband_file in PBAND_FILES:
        if not os.path.isfile(pband_file):
            continue
        
        print(f"正在处理 {pband_file} ...")
        orbital_names, kpts, energies, weights = read_pband(pband_file)
        
        # 计算全局最大权重
        global_max = 0.0
        for orb in orbital_names:
            max_w = np.max(weights[orb])
            if max_w > global_max:
                global_max = max_w
        
        # 为每个轨道单独生成图片
        for orb in orbital_names:
            fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
            
            # 创建颜色映射
            cmap = create_colormap()
            
            # 绘制完整能带背景
            if band_data is not None:
                plot_band_background(ax, band_data[0], band_data[1], YLIM)
            
            # 绘制单个轨道
            plot_fatband_color(ax, kpts, energies, weights[orb],
                              ylim=YLIM,
                              high_sym_coords=sym_coords,
                              high_sym_labels=sym_labels,
                              title=f"{orb}",
                              global_max_weight=global_max,
                              cmap=cmap)
            
            # 设置坐标轴
            setup_axes(ax, YLIM, sym_coords, sym_labels)
            
            # 添加颜色条
            add_colorbar(fig, ax, cmap, global_max)
            
            # 设置标题
            fig.suptitle(f"Projected Band Structure: {orb} orbital\n{pband_file}",
                        fontsize=10, y=1.02)
            
            # 调整布局
            fig.tight_layout()
            
            # 保存图片
            if SAVE_FIG:
                base = os.path.splitext(os.path.basename(pband_file))[0]
                out_name = f"{OUTPUT_PREFIX}_{base}_{orb}.{FIG_FORMAT}"
                plt.savefig(out_name, format=FIG_FORMAT, dpi=DPI,
                           bbox_inches='tight', pad_inches=0.1)
                print(f"  图片已保存至 {out_name}")
                plt.close(fig)
            else:
                plt.show()

if __name__ == "__main__":
    print("=" * 60)
    print("Fatband Color Plot Script (黄-橙-红渐变，带完整能带背景)")
    print("=" * 60)
    print(f"能量范围: {YLIM}")
    print(f"Y 轴刻度步长: {Y_TICK_STEP} eV")
    print(f"输出格式: {FIG_FORMAT}")
    print("=" * 60)
    
    # 1. 默认模式：均匀排列所有轨道
    print("\n1. 绘制均匀排列的子图...")
    main()
    
    # 2. 可选：单独绘制每个轨道（取消注释以启用）
    # print("\n2. 单独绘制每个轨道...")
    # plot_each_orbital_separately()
    
    print("\n" + "=" * 60)
    print("完成！")
    print("=" * 60)
