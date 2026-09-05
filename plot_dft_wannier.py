#!/usr/bin/env python3
"""
同时绘制自旋向上和自旋向下的 DFT 与 Wannier 能带对比图。
- DFT：实线，不透明，加粗
- Wannier：虚线，半透明，较细
颜色：红色（up），蓝色（dn）
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os

# ========================== 用户配置 ==========================
# -------- 数据路径 ----------
UP_DIR = "./up"                 # 自旋向上文件夹路径
DN_DIR = "./dn"                 # 自旋向下文件夹路径

# -------- 能量参数 ----------
FERMI_ENERGY = -2.2883          # 费米能级 (eV)，用于平移能量
YLIM = (-3.5, 3)                  # Y 轴范围 (eV)，设为 None 则自动

# -------- 输出设置 ----------
OUTPUT_IMAGE = "dft_vs_wannier_both_spins.png"
DPI = 300

# -------- 绘图样式（增强对比） ----------
# DFT 样式
DFT_UP_COLOR = 'red'
DFT_DN_COLOR = 'blue'
DFT_LINESTYLE = '-'             # 实线
DFT_LINEWIDTH = 2.5             # 加粗
DFT_ALPHA = 0.3                 # 不透明

# Wannier 样式
WANNIER_UP_COLOR = 'red'
WANNIER_DN_COLOR = 'blue'
WANNIER_LINESTYLE = '--'        # 虚线
WANNIER_LINEWIDTH = 1.0         # 较细
WANNIER_ALPHA = 0.9            # 半透明

# =============================================================

def read_BAND(filename, spin='up'):
    """从 BAND.dat 读取 DFT 能带（按 # Band-Index 分割）"""
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
            e_list.append(e_arr[sort_idx] - FERMI_ENERGY)
    return k_list, e_list

def read_wannier_band(filename):
    """
    读取 wannier90_band.dat，支持两种格式：
    1. 多列格式 (k, e1, e2, ...)
    2. 两列顺序格式 (k, e) 首尾相连
    返回 (k, energies_matrix)，矩阵每列是一条能带。
    """
    data = np.loadtxt(filename)

    if data.ndim == 2 and data.shape[1] > 2:
        k = data[:, 0]
        energies = data[:, 1:] - FERMI_ENERGY
        return k, energies

    elif data.ndim == 2 and data.shape[1] == 2:
        k_all = data[:, 0]
        e_all = data[:, 1]

        split_idx = np.where(np.diff(k_all) < 0)[0] + 1
        if len(split_idx) == 0:
            return k_all, (e_all - FERMI_ENERGY).reshape(-1, 1)

        k_segments = []
        e_segments = []
        start = 0
        for end in split_idx:
            k_segments.append(k_all[start:end])
            e_segments.append(e_all[start:end])
            start = end
        k_segments.append(k_all[start:])
        e_segments.append(e_all[start:])

        min_len = min(len(k) for k in k_segments)
        energies = np.column_stack([e[:min_len] - FERMI_ENERGY for e in e_segments])
        k = k_segments[0][:min_len]
        return k, energies
    else:
        raise ValueError(f"无法解析 {filename} 格式")

def parse_gnu_xtics(gnu_file):
    """从 gnu 文件中解析 set xtics"""
    with open(gnu_file, 'r') as f:
        content = f.read()
    pattern = r'set\s+xtics\s*\(\s*(.*?)\s*\)'
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        raise ValueError("未找到 xtics")
    tics_str = match.group(1)
    pairs = re.findall(r'\"([^\"]*)\"\s*:\s*([\d.]+)', tics_str)
    if not pairs:
        pairs = re.findall(r'\"([^\"]*)\"\s+([\d.]+)', tics_str)
    labels = [p[0] for p in pairs]
    positions = [float(p[1]) for p in pairs]
    return labels, positions

def parse_gnu_xrange(gnu_file):
    """从 gnu 文件中解析 xrange"""
    with open(gnu_file, 'r') as f:
        content = f.read()
    match = re.search(r'set\s+xrange\s*\[\s*([\d.-]+)\s*:\s*([\d.-]+)\s*\]', content)
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None

def main():
    # 检查文件夹是否存在
    for d in [UP_DIR, DN_DIR]:
        if not os.path.isdir(d):
            print(f"错误：目录 {d} 不存在")
            return

    # 读取 DFT 数据
    print("读取 DFT 能带...")
    k_up_dft, e_up_dft = read_BAND(os.path.join(UP_DIR, "BAND.dat"), spin='up')
    k_dn_dft, e_dn_dft = read_BAND(os.path.join(DN_DIR, "BAND.dat"), spin='down')
    print(f"  Up DFT: {len(k_up_dft)} 条能带")
    print(f"  Dn DFT: {len(k_dn_dft)} 条能带")

    # 读取 Wannier 数据
    print("读取 Wannier 插值能带...")
    k_up_wann, e_up_wann = read_wannier_band(os.path.join(UP_DIR, "wannier90_band.dat"))
    k_dn_wann, e_dn_wann = read_wannier_band(os.path.join(DN_DIR, "wannier90_band.dat"))
    print(f"  Up Wannier: {e_up_wann.shape[1]} 条能带")
    print(f"  Dn Wannier: {e_dn_wann.shape[1]} 条能带")

    # 从 up 的 gnu 文件解析高对称点（如果存在）
    gnu_file_up = os.path.join(UP_DIR, "wannier90_band.gnu")
    labels, positions = [], []
    if os.path.isfile(gnu_file_up):
        try:
            labels, positions = parse_gnu_xtics(gnu_file_up)
            print(f"高对称点：{list(zip(labels, positions))}")
        except Exception as e:
            print(f"解析 xtics 失败：{e}")
            labels, positions = [], []
    else:
        # 尝试从 dn 读取
        gnu_file_dn = os.path.join(DN_DIR, "wannier90_band.gnu")
        if os.path.isfile(gnu_file_dn):
            try:
                labels, positions = parse_gnu_xtics(gnu_file_dn)
                print(f"高对称点：{list(zip(labels, positions))}")
            except:
                pass

    # 获取 xrange（从 up 的 gnu 读取，若没有则自动）
    xmin, xmax = None, None
    if os.path.isfile(gnu_file_up):
        xmin, xmax = parse_gnu_xrange(gnu_file_up)

    # 绘图
    fig, ax = plt.subplots(figsize=(9, 6))

    # ---- 绘制 DFT 能带 (实线，不透明，加粗) ----
    for i, (k, e) in enumerate(zip(k_up_dft, e_up_dft)):
        ax.plot(k, e, color=DFT_UP_COLOR, linestyle=DFT_LINESTYLE,
                linewidth=DFT_LINEWIDTH, alpha=DFT_ALPHA,
                label='DFT up' if i == 0 else "")

    for i, (k, e) in enumerate(zip(k_dn_dft, e_dn_dft)):
        ax.plot(k, e, color=DFT_DN_COLOR, linestyle=DFT_LINESTYLE,
                linewidth=DFT_LINEWIDTH, alpha=DFT_ALPHA,
                label='DFT dn' if i == 0 else "")

    # ---- 绘制 Wannier 能带 (虚线，半透明，较细) ----
    for i in range(e_up_wann.shape[1]):
        ax.plot(k_up_wann, e_up_wann[:, i], color=WANNIER_UP_COLOR,
                linestyle=WANNIER_LINESTYLE, linewidth=WANNIER_LINEWIDTH,
                alpha=WANNIER_ALPHA, label='Wannier up' if i == 0 else "")

    for i in range(e_dn_wann.shape[1]):
        ax.plot(k_dn_wann, e_dn_wann[:, i], color=WANNIER_DN_COLOR,
                linestyle=WANNIER_LINESTYLE, linewidth=WANNIER_LINEWIDTH,
                alpha=WANNIER_ALPHA, label='Wannier dn' if i == 0 else "")

    # 费米能级参考线
    ax.axhline(y=0, color='black', linestyle=':', linewidth=0.8, alpha=0.7)

    # 坐标轴范围
    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    else:
        # 自动：合并所有 k 点
        all_k = np.concatenate([k_up_dft[0], k_dn_dft[0], k_up_wann, k_dn_wann])
        ax.set_xlim(np.min(all_k), np.max(all_k))

    if YLIM is not None:
        ax.set_ylim(YLIM[0], YLIM[1])
        print(f"使用手动 Y 轴范围：{YLIM[0]:.2f} ~ {YLIM[1]:.2f} eV")
    else:
        # 自动
        all_e = np.concatenate([e.flatten() for e in e_up_dft] +
                               [e.flatten() for e in e_dn_dft] +
                               [e_up_wann.flatten(), e_dn_wann.flatten()])
        margin = 0.1 * (np.max(all_e) - np.min(all_e))
        ax.set_ylim(np.min(all_e) - margin, np.max(all_e) + margin)
        print("自动确定 Y 轴范围")

    # 高对称点垂直线和刻度
    if labels and positions:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        for pos in positions:
            ax.axvline(x=pos, color='gray', linestyle='--', linewidth=0.5, alpha=0.6)

    ax.set_xlabel(r'$k$-path')
    ax.set_ylabel('Energy - E_F (eV)')
    ax.set_title('DFT vs Wannier: Both Spins')

    # 图例（只显示一次，避免重复）
    handles, labels_leg = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_leg, handles))  # 去重
    ax.legend(by_label.values(), by_label.keys(), loc='best')

    ax.grid(True, alpha=0.15)

    plt.savefig(OUTPUT_IMAGE, dpi=DPI, bbox_inches='tight')
    print(f"图片已保存至 {OUTPUT_IMAGE}")
    plt.close()

if __name__ == "__main__":
    main()
