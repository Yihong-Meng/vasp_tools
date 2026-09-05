#!/usr/bin/env python3
"""
DFT 能带与 Wannier90 插值能带对比绘图（单自旋通道）
自动适配 Wannier90 的两种输出格式：
1. 多列格式：k, band1, band2, ...
2. 两列顺序格式：k, e (所有能带首尾相连)

能量自动减去费米能级 (FERMI_ENERGY)，使费米能级为 0 eV。
可手动指定 Y 轴范围（YLIM），方便与 Fatband 图对比。
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os

# ========================== 用户配置 ==========================
SPIN_CHANNEL = 'down'               # 选择 'up' 或 'down'
FERMI_ENERGY = 0.2785             # 费米能级 (eV)，用于将能量平移至 E_F=0

# 手动指定 Y 轴范围，例如 YLIM = (-2, 2) ；若设为 None 则自动确定
YLIM = (-2, 2)                    # 自定义能量范围，与 Fatband 图保持一致

OUTPUT_IMAGE = 'band_comparison.png'
DPI = 300

DFT_COLOR = 'black'
DFT_LINEWIDTH = 2.0
DFT_LINESTYLE = '-'

WANNIER_COLOR = 'red'
WANNIER_LINEWIDTH = 0.8
WANNIER_LINESTYLE = '-'
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
        raise ValueError("无法解析 NKPTS/NBANDS")

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
    智能读取 wannier90_band.dat：
    - 如果列数 > 2：按多列格式解析 (k, band1, band2, ...)
    - 如果列数 == 2：按顺序格式解析，自动检测 k 点重置位置并分割能带
    返回的能带能量已减去 FERMI_ENERGY。
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

        print(f"  检测到两列顺序格式，自动分割为 {energies.shape[1]} 条能带")
        return k, energies

    else:
        raise ValueError("无法解析 wannier90_band.dat 格式")

def parse_gnu_xtics(gnu_file):
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

def parse_gnu_range(gnu_file):
    with open(gnu_file, 'r') as f:
        content = f.read()
    xmatch = re.search(r'set\s+xrange\s*\[\s*([\d.-]+)\s*:\s*([\d.-]+)\s*\]', content)
    ymatch = re.search(r'set\s+yrange\s*\[\s*([\d.-]+)\s*:\s*([\d.-]+)\s*\]', content)
    xmin = xmax = ymin = ymax = None
    if xmatch:
        xmin, xmax = float(xmatch.group(1)), float(xmatch.group(2))
    if ymatch:
        ymin, ymax = float(ymatch.group(1)), float(ymatch.group(2))
    return xmin, xmax, ymin, ymax

def main():
    for f in ['BAND.dat', 'wannier90_band.dat']:
        if not os.path.isfile(f):
            print(f"错误：{f} 不存在")
            return

    spin = SPIN_CHANNEL.lower()
    if spin not in ['up', 'down']:
        print("SPIN_CHANNEL 必须是 'up' 或 'down'")
        return

    print(f"费米能级: {FERMI_ENERGY:.4f} eV，所有能量已平移至 E_F=0")

    # DFT
    k_dft_list, e_dft_list = read_BAND('BAND.dat', spin)
    print(f"DFT 能带数：{len(k_dft_list)}")

    # Wannier
    k_wann, e_wann = read_wannier_band('wannier90_band.dat')
    print(f"Wannier 能带数：{e_wann.shape[1]}")

    # 解析 gnu 信息
    gnu_file = 'wannier90_band.gnu'
    labels, positions = [], []
    xmin = xmax = ymin = ymax = None
    if os.path.isfile(gnu_file):
        try:
            labels, positions = parse_gnu_xtics(gnu_file)
            print(f"高对称点：{list(zip(labels, positions))}")
        except:
            pass
        try:
            xmin, xmax, ymin, ymax = parse_gnu_range(gnu_file)
        except:
            pass

    fig, ax = plt.subplots(figsize=(8, 6))

    # 绘制 DFT
    for i, (kd, ed) in enumerate(zip(k_dft_list, e_dft_list)):
        ax.plot(kd, ed, color=DFT_COLOR, lw=DFT_LINEWIDTH, ls=DFT_LINESTYLE,
                label='DFT' if i == 0 else "")

    # 绘制 Wannier
    for i in range(e_wann.shape[1]):
        ax.plot(k_wann, e_wann[:, i], color=WANNIER_COLOR,
                lw=WANNIER_LINEWIDTH, ls=WANNIER_LINESTYLE,
                label='Wannier' if i == 0 else "")

    # 费米能级参考线
    ax.axhline(y=0, color='gray', ls='--', lw=0.5, alpha=0.7)

    # ----- 设置坐标轴范围 -----
    # X 轴：优先使用 gnu 中的 xrange，否则自动
    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    else:
        all_k = np.concatenate([kd for kd in k_dft_list] + [k_wann])
        ax.set_xlim(np.min(all_k), np.max(all_k))

    # Y 轴：如果用户设置了 YLIM，则强制使用；否则从 gnu 读取 yrange；否则自动
    if YLIM is not None:
        ax.set_ylim(YLIM[0], YLIM[1])
        print(f"使用手动 Y 轴范围：{YLIM[0]:.2f} ~ {YLIM[1]:.2f} eV")
    elif ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)
        print(f"从 gnu 文件读取 Y 轴范围：{ymin:.2f} ~ {ymax:.2f} eV")
    else:
        all_e = np.concatenate([ed for ed in e_dft_list] + [e_wann.flatten()])
        margin = 0.1 * (np.max(all_e) - np.min(all_e))
        ax.set_ylim(np.min(all_e) - margin, np.max(all_e) + margin)
        print("自动确定 Y 轴范围")

    # 高对称点刻度
    if labels and positions:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        for pos in positions:
            ax.axvline(x=pos, color='gray', ls='--', lw=0.5, alpha=0.7)

    ax.set_xlabel(r'$k$-path')
    ax.set_ylabel('Energy - E_F (eV)')
    ax.set_title(f'DFT vs Wannier (Spin-{spin.upper()})')
    ax.legend()
    ax.grid(True, alpha=0.2)

    plt.savefig(OUTPUT_IMAGE, dpi=DPI, bbox_inches='tight')
    print(f"图片已保存至 {OUTPUT_IMAGE}")
    plt.close()

if __name__ == "__main__":
    main()
