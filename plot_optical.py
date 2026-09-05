#!/usr/bin/env python3
"""
Mn2BrI 光学性质分析与绘图
=========================
基于 postw90 输出的 kubo 光学导电率、JDOS 和 Wannier 能带，
绘制四张核心图：

  1) 介电函数 ε(ω) — 吸收边和介电响应
  2) Wannier 能带 vs DFT 能带 — 验证拟合在光学能量范围是否可靠
  3) JDOS + 光学导电率联合分析 — 跃迁归属
  4) 光学性质与可见光波段对比 — 显示吸收峰位于何种颜色

用法:
  python3 plot_optical.py              # 当前目录（必须有 kubo/jdos/band 文件）
  python3 plot_optical.py /path/to/dir

"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============== 用户配置 ==============
FERMI_ENERGY = -2.2883   # 费米能级 (eV)，用于能带平移
# =====================================

# 获取工作目录
work_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
os.chdir(work_dir)   # 切换到目标目录

print(f"Working directory: {os.getcwd()}")

# ====================================================================
# 1. 读 kubo 导电率数据
# ====================================================================
def read_kubo(filename):
    """
    读 wannier90-kubo_S_*.dat 格式:
    # omega  sigma_re  sigma_im   (可能还有其它列)
    返回 (omega, sigma_re, sigma_im)
    """
    if not os.path.exists(filename):
        return None, None, None
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 2:
                vals = [float(parts[0]), float(parts[1])]
                if len(parts) >= 3:
                    vals.append(float(parts[2]))
                data.append(vals)
    if not data:
        return None, None, None
    data = np.array(data)
    omega = data[:,0]
    sigma_re = data[:,1] if data.shape[1] > 1 else None
    sigma_im = data[:,2] if data.shape[1] > 2 else None
    return omega, sigma_re, sigma_im

# 读取 xx 分量（面内）
omega, sxx_re, sxx_im = read_kubo('wannier90-kubo_S_xx.dat')
if omega is None:
    print("ERROR: 找不到 wannier90-kubo_S_xx.dat 或数据格式错误")
    sys.exit(1)

print(f"  ω 范围: [{omega.min():.2f}, {omega.max():.2f}] eV, {len(omega)} 个点")
print(f"  σ_xx 峰值: {np.max(np.abs(sxx_re)):.4f}")

# 可选读取 yy 分量
omega_yy, syy_re, _ = read_kubo('wannier90-kubo_S_yy.dat')

# ====================================================================
# 2. 读 JDOS
# ====================================================================
jdos = None
jdos_file = 'wannier90-jdos.dat'
if os.path.exists(jdos_file):
    jdos = []
    with open(jdos_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 2:
                jdos.append([float(parts[0]), float(parts[1])])
    jdos = np.array(jdos)
    print(f"  JDOS 范围: [{jdos[:,0].min():.2f}, {jdos[:,0].max():.2f}] eV")
else:
    print("  警告: 未找到 JDOS 文件，将跳过相关绘图")

# ====================================================================
# 3. 图1：介电函数 ε(ω)
# ====================================================================
fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(14, 5))

ax1a.plot(omega, sxx_re, 'b-', lw=1.5, label=r'$\mathrm{Re}\,\sigma_{xx}$')
if syy_re is not None:
    ax1a.plot(omega_yy, syy_re, 'g--', lw=1.2, alpha=0.7, label=r'$\mathrm{Re}\,\sigma_{yy}$')
ax1a.axvline(x=0, color='gray', ls='--', lw=0.5)
ax1a.set_xlabel(r'$\omega$ (eV)', fontsize=12)
ax1a.set_ylabel(r'$\mathrm{Re}\,\sigma$', fontsize=12)
ax1a.set_title('Optical conductivity', fontsize=13)
ax1a.legend(fontsize=10)
ax1a.grid(True, alpha=0.3)

# ε₂(ω) = 4π σ_xx / ω  (CGS 单位)
eps2 = np.zeros_like(sxx_re)
mask = np.abs(omega) > 0.01
eps2[mask] = 4 * np.pi * sxx_re[mask] / omega[mask]

ax1b.plot(omega[mask], eps2[mask], 'r-', lw=1.5, label=r'$\varepsilon_2(\omega)$')
ax1b.axvline(x=0, color='gray', ls='--', lw=0.5)
ax1b.set_xlabel(r'$\omega$ (eV)', fontsize=12)
ax1b.set_ylabel(r'$\varepsilon_2$', fontsize=12)
ax1b.set_title(r'Dielectric function $\varepsilon_2$', fontsize=13)
ax1b.legend(fontsize=10)
ax1b.grid(True, alpha=0.3)

plt.tight_layout()
fig1.savefig('optical_epsilon.png', dpi=150)
fig1.savefig('optical_epsilon.pdf', dpi=150)
print(f"\n  已保存: optical_epsilon.png/pdf")

# ====================================================================
# 4. 图2：能带对比 (Wannier vs DFT) 带费米能级平移
# ====================================================================
def read_bands(filename, shift_fermi=True):
    """
    读能带文件，返回每条能带的 (k, E) 数组。
    若 shift_fermi=True，则 E 减去 FERMI_ENERGY。
    支持 BAND.dat 格式（多段 # Band-Index）和 wannier90_band.dat 格式。
    """
    if not os.path.exists(filename):
        return []
    bands = []
    cur_band = []
    with open(filename) as f:
        for line in f:
            line_strip = line.strip()
            if not line_strip or line_strip.startswith('#'):
                if cur_band:
                    bands.append(np.array(cur_band))
                    cur_band = []
                continue
            parts = line_strip.split()
            if len(parts) >= 2:
                try:
                    k = float(parts[0])
                    e = float(parts[1])
                except ValueError:
                    continue
                if shift_fermi:
                    e -= FERMI_ENERGY
                cur_band.append([k, e])
        if cur_band:
            bands.append(np.array(cur_band))
    return bands

wannier_file = 'wannier90_band.dat'
dft_file = 'BAND.dat'

if os.path.exists(wannier_file) and os.path.exists(dft_file):
    w_bands = read_bands(wannier_file, shift_fermi=True)
    d_bands = read_bands(dft_file, shift_fermi=True)

    fig2, ax2 = plt.subplots(figsize=(8, 6))

    # DFT 能带（灰色）
    for b in d_bands:
        ax2.plot(b[:,0], b[:,1], 'gray', lw=0.5, alpha=0.5)

    # Wannier 能带（蓝色）
    for b in w_bands:
        ax2.plot(b[:,0], b[:,1], 'b-', lw=1.0, alpha=0.8)

    # 费米能级参考线
    ax2.axhline(y=0, color='black', ls='--', lw=0.8, alpha=0.7, label='E_F')

    ax2.set_xlabel('k-path', fontsize=12)
    ax2.set_ylabel('Energy - E_F (eV)', fontsize=12)
    ax2.set_title('Wannier (blue) vs DFT (gray) band structure', fontsize=13)
    ax2.legend(handles=[
        Line2D([0],[0], color='blue', lw=1.5, label='Wannier'),
        Line2D([0],[0], color='gray', lw=1.5, label='DFT'),
        Line2D([0],[0], color='black', ls='--', lw=1, label='E_F')
    ], fontsize=10)
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    fig2.savefig('band_comparison_new.png', dpi=150)
    fig2.savefig('band_comparison_new.pdf', dpi=150)
    print(f"  已保存: band_comparison_new.png/pdf")
else:
    print(f"  警告: 找不到 band 数据文件，跳过能带对比")

# ====================================================================
# 5. 图3：JDOS + 光学导电率联合分析
# ====================================================================
if jdos is not None:
    fig3, ax3 = plt.subplots(figsize=(10, 5))

    ax3.plot(omega, sxx_re, 'b-', lw=1.5, label=r'$\mathrm{Re}\,\sigma_{xx}(\omega)$')
    ax3.set_xlabel(r'$\omega$ (eV)', fontsize=12)
    ax3.set_ylabel(r'$\mathrm{Re}\,\sigma_{xx}$', fontsize=12, color='blue')
    ax3.tick_params(axis='y', labelcolor='blue')

    ax3b = ax3.twinx()
    ax3b.plot(jdos[:,0], jdos[:,1], 'r--', lw=1.2, alpha=0.7, label='JDOS')
    ax3b.set_ylabel('JDOS', fontsize=12, color='red')
    ax3b.tick_params(axis='y', labelcolor='red')

    ax3.set_title('Optical conductivity + Joint DOS', fontsize=13)
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3b.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    fig3.savefig('sigma_jdos.png', dpi=150)
    fig3.savefig('sigma_jdos.pdf', dpi=150)
    print(f"  已保存: sigma_jdos.png/pdf")

# ====================================================================
# 6. 图4：光学性质 + 可见光波段标识（修正版）
# ====================================================================
def add_visible_spectrum(ax, alpha=0.15):
    """
    在坐标轴 ax 上添加可见光波段（1.65-3.26 eV）的彩虹色背景。
    使用 axvspan 在 x 方向绘制垂直色带。
    """
    # 能量边界 (eV) 对应颜色
    energy_bounds = [1.65, 2.00, 2.10, 2.18, 2.50, 2.76, 3.26]
    colors = ['orange', 'yellow', 'limegreen', 'cyan', 'blue', 'violet']
    for i in range(len(energy_bounds)-1):
        xmin, xmax = energy_bounds[i], energy_bounds[i+1]
        ax.axvspan(xmin, xmax, facecolor=colors[i], alpha=alpha, edgecolor='none')
    # 添加标记文本（放在左下角或顶部）
    ylim = ax.get_ylim()
    ypos = ylim[0] + 0.05 * (ylim[1] - ylim[0])  # 放在底部上方一点点
    ax.text(1.65, ypos, 'Visible', ha='left', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

# 创建新图
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 5))

# 子图1：σ_xx 与可见光
ax4a.plot(omega, sxx_re, 'b-', lw=2, label=r'$\mathrm{Re}\,\sigma_{xx}$')
if syy_re is not None:
    ax4a.plot(omega_yy, syy_re, 'g--', lw=1.5, alpha=0.7, label=r'$\mathrm{Re}\,\sigma_{yy}$')
ax4a.axvline(x=0, color='gray', ls='--', lw=0.5)
ax4a.set_xlabel(r'$\omega$ (eV)', fontsize=12)
ax4a.set_ylabel(r'$\mathrm{Re}\,\sigma$', fontsize=12)
ax4a.set_title('Optical conductivity with visible spectrum', fontsize=13)
ax4a.legend(fontsize=10)
ax4a.grid(True, alpha=0.3)
# 添加可见光色带
add_visible_spectrum(ax4a, alpha=0.2)
# 调整x轴范围以包含可见光区域
ax4a.set_xlim(left=min(0, omega.min()), right=max(4.0, omega.max()))

# 子图2：ε₂ 与可见光
mask2 = mask  # 避免除以零
ax4b.plot(omega[mask2], eps2[mask2], 'r-', lw=2, label=r'$\varepsilon_2(\omega)$')
ax4b.axvline(x=0, color='gray', ls='--', lw=0.5)
ax4b.set_xlabel(r'$\omega$ (eV)', fontsize=12)
ax4b.set_ylabel(r'$\varepsilon_2$', fontsize=12)
ax4b.set_title(r'Dielectric function $\varepsilon_2$ with visible spectrum', fontsize=13)
ax4b.legend(fontsize=10)
ax4b.grid(True, alpha=0.3)
add_visible_spectrum(ax4b, alpha=0.2)
ax4b.set_xlim(left=min(0, omega.min()), right=max(4.0, omega.max()))

plt.tight_layout()
fig4.savefig('optical_visible.png', dpi=150)
fig4.savefig('optical_visible.pdf', dpi=150)
print(f"  已保存: optical_visible.png/pdf")

# ====================================================================
# 汇总分析提示
# ====================================================================
print("\n===== 分析要点 =====")
print("1. ε₂(ω) 出现第一个峰/台阶 → 对应带间跃迁阈值（光学带隙）")
print("2. 峰位 → 对应 k 空间中贡献最大的跃迁（查能带图定位）")
print("3. JDOS 峰位置应和 σ_xx 峰基本一致（但不完全重合，因为矩阵元权重不同）")
print("4. 能带对比：Wannier 在光学能量范围（0~6 eV）内必须和 DFT 一致")
print("5. σ_yy ≈ σ_xx 表示面内各向同性；σ_xy ≈ 0 表示无反常霍尔效应")
print("6. 费米能级 E_F = {:.4f} eV 已作为能量零点".format(FERMI_ENERGY))
print("7. 可见光范围 (380-750 nm) 已用彩虹色条标注在 optical_visible 图中")
