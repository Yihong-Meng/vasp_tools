#!/usr/bin/env python3
"""
分析 VASP EIGENVAL 文件，提取能带能量信息
"""
import numpy as np
import sys

def read_eigenval(filename):
    """读取 EIGENVAL 文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 解析头文件（前 5 行）
    # 第 1 行：体系信息
    # 第 2 行：nkpts, nbands, 等
    # 第 3-5 行：晶格常数等
    
    # 查找 nkpts 和 nbands
    for i in range(10):
        parts = lines[i].split()
        if len(parts) >= 2:
            try:
                nkpts = int(parts[0])
                nbands = int(parts[1])
                if 1 <= nkpts <= 500000 and 1 <= nbands <= 50000:
                    header_lines = i
                    break
            except:
                continue
    
    print(f"nkpts = {nkpts}")
    print(f"nbands = {nbands}")
    
    # 解析能量数据
    energies = []
    i = header_lines + 1
    
    while i < len(lines) and len(energies) < nkpts:
        # 跳过空行
        if not lines[i].strip():
            i += 1
            continue
        
        parts = lines[i].split()
        if len(parts) == 4:  # k 点信息行
            k_energies = []
            i += 1
            # 读取能带能量
            while i < len(lines) and len(k_energies) < nbands:
                if not lines[i].strip():
                    i += 1
                    continue
                band_parts = lines[i].split()
                if len(band_parts) >= 2:
                    k_energies.append(float(band_parts[1]))  # 第二列是能量
                i += 1
                # 检查是否到下一个 k 点
                if len(band_parts) == 4 and len(k_energies) > 0:
                    break
            
            if len(k_energies) == nbands:
                energies.append(k_energies)
    
    return np.array(energies), nbands, nkpts

def analyze_eigenval(filename, fermi_energy=0.0):
    """分析 EIGENVAL 文件"""
    print(f"分析文件: {filename}")
    print("=" * 50)
    
    energies, nbands, nkpts = read_eigenval(filename)
    
    print(f"\n总能带数: {nbands}")
    print(f"总 k 点数: {nkpts}")
    
    # 计算费米能级（如果有）
    # VASP 的 EIGENVAL 中费米能级通常在头文件中
    # 这里假设已经从 OUTCAR 获取并传入
    
    print(f"\n费米能级: {fermi_energy:.4f} eV")
    
    # 分析每条能带
    print(f"\n{'能带':>6} | {'最小能量':>12} | {'最大能量':>12} | {'宽度':>10} | {'费米能级以下':>12}")
    print("-" * 70)
    
    band_analysis = []
    for band_idx in range(nbands):
        band_energies = energies[:, band_idx]
        emin = np.min(band_energies)
        emax = np.max(band_energies)
        width = emax - emin
        below_fermi = np.sum(band_energies < fermi_energy)
        
        band_analysis.append({
            'band': band_idx + 1,
            'emin': emin,
            'emax': emax,
            'width': width,
            'below_fermi': below_fermi
        })
        
        # 只打印前 60 条能带
        if band_idx < 60:
            print(f"{band_idx+1:6d} | {emin:12.4f} | {emax:12.4f} | {width:10.4f} | {below_fermi:12d}")
    
    # 统计费米能级附近的能带
    print(f"\n{'='*50}")
    print("费米能级附近的能带统计:")
    
    energy_ranges = [
        (-1.0, 1.0),
        (-2.0, 2.0),
        (-5.0, 5.0),
        (-10.0, 10.0)
    ]
    
    for emin, emax in energy_ranges:
        count = 0
        for band_idx in range(nbands):
            band_energies = energies[:, band_idx]
            band_center = np.mean(band_energies)
            if emin <= band_center <= emax:
                count += 1
        print(f"  {emin:6.1f} ~ {emax:6.1f} eV: {count} 条能带")
    
    # 能量范围统计
    all_energies = energies.flatten()
    print(f"\n整体能量范围: {np.min(all_energies):.4f} ~ {np.max(all_energies):.4f} eV")
    print(f"费米能级以下总能带数: {np.sum(all_energies < fermi_energy)}")
    
    return band_analysis

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python3 analyze_eigenval.py EIGENVAL [费米能级]")
        print("示例: python3 analyze_eigenval.py EIGENVAL -3.5")
        sys.exit(1)
    
    filename = sys.argv[1]
    fermi_energy = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    
    analyze_eigenval(filename, fermi_energy)
