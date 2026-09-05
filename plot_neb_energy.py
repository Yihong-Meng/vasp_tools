#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NEB能量图绘制脚本
横坐标：归一化的反应坐标（0到1）
纵坐标：相对能量（meV）
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 设置期刊风格
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

def main():
    # NEB数据
    # 第一列：序号，第二列：最大受力，第三列：能量，第四列：相对能量(eV)
    data = np.array([
        [0,  0.000000, -479.624300, 0.000000],
        [1,  0.006897, -479.619600, 0.004700],
        [2,  0.011936, -479.618900, 0.005400],
        [3,  0.013341, -479.613200, 0.011100],
        [4,  0.015041, -479.604500, 0.019800],
        [5,  0.019190, -479.602900, 0.021400],
        [6,  0.025532, -479.599300, 0.025000],
        [7,  0.028779, -479.594900, 0.029400],
        [8,  0.032768, -479.599100, 0.025200],
        [9,  0.038516, -479.601800, 0.022500],
        [10, 0.039987, -479.603300, 0.021000],
        [11, 0.041062, -479.611800, 0.012500],
        [12, 0.051288, -479.616700, 0.007600],
        [13, 0.057520, -479.617400, 0.006900],
        [14, 0.000000, -479.622100, 0.002200]
    ])
    
    # 提取数据
    n_images = len(data)
    reaction_coordinate = data[:, 0] / (n_images - 1)  # 归一化到0-1
    energy_eV = data[:, 3]  # 相对能量(eV)
    energy_meV = energy_eV * 1000  # 转换为meV
    max_force = data[:, 1]  # 最大受力
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 绘制能量曲线
    ax.plot(reaction_coordinate, energy_meV, 'o-', color='blue', 
            markersize=8, linewidth=2, markerfacecolor='red', 
            markeredgecolor='black', markeredgewidth=1.5)
    
    # 找出过渡态（能量最高点）
    ts_idx = np.argmax(energy_meV)
    ts_energy = energy_meV[ts_idx]
    ts_coordinate = reaction_coordinate[ts_idx]
    
    # 标注过渡态
    ax.plot(ts_coordinate, ts_energy, 's', color='red', markersize=12, 
            markeredgecolor='black', markeredgewidth=2, label='Transition State')
    ax.annotate(f'TS: {ts_energy:.1f} meV', (ts_coordinate, ts_energy), 
               textcoords="offset points", xytext=(0, 15), ha='center', 
               fontsize=11, fontweight='bold')
    
    # 设置坐标轴
    ax.set_xlabel('Reaction Coordinate', fontsize=14)
    ax.set_ylabel('Relative Energy (meV)', fontsize=14)
    ax.set_title('NEB Energy Profile', fontsize=16, fontweight='bold')
    
    # 设置坐标轴范围
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(min(energy_meV) - 5, max(energy_meV) + 10)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # 添加图例
    ax.legend(loc='upper left', fontsize=11)
    
    # 设置坐标轴刻度
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    
    # 保存图片
    plt.savefig('neb_energy_profile.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    print("NEB能量图已保存至: neb_energy_profile.png")
    
    # 计算活化能
    activation_energy = ts_energy - energy_meV[0]
    
    # 打印关键信息
    print("\nNEB计算结果:")
    print(f"  图像数量: {n_images}")
    print(f"  起点能量: {energy_meV[0]:.1f} meV")
    print(f"  终点能量: {energy_meV[-1]:.1f} meV")
    print(f"  过渡态能量: {ts_energy:.1f} meV")
    print(f"  活化能: {activation_energy:.1f} meV")
    print(f"  反应能: {energy_meV[-1] - energy_meV[0]:.1f} meV")
    
    plt.close()

if __name__ == "__main__":
    main()
