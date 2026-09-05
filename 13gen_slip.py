#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_slip.py — 生成 8×8 = 64 个滑移结构 POSCAR
=================================================

结构（68 原子，自下而上）:
  底部 Mn2BrI (8 Mn + 4 Br + 4 I = 16 原子)  ← 固定
  BN 双层    (B9N9 + B9N9 = 36 原子)          ← 固定
  顶部 Mn2BrI (8 Mn + 4 Br + 4 I = 16 原子)  ← 滑移

步骤：
  1) 自动识别 Z 轴（堆叠方向）：
     三条晶格矢量中模长最大者 = 层间堆叠矢量（本例为 c 轴）；
     另外两条为面内晶格方向 a1、a2，即滑移方向。
  2) 自动识别最上层 Mn2BrI：
     沿堆叠方向分数坐标排序，找层间最大间隙，取最上层并校验成分为
     Mn2BrI (8 Mn + 4 Br + 4 I)。若自动识别失败可用 --zthresh 手动指定。
  3) 沿 a1、a2 两个面内方向分别以 1/7 为步长滑移最上层（含 0，共 8×8 = 64 个）。
  4) 所有 68 个原子 Selective Dynamics 标记均为 F F T
     （x、y 固定不动，只优化 z 轴）。

输出：
  structures/POSCAR_{y}_{x}   （x = 沿 a1 的步数, y = 沿 a2 的步数）

用法：
  python3 gen_slip.py                 # 使用当前目录 POSCAR
  python3 gen_slip.py -p POSCAR_7.55  # 指定输入 POSCAR
  python3 gen_slip.py -d 7            # 晶格等分数（步长=1/7），默认 7
  python3 gen_slip.py --zthresh 0.62  # 手动指定顶层 z 阈值
"""

import os
import sys
import argparse
from collections import Counter

import numpy as np

# ================= 参数设置 =================
DEFAULT_NDIV = 7          # 晶格等分数：步长 = a/NDIV（默认1/7）
                          # 每方向取 0,1/7,...,7/7 共 NDIV+1 = 8 个点 → 8×8=64
TOP_EXPECT   = {'Mn': 8, 'Br': 4, 'I': 4}   # 顶层 Mn2BrI 期望成分
# ===========================================


def read_poscar(poscar_path):
    """读取 POSCAR（支持 Direct/Cartesian、Selective Dynamics），自动跳过末尾零行。"""
    with open(poscar_path, 'r') as f:
        lines = f.readlines()

    comment = lines[0].strip()
    scale = float(lines[1].strip())
    lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
    lattice *= scale

    # 元素行与数量行
    line5 = lines[5].split()
    line6 = lines[6].split()
    try:
        atom_counts = list(map(int, line5))
        elements = None
        coord_start_line = 6
    except ValueError:
        elements = line5
        atom_counts = list(map(int, line6))
        coord_start_line = 7

    # 跳过 Selective Dynamics 标记行
    if 'selective' in lines[coord_start_line - 1].lower():
        coord_start_line += 1

    expected = sum(atom_counts)
    coords = []
    for line in lines[coord_start_line:]:
        if line.strip() == '':
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        except ValueError:
            continue
        # 跳过 CONTCAR 末尾填充的 0 行
        if abs(x) < 1e-10 and abs(y) < 1e-10 and abs(z) < 1e-10:
            continue
        coords.append([x, y, z])
        if len(coords) == expected:
            break

    coords = np.array(coords)
    if len(coords) != expected:
        raise ValueError(f"期望 {expected} 个原子，读取到 {len(coords)} 个")

    # Cartesian → 分数坐标
    is_cart = lines[coord_start_line - 1].strip().lower().startswith(('c', 'k'))
    if is_cart:
        inv_lattice = np.linalg.inv(lattice)
        coords = coords @ inv_lattice.T

    # 原子种类列表
    atom_types = []
    if elements is not None:
        for el, cnt in zip(elements, atom_counts):
            atom_types.extend([el] * cnt)
    else:
        atom_types = ['Atom'] * len(coords)

    return lattice, elements, atom_counts, coords, atom_types, comment


def detect_stacking_axis(lattice):
    """
    返回堆叠方向对应的晶格矢量索引 (0/1/2)。
    判据：模长最大的矢量即为层间堆叠矢量；同时检查它是否近似垂直于另外两条。
    """
    norms = np.linalg.norm(lattice, axis=1)
    k = int(np.argmax(norms))
    others = [i for i in range(3) if i != k]
    # 垂直性检查（点积与特征长度乘积之比 < 5%）
    for o in others:
        cos_angle = abs(np.dot(lattice[k], lattice[o])) / (norms[k] * norms[o])
        if cos_angle > 0.05:
            print(f"  [警告] 晶格矢量{k+1}与矢量{o+1}不严格垂直 (|cosθ|={cos_angle:.3f})")
            print("         请人工确认 Z 轴方向是否正确。")
    return k


def detect_top_layer(coords, atom_types, col, expected=None, zthresh=None):
    """
    找出最上层（滑移层）原子索引。
    - 若给定 zthresh，直接取 coords[:, col] > zthresh。
    - 否则按 col 列排序，从最大的层间间隙开始尝试，取“间隙上方原子”
      且成分符合 expected 的那一层。
    返回 (indices, boundary_z)
    """
    n = len(coords)
    order = np.argsort(coords[:, col])
    sorted_z = coords[order, col]

    if zthresh is not None:
        idx = np.where(coords[:, col] > zthresh)[0]
        return idx, zthresh

    # 相邻原子之间的间隙
    gaps = sorted_z[1:] - sorted_z[:-1]
    for gi in np.argsort(gaps)[::-1]:           # 从最大间隙开始尝试
        boundary = 0.5 * (sorted_z[gi] + sorted_z[gi + 1])
        top_idx = np.where(coords[:, col] > boundary)[0]
        if len(top_idx) == 0:
            continue
        if expected is not None:
            cnt = Counter(atom_types[i] for i in top_idx)
            # 成分完全匹配才接受
            if (sum(cnt.values()) == sum(expected.values())
                    and all(cnt.get(el, 0) == exp for el, exp in expected.items())):
                return top_idx, boundary
        else:
            return top_idx, boundary

    raise ValueError("自动识别最上层失败！请用 --zthresh 手动指定 z 阈值。")


def write_poscar(lattice, elements, atom_counts, coords_frac, comment,
                 out_path, flags):
    """写入带 Selective Dynamics 的 POSCAR（分数坐标）。"""
    with open(out_path, 'w') as f:
        f.write(comment + '\n')
        f.write('1.0\n')
        for vec in lattice:
            f.write('  {:20.16f}  {:20.16f}  {:20.16f}\n'.format(*vec))
        if elements is not None:
            f.write('  ' + '  '.join(elements) + '\n')
        f.write('  ' + '  '.join(map(str, atom_counts)) + '\n')
        f.write('Selective dynamics\n')
        f.write('Direct\n')
        for coord, flag in zip(coords_frac, flags):
            f.write('  {:20.16f}  {:20.16f}  {:20.16f}  {}\n'.format(
                coord[0], coord[1], coord[2], flag))


def main():
    parser = argparse.ArgumentParser(description='生成 8×8=64 个滑移结构 POSCAR')
    parser.add_argument('-p', '--poscar', default='POSCAR', help='输入 POSCAR 路径')
    parser.add_argument('-d', '--ndiv', type=int, default=DEFAULT_NDIV,
                        help=f'晶格等分数（步长=1/NDIV），默认 {DEFAULT_NDIV}，'
                             f'即每方向 {DEFAULT_NDIV+1} 个点')
    parser.add_argument('-o', '--outdir', default='structures',
                        help='输出目录，默认 structures/')
    parser.add_argument('--zthresh', type=float, default=None,
                        help='手动指定顶层 z 分数坐标阈值（默认自动识别）')
    args = parser.parse_args()

    if not os.path.exists(args.poscar):
        sys.exit(f"错误：未找到输入文件 {args.poscar}")

    # ---------- 读取 ----------
    lattice, elements, atom_counts, coords, atom_types, comment = read_poscar(args.poscar)
    n = len(coords)
    print(f"读取 {args.poscar}: {n} 个原子")
    if elements is not None:
        print(f"  元素: {elements}  数量: {atom_counts}")

    # ---------- 1. 识别 Z 轴 ----------
    z_axis = detect_stacking_axis(lattice)
    inplane = [i for i in range(3) if i != z_axis]
    print(f"[1] Z 轴 = 第 {z_axis+1} 条晶格矢量，|c| = {np.linalg.norm(lattice[z_axis]):.4f} Å")
    print(f"    面内滑移方向 = 第 {inplane[0]+1} 条 (a1, |a1|={np.linalg.norm(lattice[inplane[0]]):.4f} Å) "
          f"和第 {inplane[1]+1} 条 (a2, |a2|={np.linalg.norm(lattice[inplane[1]]):.4f} Å)")

    # ---------- 2. 识别最上层 Mn2BrI ----------
    top_idx, boundary = detect_top_layer(coords, atom_types, z_axis,
                                         TOP_EXPECT, args.zthresh)
    top_cnt = Counter(atom_types[i] for i in top_idx)
    print(f"[2] 最上层（滑移层）: {len(top_idx)} 个原子, z > {boundary:.4f}")
    print(f"    成分: {dict(top_cnt)}")

    # ---------- 3. 生成滑移结构 ----------
    os.makedirs(args.outdir, exist_ok=True)
    ndiv = args.ndiv
    npts = ndiv + 1                 # 每方向点数: 0,1/ndiv,...,ndiv/ndiv
    count = 0
    for y in range(npts):                      # 沿 a2 方向 (对应文件夹名首位)
        dy_frac = y / ndiv
        for x in range(npts):                  # 沿 a1 方向 (对应文件夹名末位)
            dx_frac = x / ndiv

            shift = np.zeros(3)
            shift[inplane[0]] = dx_frac        # a1 方向滑移
            shift[inplane[1]] = dy_frac        # a2 方向滑移

            new_coords = coords.copy()
            new_coords[top_idx] += shift
            new_coords[top_idx] -= np.floor(new_coords[top_idx])   # 折叠回 [0,1)

            flags = ['F F T'] * n              # 所有原子: 固定 x,y, 只优化 z

            folder = f"{y}_{x}"
            out_path = os.path.join(args.outdir, f"POSCAR_{folder}")
            new_comment = (f"{comment} | slide: a1*{x}/{ndiv} + a2*{y}/{ndiv}"
                           f" | z-thresh={boundary:.4f}")
            write_poscar(lattice, elements, atom_counts, new_coords,
                         new_comment, out_path, flags)
            count += 1

    print(f"[3] 生成完成！共 {count} 个结构 → {args.outdir}/POSCAR_*")
    print(f"    例如: {os.path.join(args.outdir, 'POSCAR_0_0')}  和  "
          f"{os.path.join(args.outdir, f'POSCAR_{ndiv}_{ndiv}')}")
    print(f"注意: 步长为 1/{ndiv}；{ndiv}/{ndiv}=1 相当于回到原位"
          f"(与 0_0 等价)，按要求保留全部 {npts}×{npts} 个结构。")


if __name__ == '__main__':
    main()
