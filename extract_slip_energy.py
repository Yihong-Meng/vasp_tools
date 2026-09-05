#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_slip_energy.py — 提取滑移扫描(64 个 y_x 结构)的最终能量
=================================================

遍历当前目录下 y_x 文件夹（y,x = 0..7，共 8×8 = 64 个），
从每个文件夹的 OUTCAR 中提取最后一次离子步的自由能 TOTEN，
计算相对最低能量结构的能量差 E_rel (meV)。

输出：
  1) 控制台汇总 + 8×8 相对能量矩阵（便于直接看势能面形状）
  2) CSV 文件 slip_energy.csv，列：
     y, x, dy(a2步数/ndiv), dx(a1步数/ndiv), E_TOTEN(eV), E_without_entropy(eV),
     E_rel(meV), converged

说明：
  - 能量取最后一次出现的 "free  energy   TOTEN"（结构优化最后一步），
    与 summary_tru.py 一致；同时记录 "energy without entropy" 供参考。
  - converged 列：OUTCAR 中离子步收敛("reached required accuracy") 为 yes，
    否则为 no（可能没算完或 NSW 到上限）。
  - 缺少 OUTCAR 的文件夹标记为 missing，不参与最低能量比较。

用法：
  python3 extract_slip_energy.py                  # 在包含 y_x 文件夹的目录运行
  python3 extract_slip_energy.py --calcdir /path  # 指定计算目录
  python3 extract_slip_energy.py -d 7             # 与 gen_slip.py 一致
  python3 extract_slip_energy.py -o energy.csv    # 自定义输出文件名
"""

import os
import re
import sys
import csv
import argparse


def read_outcar(outcar_path):
    """
    读取 OUTCAR，返回 dict:
      toten          : 最后一次 TOTEN (float 或 None)
      without_entropy: 最后一次 energy without entropy (float 或 None)
      ionic_converged: 是否出现 "reached required accuracy"
      clean_end       : 是否出现 "General timing and accounting"
    """
    r = {'toten': None, 'without_entropy': None,
         'ionic_converged': False, 'clean_end': False}
    if not os.path.isfile(outcar_path):
        return r

    with open(outcar_path, 'r', errors='ignore') as f:
        text = f.read()

    toten = re.findall(r'free\s+energy\s+TOTEN\s*=\s*([-\d.]+)', text)
    if toten:
        r['toten'] = float(toten[-1])

    wo = re.findall(r'energy\s+without\s+entropy\s*=\s*([-\d.]+)', text)
    if wo:
        r['without_entropy'] = float(wo[-1])

    r['ionic_converged'] = 'reached required accuracy' in text
    r['clean_end'] = 'General timing and accounting' in text

    return r


def main():
    parser = argparse.ArgumentParser(description='提取滑移扫描64个结构的能量')
    parser.add_argument('--calcdir', default='.',
                        help='包含 y_x 文件夹的目录，默认当前目录')
    parser.add_argument('-d', '--ndiv', type=int, default=7,
                        help='晶格等分数（与 gen_slip.py 一致），默认 7 → 8×8=64')
    parser.add_argument('-o', '--out', default='slip_energy.csv',
                        help='输出 CSV 文件名，默认 slip_energy.csv')
    args = parser.parse_args()

    ndiv = args.ndiv
    npts = ndiv + 1
    base = args.calcdir

    if not os.path.isdir(base):
        sys.exit(f"错误：找不到目录 {base}")

    # rows 每个元素: [y, x, dy_frac, dx_frac, toten, without_entropy, e_rel, converged]
    rows = []
    missing = []

    for y in range(npts):
        for x in range(npts):
            folder = f"{y}_{x}"
            outcar = os.path.join(base, folder, 'OUTCAR')
            r = read_outcar(outcar)
            row = [y, x, y / ndiv, x / ndiv,
                   r['toten'], r['without_entropy'], None, r['ionic_converged']]
            if r['toten'] is None:
                missing.append(folder)
            rows.append(row)

    # ---- 计算相对能量 ----
    energies = [r[4] for r in rows if r[4] is not None]
    if not energies:
        print("没有找到任何 OUTCAR 能量数据。")
        sys.exit(1)

    e_min = min(energies)
    min_row = next(r for r in rows if r[4] == e_min)
    min_folder = f"{min_row[1]}_{min_row[0]}"          # x_y 命名
    for r in rows:
        if r[4] is not None:
            r[6] = (r[4] - e_min) * 1000.0

    # ---- 写 CSV ----
    with open(args.out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['y', 'x', 'dy_a2', 'dx_a1', 'E_TOTEN_eV',
                         'E_without_entropy_eV', 'E_rel_meV', 'converged'])
        for r in rows:
            if r[4] is None:
                writer.writerow([r[0], r[1], f"{r[2]:.6f}", f"{r[3]:.6f}",
                                 'NaN', 'NaN', 'NaN', 'missing'])
            else:
                wo = f"{r[5]:.8f}" if r[5] is not None else 'NaN'
                conv = 'yes' if r[7] else 'no'
                writer.writerow([r[0], r[1], f"{r[2]:.6f}", f"{r[3]:.6f}",
                                 f"{r[4]:.8f}", wo, f"{r[6]:.3f}", conv])

    # ---- 控制台汇总 ----
    print(f"共 {len(rows)} 个结构，其中 {len(energies)} 个有能量，"
          f"{len(missing)} 个缺失 OUTCAR")
    print(f"最低能量: {min_folder}  E_min = {e_min:.8f} eV\n")

    print("相对能量矩阵 E_rel (meV)，行=y(a2方向)，列=x(a1方向)：")
    print("     " + "".join(f"{x:9d}" for x in range(npts)))
    print("     " + "-" * (9 * npts))
    for y in range(npts):
        line = f"y={y}  "
        for x in range(npts):
            r = rows[y * npts + x]
            line += f"{r[6]:9.3f}" if r[6] is not None else f"{'NaN':>9s}"
        print(line)
    print()

    # ---- 缺失/未收敛提示 ----
    if missing:
        print(f"警告：以下 {len(missing)} 个文件夹没有能量数据: {missing}")
    notconv = [f"{r[1]}_{r[0]}" for r in rows
               if r[4] is not None and not r[7]]
    if notconv:
        print(f"警告：以下 {len(notconv)} 个结构离子步未收敛: {notconv}")

    print(f"\n结果已保存至 {args.out}")


if __name__ == '__main__':
    main()
