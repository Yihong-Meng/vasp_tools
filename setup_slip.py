#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_slip.py — 为 64 个滑移结构创建计算文件夹
=================================================

对 gen_slip.py 生成的 structures/POSCAR_{y}_{x}（y,x = 0..7），
在当前目录创建 y_x 计算文件夹，并放入 VASP 输入文件：

  structures/POSCAR_{y}_{x}  →  y_x/POSCAR         （复制）
  INCAR                      →  y_x/INCAR          （复制）
  vasp_stru.sh               →  y_x/vasp_stru.sh   （复制，保持可执行）
  POTCAR                     →  y_x/POTCAR         （软链接）
  KPOINTS                    →  y_x/KPOINTS        （软链接）

用法：
  python3 setup_slip.py                 # 默认 8×8 = 64 个
  python3 setup_slip.py -d 7            # 与 gen_slip.py 的 -d 保持一致
  python3 setup_slip.py --src structures --calcdir .
"""

import os
import sys
import argparse
import shutil
import stat


def main():
    parser = argparse.ArgumentParser(description='创建滑移计算文件夹')
    parser.add_argument('-d', '--ndiv', type=int, default=7,
                        help='晶格等分数（与 gen_slip.py 一致），默认 7 → 8×8=64')
    parser.add_argument('--src', default='structures',
                        help='POSCAR 源目录，默认 structures/')
    parser.add_argument('--calcdir', default='.',
                        help='计算文件夹创建位置，默认当前目录')
    args = parser.parse_args()

    ndiv = args.ndiv
    npts = ndiv + 1

    # ---- 检查模板文件 ----
    required = ['INCAR', 'POTCAR', 'KPOINTS', 'vasp_stru.sh']
    missing = [f for f in required if not os.path.isfile(f)]
    if missing:
        sys.exit(f"错误：当前目录缺少文件: {', '.join(missing)}"
                 f"\n请先在当前目录准备好 INCAR / POTCAR / KPOINTS / vasp_stru.sh")

    src_dir = args.src
    if not os.path.isdir(src_dir):
        sys.exit(f"错误：找不到源目录 {src_dir}，请先运行 gen_slip.py")

    os.makedirs(args.calcdir, exist_ok=True)

    total = npts * npts
    created = 0
    skipped = 0

    print(f"源目录: {src_dir}/")
    print(f"目标目录: {args.calcdir}/  (共 {total} 个 y_x 文件夹)")

    for y in range(npts):
        for x in range(npts):
            folder = f"{y}_{x}"
            src_poscar = os.path.join(src_dir, f"POSCAR_{folder}")
            if not os.path.isfile(src_poscar):
                print(f"  [跳过] {folder}: 缺少 {src_poscar}")
                skipped += 1
                continue

            calc_dir = os.path.join(args.calcdir, folder)
            os.makedirs(calc_dir, exist_ok=True)

            # 复制 POSCAR
            shutil.copy2(src_poscar, os.path.join(calc_dir, 'POSCAR'))

            # 复制 INCAR
            shutil.copy2('INCAR', os.path.join(calc_dir, 'INCAR'))

            # 复制 vasp_stru.sh 并确保可执行
            dst_sh = os.path.join(calc_dir, 'vasp_stru.sh')
            shutil.copy2('vasp_stru.sh', dst_sh)
            st = os.stat(dst_sh)
            os.chmod(dst_sh, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

            # 软链接 POTCAR、KPOINTS（若尚不存在）
            for fname in ('POTCAR', 'KPOINTS'):
                dst = os.path.join(calc_dir, fname)
                if not os.path.lexists(dst):
                    rel = os.path.relpath(os.path.abspath(fname), calc_dir)
                    os.symlink(rel, dst)

            created += 1

    print(f"\n完成！创建 {created} 个文件夹，跳过 {skipped} 个。")
    print("下一步: nohup bash submit_slip.sh > submit.log 2>&1 &")


if __name__ == '__main__':
    main()
