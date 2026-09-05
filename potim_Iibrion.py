#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
potim_Iibrion.py — 批量修改未收敛/未计算文件夹的 INCAR
=======================================================

逻辑：
  1. 扫描所有 y_x 文件夹
  2. 检查 OUTCAR 中是否已收敛
  3. 已收敛 → 跳过
  4. 未收敛或从未提交 → 修改 INCAR：IBRION=1, POTIM=0.2
     - IBRION：若已有则改值，若没有则追加
     - POTIM：若已有则改值，若没有则追加

用法：
  python3 potim_Iibrion.py                     # 当前目录
  python3 potim_Iibrion.py /path/to            # 指定目录
  python3 potim_Iibrion.py -d 7                # ndiv=7 → 8×8=64
  python3 potim_Iibrion.py --dry-run           # 预览，不执行
  python3 potim_Iibrion.py --ibrion 1          # IBRION 目标值（默认1）
  python3 potim_Iibrion.py --potim 0.2         # POTIM 目标值（默认0.2）
"""

import os
import re
import sys
import argparse


def find_yx_folders(base, ndiv):
    """自动发现 y_x 文件夹。"""
    folders = []
    for y in range(ndiv + 1):
        for x in range(ndiv + 1):
            d = f"{y}_{x}"
            if os.path.isdir(os.path.join(base, d)):
                folders.append(d)
    return folders


def check_converged(outcar_path):
    """检查 OUTCAR 中离子步是否收敛。"""
    if not os.path.isfile(outcar_path):
        return None
    try:
        with open(outcar_path, 'r', errors='ignore') as f:
            return 'reached required accuracy' in f.read()
    except:
        return None


def set_incar_int(incar_path, param, value):
    """
    设置 INCAR 整数参数。
    返回 'modified'(修改已有) / 'added'(追加) / 'not_found'(INCAR不存在)。
    """
    if not os.path.isfile(incar_path):
        return 'not_found'

    with open(incar_path, 'r') as f:
        content = f.read()

    pattern = re.compile(rf'^{param}\s*[=]?\s*\d+.*$', re.MULTILINE | re.IGNORECASE)
    if pattern.search(content):
        new_content = re.sub(rf'^({param}\s*[=]?\s*)\d+', rf'\g<1>{value}',
                             content, count=1, flags=re.MULTILINE | re.IGNORECASE)
        with open(incar_path, 'w') as f:
            f.write(new_content)
        return 'modified'
    else:
        with open(incar_path, 'a') as f:
            f.write(f"\n{param} = {value}\n")
        return 'added'


def set_incar_float(incar_path, param, value):
    """
    设置 INCAR 浮点参数。
    返回 'modified' / 'added' / 'not_found'。
    """
    if not os.path.isfile(incar_path):
        return 'not_found'

    with open(incar_path, 'r') as f:
        content = f.read()

    pattern = re.compile(rf'^{param}\s*[=]?\s*[\d.eE+-]+.*$', re.MULTILINE | re.IGNORECASE)
    if pattern.search(content):
        new_content = re.sub(rf'^({param}\s*[=]?\s*)[\d.eE+-]+',
                             rf'\g<1>{value}', content, count=1,
                             flags=re.MULTILINE | re.IGNORECASE)
        with open(incar_path, 'w') as f:
            f.write(new_content)
        return 'modified'
    else:
        with open(incar_path, 'a') as f:
            f.write(f"\n{param} = {value}\n")
        return 'added'


def main():
    parser = argparse.ArgumentParser(description='批量修改未收敛/未计算文件夹的 INCAR')
    parser.add_argument('calcdir', nargs='?', default='.',
                        help='包含 y_x 文件夹的目录（默认当前目录）')
    parser.add_argument('-d', '--ndiv', type=int, default=6,
                        help='晶格等分数（默认6 → 7×7=49）')
    parser.add_argument('--ibrion', type=int, default=1,
                        help='IBRION 目标值（默认1）')
    parser.add_argument('--potim', type=float, default=0.2,
                        help='POTIM 目标值（默认0.2）')
    parser.add_argument('--dry-run', action='store_true',
                        help='预览模式，不执行任何操作')
    args = parser.parse_args()

    base = os.path.abspath(args.calcdir)
    if not os.path.isdir(base):
        sys.exit(f"错误：目录不存在 {base}")

    folders = find_yx_folders(base, args.ndiv)
    if not folders:
        sys.exit(f"未找到 y_x 文件夹（ndiv={args.ndiv}）")

    print(f"目录: {base}")
    print(f"网格: {args.ndiv+1}×{args.ndiv+1}，共 {len(folders)} 个文件夹")
    print(f"目标: IBRION={args.ibrion}, POTIM={args.potim}")
    if args.dry_run:
        print("[DRY RUN] 只预览，不执行")
    print()

    converged = 0
    modified = 0
    no_outcar = 0  # 从未提交过

    for folder in folders:
        d = os.path.join(base, folder)
        outcar = os.path.join(d, 'OUTCAR')
        incar = os.path.join(d, 'INCAR')

        # 检查收敛状态
        result = check_converged(outcar)
        if result is True:
            converged += 1
            continue

        # 未收敛或从未提交 → 改 INCAR
        has_outcar = os.path.isfile(outcar)
        tag = "[未提交]" if not has_outcar else "[未收敛]"

        if args.dry_run:
            print(f"  {tag} {folder}  →  改 INCAR: IBRION={args.ibrion}, POTIM={args.potim}")
            modified += 1
            continue

        if not has_outcar:
            no_outcar += 1

        r1 = set_incar_int(incar, 'IBRION', args.ibrion)
        r2 = set_incar_float(incar, 'POTIM', args.potim)
        modified += 1

        parts = []
        if r1 == 'modified': parts.append('IBRION已改')
        elif r1 == 'added': parts.append('IBRION已加')
        if r2 == 'modified': parts.append('POTIM已改')
        elif r2 == 'added': parts.append('POTIM已加')
        info = ', '.join(parts) if parts else '无INCAR'

        print(f"  {tag} {folder}  →  {info}")

    # ---- 汇总 ----
    print()
    print("===== 汇总 =====")
    print(f"  已收敛(跳过):       {converged} 个")
    print(f"  已修改 INCAR:       {modified} 个")
    if no_outcar:
        print(f"  其中从未提交过:     {no_outcar} 个")

    if not args.dry_run and modified > 0:
        print()
        print("下一步: bash submit_slip.sh   # 提交/重新提交")


if __name__ == '__main__':
    main()
