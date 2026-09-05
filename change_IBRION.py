#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
change_IBRION.py — 对未收敛的滑移结构换算法重启
================================================

对所有 y_x 文件夹，检查 OUTCAR 中离子步是否收敛。
未收敛的执行：
  1. INCAR 中 IBRION 2 → 1（共轭梯度 → 拟牛顿法）
  2. INCAR 中 POTIM → 0.2（降低步长，抑制震荡）
  3. CONTCAR → POSCAR（从上次断点继续）
  4. 删除 OUTCAR, vasprun.xml, OSZICAR, .jobid（让 submit_slip.sh 重新拾取）
  5. 保留 WAVECAR（加速重启电子步）

已收敛的不动。没有 OUTCAR 的不动。IBRION 本来就不是 2 的会警告。

用法：
  python3 change_IBRION.py                     # 当前目录
  python3 change_IBRION.py /path/to            # 指定目录
  python3 change_IBRION.py -d 7                # ndiv=7 → 8×8=64
  python3 change_IBRION.py --dry-run           # 预览，不执行
  python3 change_IBRION.py --target 1          # IBRION 目标值（默认1）
  python3 change_IBRION.py --potim 0.1         # POTIM 目标值（默认0.2）
"""

import os
import re
import sys
import shutil
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
        return None  # 没有 OUTCAR
    try:
        with open(outcar_path, 'r', errors='ignore') as f:
            content = f.read()
        return 'reached required accuracy' in content
    except:
        return None


def change_incar_param(incar_path, param_name, old_val, new_val):
    """
    修改 INCAR 中某个整数参数的值。
    支持格式：PARAM = 2  /  PARAM= 2  /  PARAM   2
    返回 (是否修改成功, 修改前的值)。
    """
    if not os.path.isfile(incar_path):
        return False, None

    with open(incar_path, 'r') as f:
        content = f.read()

    pattern = re.compile(rf'({param_name}\s*[=]?\s*)(\d+)(.*)', re.IGNORECASE)
    match = pattern.search(content)
    if not match:
        return False, None

    old = int(match.group(2))
    if old != old_val:
        return False, old

    new_content = pattern.sub(rf'\g<1>{new_val}\3', content, count=1)

    with open(incar_path, 'w') as f:
        f.write(new_content)

    return True, old


def change_incar_float(incar_path, param_name, new_val):
    """
    修改或添加 INCAR 中某个浮点参数。
    如果存在则替换值，不存在则在文件末尾追加。
    返回 (是否为修改已有行, 修改/添加的值)。
    """
    if not os.path.isfile(incar_path):
        return False, new_val

    with open(incar_path, 'r') as f:
        content = f.read()

    # 先检查是否已存在
    pattern = re.compile(rf'({param_name}\s*[=]?\s*)([\d.]+)(.*)', re.IGNORECASE)
    match = pattern.search(content)
    if match:
        old_val = float(match.group(2))
        new_content = pattern.sub(rf'\g<1>{new_val}\3', content, count=1)
        with open(incar_path, 'w') as f:
            f.write(new_content)
        return True, old_val
    else:
        # 不存在，追加到末尾
        with open(incar_path, 'a') as f:
            f.write(f"\n{param_name} = {new_val}\n")
        return False, None


def main():
    parser = argparse.ArgumentParser(
        description='对未收敛的文件夹：改 IBRION/POTIM + CONTCAR→POSCAR + 清理重启')
    parser.add_argument('calcdir', nargs='?', default='.',
                        help='包含 y_x 文件夹的目录（默认当前目录）')
    parser.add_argument('-d', '--ndiv', type=int, default=6,
                        help='晶格等分数（默认6 → 7×7=49）')
    parser.add_argument('--target', type=int, default=1,
                        help='IBRION 目标值（默认1：拟牛顿法）')
    parser.add_argument('--potim', type=float, default=0.2,
                        help='POTIM 目标值（默认0.2，降低步长抑制震荡）')
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
    if args.dry_run:
        print("[DRY RUN] 只预览，不执行")
    print()

    converged = 0
    restarted = 0
    skipped_no_outcar = 0
    skipped_no_contcar = 0
    skipped_ibrion_other = 0

    for folder in folders:
        d = os.path.join(base, folder)
        outcar = os.path.join(d, 'OUTCAR')
        contcar = os.path.join(d, 'CONTCAR')
        poscar = os.path.join(d, 'POSCAR')
        incar = os.path.join(d, 'INCAR')

        # 无 OUTCAR → 跳过
        if not os.path.isfile(outcar):
            skipped_no_outcar += 1
            continue

        # 已收敛 → 跳过
        result = check_converged(outcar)
        if result is True:
            converged += 1
            continue

        # 未收敛，需要重启
        # 检查 CONTCAR
        if not os.path.isfile(contcar):
            print(f"  [警告] {folder}: 未收敛但没有CONTCAR，跳过")
            skipped_no_contcar += 1
            continue

        # 修改 IBRION
        ok_ibrion, old_ibrion = change_incar_param(incar, 'IBRION', 2, args.target)

        if not ok_ibrion:
            if old_ibrion is not None:
                print(f"  [跳过] {folder}: IBRION={old_ibrion}（不是2），跳过修改")
                skipped_ibrion_other += 1
                continue
            else:
                print(f"  [警告] {folder}: INCAR中未找到IBRION行")
                skipped_ibrion_other += 1
                continue

        # 修改 POTIM
        was_existing, old_potim = change_incar_float(incar, 'POTIM', args.potim)

        # 执行重启操作
        restarted += 1
        potim_info = f"POTIM={old_potim}→{args.potim}" if was_existing else f"添加POTIM={args.potim}"

        if args.dry_run:
            print(f"  [将重启] {folder}: IBRION 2→{args.target}, {potim_info}, CONTCAR→POSCAR")
        else:
            # CONTCAR → POSCAR
            shutil.copy2(contcar, poscar)
            # 清理重启标记
            for fname in ['OUTCAR', 'vasprun.xml', 'OSZICAR', '.jobid']:
                fpath = os.path.join(d, fname)
                if os.path.isfile(fpath):
                    os.remove(fpath)
            print(f"  [已重启] {folder}: IBRION 2→{args.target}, {potim_info}, CONTCAR→POSCAR")

    # ---- 汇总 ----
    print()
    print("===== 汇总 =====")
    print(f"  已收敛(跳过):       {converged} 个")
    print(f"  需要重启(IBRION改): {restarted} 个")
    if skipped_no_outcar:
        print(f"  无OUTCAR(跳过):     {skipped_no_outcar} 个")
    if skipped_no_contcar:
        print(f"  缺CONTCAR(跳过):    {skipped_no_contcar} 个")
    if skipped_ibrion_other:
        print(f"  IBRION非2(跳过):    {skipped_ibrion_other} 个")

    if not args.dry_run and restarted > 0:
        print()
        print("下一步: bash submit_slip.sh   # 重新提交")


if __name__ == '__main__':
    main()
