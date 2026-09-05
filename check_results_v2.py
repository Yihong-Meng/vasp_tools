#!/usr/bin/env python3
"""
check_results_v2.py — 批量检查滑移扫描文件夹的计算收敛情况（改进版）
=================================================================

改进点：
1. 通过检查log文件最后一行判断收敛
2. 收敛判断依据：最后一行包含 "reached required accuracy - stopping structural energy minimisation"
3. ZBRENT错误输出特定提示
4. 其他情况输出"不收敛"

自动扫描当前目录下所有 y_x 形式的文件夹（如 0_0, 3_5, 6_6），
检查每个文件夹中 VASP 计算（log文件）的收敛状态：
  - 最后一行是否包含收敛标志
  - ZBRENT错误检测
  - 从OUTCAR读取能量

输出汇总表格 + 排序后的能量列表 + 问题文件夹清单。

用法：
  python3 check_results_v2.py                    # 当前目录
  python3 check_results_v2.py /path/to/calcs     # 指定目录
  python3 check_results_v2.py -d 7               # 指定网格大小（默认6 → 7×7）
  python3 check_results_v2.py --no-energy        # 不输出能量列表
"""

import os
import sys
import re
import argparse
from collections import OrderedDict


def check_one(calc_dir, folder):
    """检查单个计算文件夹的 log 文件状态。"""
    r = {"folder": folder, "has_log": False, "clean_end": False,
         "ionic_converged": False, "n_ionic_steps": 0,
         "final_energy": None, "wo_energy": None, "mag_moment": None,
         "error_msg": None, "slurm_err": None, "running": False}

    # 查找log文件（可能是log、slurm-*.out或*.log）
    log_file = None
    for fname in os.listdir(calc_dir):
        if fname == "log":  # 直接检查文件名是否为"log"
            log_file = os.path.join(calc_dir, fname)
            break
        elif fname.startswith("slurm-") and fname.endswith(".out"):
            log_file = os.path.join(calc_dir, fname)
            break
        elif fname.endswith(".log"):
            log_file = os.path.join(calc_dir, fname)
            break
    
    if not log_file:
        r["error_msg"] = "No log file"
        return r
    r["has_log"] = True

    try:
        with open(log_file, errors='ignore') as f:
            content = f.read()
            # 读取最后一行
            lines = content.splitlines()
            last_line = lines[-1].strip() if lines else ""
    except:
        r["error_msg"] = "Cannot read log file"
        return r

    # 根据最后一行判断收敛状态
    if "reached required accuracy - stopping structural energy minimisation" in last_line:
        r["ionic_converged"] = True
        r["clean_end"] = True
    elif "ZBRENT: fatal error in bracketing" in content:
        r["error_msg"] = "please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue"
        return r
    else:
        # 其他情况，检查是否正常结束
        r["clean_end"] = "General timing and accounting" in content
        if not r["clean_end"]:
            r["error_msg"] = "不收敛"
        return r

    # 如果收敛了，从OUTCAR读取能量
    outcar = os.path.join(calc_dir, "OUTCAR")
    if os.path.isfile(outcar):
        try:
            with open(outcar, errors='ignore') as f:
                outcar_content = f.read()
            
            # 离子步数 & 能量
            energies = re.findall(r"energy\s+without\s+entropy\s*=\s*([-\d.]+)", outcar_content)
            if energies:
                r["n_ionic_steps"] = len(energies)
                r["wo_energy"] = float(energies[-1])

            toten = re.findall(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)", outcar_content)
            if toten:
                r["final_energy"] = float(toten[-1])

            # 如果没有找到能量，尝试其他模式
            if not r["final_energy"] and not r["wo_energy"]:
                # 尝试从F=行提取能量
                f_energy = re.findall(r"F\s*=\s*([-\d.]+)", outcar_content)
                if f_energy:
                    r["wo_energy"] = float(f_energy[-1])
                    r["n_ionic_steps"] = len(f_energy)
                
                # 尝试从TOTEN行提取
                toten2 = re.findall(r"TOTEN\s*=\s*([-\d.]+)", outcar_content)
                if toten2:
                    r["final_energy"] = float(toten2[-1])

            # 磁矩（取最后一个）
            mags = re.findall(r"total magnetization\s*\(x\)\s*:\s*([-\d.]+)", outcar_content)
            if mags:
                r["mag_moment"] = float(mags[-1])
                
        except:
            pass

    # SLURM .err
    for fname in os.listdir(calc_dir):
        if fname.endswith(".err") and not fname.startswith("slurm-"):
            try:
                with open(os.path.join(calc_dir, fname)) as fh:
                    err = fh.read().strip()
                if err:
                    # 过滤掉常见的非关键错误
                    if "whoami: cannot find name for user ID" in err:
                        # 这个错误通常不影响计算，忽略它
                        pass
                    else:
                        r["slurm_err"] = err[:300]
            except:
                pass
            break

    return r


def print_result(r):
    """格式化单个结果行。"""
    if not r["has_log"]:
        tag = "[运行中]" if r["running"] else "[无log文件]"
        return f"{r['folder']:6s}  {tag}"

    clean = "OK" if r["clean_end"] else "异常"
    ionic = "OK" if r["ionic_converged"] else "FAIL"
    energy = f"{r['final_energy']:.6f}" if r["final_energy"] else (
             f"{r['wo_energy']:.6f}" if r['wo_energy'] else "--")
    nsteps = f"{r['n_ionic_steps']}" if r["n_ionic_steps"] else "--"

    issues = []
    if not r["clean_end"]:
        issues.append("异常结束")
    if not r["ionic_converged"]:
        issues.append("未收敛")
    if r["error_msg"] and r["error_msg"] != "计算仍在运行":
        issues.append(r["error_msg"][:80])
    if r["slurm_err"]:
        issues.append(f"SLURM: {r['slurm_err'][:60]}")

    status = "OK" if not issues else "; ".join(issues)

    return (f"{r['folder']:6s}  结束:{clean:4s}  离子步:{ionic:4s}  "
            f"E={energy:>15s}  {nsteps:>3s}步  {status}")


def discover_folders(base, ndiv=6):
    """自动发现 y_x 文件夹，按 y_x 顺序排列。"""
    npts = ndiv + 1
    folders = []
    for y in range(npts):
        for x in range(npts):
            d = f"{y}_{x}"
            if os.path.isdir(os.path.join(base, d)):
                folders.append(d)
    
    # 如果没有找到y_x格式的文件夹，尝试查找所有子文件夹
    if not folders:
        for item in os.listdir(base):
            item_path = os.path.join(base, item)
            if os.path.isdir(item_path):
                folders.append(item)
    
    return folders


def main():
    parser = argparse.ArgumentParser(description='批量检查滑移扫描计算收敛情况（改进版）')
    parser.add_argument('calcdir', nargs='?', default='.',
                        help='包含 y_x 文件夹的目录（默认当前目录）')
    parser.add_argument('-d', '--ndiv', type=int, default=6,
                        help='晶格等分数（默认6 → 7×7=49 个文件夹）')
    parser.add_argument('--no-energy', action='store_true',
                        help='不输出能量排序列表')
    args = parser.parse_args()

    base = os.path.abspath(args.calcdir)
    if not os.path.isdir(base):
        print(f"错误：目录不存在 {base}")
        sys.exit(1)

    folders = discover_folders(base, args.ndiv)
    if not folders:
        print(f"在 {base} 下未找到 y_x 文件夹（ndiv={args.ndiv}）")
        sys.exit(1)

    npts = args.ndiv + 1
    print(f"扫描目录: {base}")
    print(f"网格: {npts}×{npts} (ndiv={args.ndiv})，找到 {len(folders)} 个文件夹\n")

    # ---- 逐个检查 ----
    print(f"{'文件夹':6s}  {'结束':4s}  {'离子步':4s}  {'Energy(eV)':>15s}  {'步数':>3s}  状态")
    print("=" * 90)

    results = []
    for f in folders:
        d = os.path.join(base, f)
        r = check_one(d, f)
        results.append(r)
        print(print_result(r))

    # ---- 汇总统计 ----
    total = len(results)
    has_log = sum(1 for r in results if r["has_log"])
    converged = sum(1 for r in results if r["ionic_converged"])
    clean = sum(1 for r in results if r["clean_end"])
    running = sum(1 for r in results if r["running"])
    missing = sum(1 for r in results if not r["has_log"] and not r["running"])

    print("=" * 90)
    print(f"汇总: {total} 个文件夹")
    print(f"  有 log 文件:  {has_log} 个")
    print(f"  正常结束:     {clean} 个")
    print(f"  离子步收敛:   {converged} 个")
    if running:
        print(f"  仍在运行:     {running} 个")
    if missing:
        print(f"  缺少 log 文件: {missing} 个")

    # ---- 能量排序 ----
    if not args.no_energy:
        valid = [r for r in results if r["final_energy"] is not None]
        if valid:
            valid.sort(key=lambda r: r["final_energy"])
            e_min = valid[0]["final_energy"]
            print(f"\n能量排序（E_min = {e_min:.6f} eV）：")
            print(f"{'排名':>4s}  {'文件夹':6s}  {'E(eV)':>15s}  {'dE(meV)':>10s}")
            print("-" * 45)
            for i, r in enumerate(valid):
                dE = (r["final_energy"] - e_min) * 1000
                print(f"{i+1:>4d}  {r['folder']:6s}  {r['final_energy']:>15.6f}  {dE:>10.3f}")

    # ---- 问题文件夹清单 ----
    problems = [r for r in results if not r["ionic_converged"] or r["error_msg"]
                or r["slurm_err"] or r["running"]]
    if problems:
        print(f"\n有问题的文件夹 ({len(problems)} 个)：")
        for r in problems:
            reasons = []
            if r["running"]:
                reasons.append("运行中")
            if not r["has_log"] and not r["running"]:
                reasons.append("无log文件")
            if not r["clean_end"]:
                reasons.append("异常结束")
            if not r["ionic_converged"] and r["has_log"]:
                reasons.append(f"未收敛({r['n_ionic_steps']}步)")
            if r["error_msg"] and r["error_msg"] != "计算仍在运行":
                reasons.append(r["error_msg"][:60])
            if r["slurm_err"]:
                reasons.append(r["slurm_err"][:40])
            print(f"  {r['folder']}: {', '.join(reasons)}")


if __name__ == "__main__":
    main()