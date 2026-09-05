#!/usr/bin/env python3
"""
check_results.py — 批量检查滑移扫描文件夹的计算收敛情况
======================================================

自动扫描当前目录下所有 y_x 形式的文件夹（如 0_0, 3_5, 6_6），
检查每个文件夹中 VASP 计算（OUTCAR）的收敛状态：
  - 是否正常结束（General timing and accounting）
  - 离子步是否收敛（reached required accuracy）
  - VASP 运行错误
  - SLURM 提交错误（.err 文件）
  - 最终能量、磁矩、离子步数

输出汇总表格 + 排序后的能量列表 + 问题文件夹清单。

用法：
  python3 check_results.py                    # 当前目录
  python3 check_results.py /path/to/calcs     # 指定目录
  python3 check_results.py -d 7               # 指定网格大小（默认6 → 7×7）
  python3 check_results.py --no-energy        # 不输出能量列表
"""

import os
import sys
import re
import argparse
from collections import OrderedDict


def check_one(calc_dir, folder):
    """检查单个计算文件夹的 OUTCAR 状态。"""
    r = {"folder": folder, "has_outcar": False, "clean_end": False,
         "ionic_converged": False, "n_ionic_steps": 0,
         "final_energy": None, "wo_energy": None, "mag_moment": None,
         "error_msg": None, "slurm_err": None, "running": False}

    outcar = os.path.join(calc_dir, "OUTCAR")

    if not os.path.isfile(outcar):
        # 检查是否任务还在运行（slurm-*.out 存在且 OUTCAR 未完成）
        slurm_outs = [f for f in os.listdir(calc_dir) if f.startswith("slurm-") and f.endswith(".out")]
        if slurm_outs:
            r["running"] = True
            r["error_msg"] = "计算仍在运行"
        else:
            r["error_msg"] = "No OUTCAR"
        return r
    r["has_outcar"] = True

    try:
        with open(outcar, errors='ignore') as f:
            content = f.read()
    except:
        r["error_msg"] = "Cannot read OUTCAR"
        return r

    r["clean_end"] = "General timing and accounting" in content
    r["ionic_converged"] = "reached required accuracy" in content

    # VASP 错误
    for pat in [r"ZBRENT\s+AND\s+FUNC\s+ABORTING", r"BRENT\s+ABORTING",
                r"TOO\s+MANY\s+STEPS", r"internal\s+error",
                r"gradient\s+not\s+finite", r"I\s+REFUSE\s+TO\s+CONTINUE"]:
        if re.search(pat, content, re.IGNORECASE):
            idx = content.lower().find(pat.split()[0].lower()[:10])
            r["error_msg"] = f"VASP error: {content[max(0,idx-30):idx+80].strip()}"
            return r

    # 离子步数 & 能量
    energies = re.findall(r"energy\s+without\s+entropy\s*=\s*([-\d.]+)", content)
    if energies:
        r["n_ionic_steps"] = len(energies)
        r["wo_energy"] = float(energies[-1])

    toten = re.findall(r"free\s+energy\s+TOTEN\s*=\s*([-\d.]+)", content)
    if toten:
        r["final_energy"] = float(toten[-1])

    # 磁矩（取最后一个）
    mags = re.findall(r"total magnetization\s*\(x\)\s*:\s*([-\d.]+)", content)
    if mags:
        r["mag_moment"] = float(mags[-1])

    # SLURM .err
    for fname in os.listdir(calc_dir):
        if fname.endswith(".err") and not fname.startswith("slurm-"):
            try:
                with open(os.path.join(calc_dir, fname)) as fh:
                    err = fh.read().strip()
                if err:
                    r["slurm_err"] = err[:300]
            except:
                pass
            break

    return r


def print_result(r):
    """格式化单个结果行。"""
    if not r["has_outcar"]:
        tag = "[运行中]" if r["running"] else "[无OUTCAR]"
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
        issues.append(r["error_msg"][:60])
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
    return folders


def main():
    parser = argparse.ArgumentParser(description='批量检查滑移扫描计算收敛情况')
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
    has_outcar = sum(1 for r in results if r["has_outcar"])
    converged = sum(1 for r in results if r["ionic_converged"])
    clean = sum(1 for r in results if r["clean_end"])
    running = sum(1 for r in results if r["running"])
    missing = sum(1 for r in results if not r["has_outcar"] and not r["running"])

    print("=" * 90)
    print(f"汇总: {total} 个文件夹")
    print(f"  有 OUTCAR:    {has_outcar} 个")
    print(f"  正常结束:     {clean} 个")
    print(f"  离子步收敛:   {converged} 个")
    if running:
        print(f"  仍在运行:     {running} 个")
    if missing:
        print(f"  缺少 OUTCAR:  {missing} 个")

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
            if not r["has_outcar"] and not r["running"]:
                reasons.append("无OUTCAR")
            if not r["clean_end"]:
                reasons.append("异常结束")
            if not r["ionic_converged"] and r["has_outcar"]:
                reasons.append(f"未收敛({r['n_ionic_steps']}步)")
            if r["error_msg"] and r["error_msg"] != "计算仍在运行":
                reasons.append(r["error_msg"][:40])
            if r["slurm_err"]:
                reasons.append(r["slurm_err"][:40])
            print(f"  {r['folder']}: {', '.join(reasons)}")


if __name__ == "__main__":
    main()
