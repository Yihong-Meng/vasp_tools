#!/bin/bash
# ============================================================
# submit_slip.sh — 批量提交滑移结构计算 (SLURM, 最多同时运行 4 个)
#
# 用途: 自动提交 y_x 文件夹（y,x = 0..7，共 64 个）中的计算任务。
#       同一时刻最多只有 MAX_JOBS 个任务在运行，算完一个自动补交下一个。
#
# 用法:
#   bash submit_slip.sh                  # 前台运行
#   nohup bash submit_slip.sh > submit.log 2>&1 &   # 后台运行(推荐)
#   tail -f submit.log
#
# 说明:
#   - 收集当前目录下所有 y_x 形式的文件夹（如 0_0, 3_5, 7_7）。
#   - 断点续跑：重新运行时自动跳过
#       (a) 已有 vasprun.xml 的文件夹（已完成）；
#       (b) 通过文件夹内 .jobid 查到仍在 SLURM 队列中运行的任务（避免重复提交）。
#     失败/被取消的任务（无 vasprun.xml 且不在队列）会被重新提交。
#   - 任务提交脚本使用各文件夹内的 vasp_stru.sh（sbatch vasp_stru.sh）。
# ============================================================

MAX_JOBS=4           # 同时运行的最大任务数
CHECK_INTERVAL=30    # 状态检查间隔(秒)

echo "===== 批量提交滑移结构计算 (最多 $MAX_JOBS 个并行) ====="
echo "检查间隔: ${CHECK_INTERVAL}s"

# ---- 1. 收集所有 y_x 文件夹 ----
mapfile -t all_dirs < <(ls -d [0-9]_[0-9] 2>/dev/null | sort -t_ -k1,1n -k2,2n)

# ---- 2. 过滤已完成 / 已在队列中的文件夹 ----
pending=()
for d in "${all_dirs[@]}"; do
    if [ -f "$d/vasprun.xml" ]; then
        echo "  [跳过] $d 已有 vasprun.xml（已完成）"
    elif [ -f "$d/.jobid" ]; then
        # 上次提交过：查一下任务是否还在队列里
        oldjid=$(cat "$d/.jobid")
        if squeue -j "$oldjid" --noheader 2>/dev/null | grep -q .; then
            echo "  [跳过] $d 已有任务在运行 (JobID=$oldjid)"
        else
            # 已离开队列但没有 vasprun.xml → 之前失败/被取消，重新提交
            rm -f "$d/.jobid"
            pending+=("$d")
        fi
    else
        pending+=("$d")
    fi
done
all_dirs=("${pending[@]}")

total=${#all_dirs[@]}
if [ $total -eq 0 ]; then
    echo "[提示] 没有待提交的 y_x 任务。"
    echo "       请先运行: python3 gen_slip.py && python3 setup_slip.py"
    exit 0
fi
echo "共 $total 个任务待计算"
echo ""

# ---- 3. 状态跟踪 ----
declare -A running    # jobid -> dir
next_idx=0
success_count=0
fail_count=0

submit_one() {
    # 没有更多任务可提交
    if [ $next_idx -ge $total ]; then return 1; fi

    local dir="${all_dirs[$next_idx]}"
    next_idx=$((next_idx + 1))

    # 检查 vasp_stru.sh 是否存在
    if [ ! -f "$dir/vasp_stru.sh" ]; then
        printf "  [跳过] %-10s (缺少 vasp_stru.sh)\n" "$dir"
        fail_count=$((fail_count + 1))
        return 0
    fi

    # sbatch 提交
    local out jid
    out=$(cd "$dir" && sbatch vasp_stru.sh 2>&1)
    jid=$(echo "$out" | grep -oE '[0-9]+$')

    if [ -n "$jid" ]; then
        running[$jid]="$dir"
        echo "$jid" > "$dir/.jobid"        # 记录 jobid，供断点续跑时查队列
        printf "  [提交] %-10s → JobID=%-8s  (%s)\n" "$dir" "$jid" "$(date '+%H:%M:%S')"
    else
        printf "  [失败] %-10s  %s\n" "$dir" "$out"
        fail_count=$((fail_count + 1))
    fi
}

# ---- 4. 提交初始批次 (最多 MAX_JOBS 个) ----
echo "--- 提交初始批次 ---"
while [ ${#running[@]} -lt $MAX_JOBS ] && [ $next_idx -lt $total ]; do
    submit_one
done
echo ""

# ---- 5. 监控 + 补充提交 ----
echo "--- 开始监控 (Ctrl+C 可中断，已提交的任务不受影响) ---"
echo ""

while [ ${#running[@]} -gt 0 ]; do
    sleep $CHECK_INTERVAL

    # 先快照当前 running 的 key（避免遍历时修改）
    keys=("${!running[@]}")

    for jid in "${keys[@]}"; do
        # 检查 job 是否还在 SLURM 队列中
        status=$(squeue -j "$jid" --noheader 2>/dev/null)
        rc=$?

        # squeue 命令本身出错（如调度器暂时无响应），跳过本次检查
        if [ $rc -ne 0 ]; then
            continue
        fi

        # 输出为空 = job 已离开队列 (完成或失败)
        if [ -z "$status" ]; then
            dir="${running[$jid]}"
            unset running[$jid]
            rm -f "$dir/.jobid"            # 任务已离开队列，清除标记

            if [ -f "$dir/vasprun.xml" ]; then
                printf "  [完成] %-10s (JobID=%-8s) %s  [进度: %d/%d]\n" \
                    "$dir" "$jid" "$(date '+%H:%M:%S')" $((success_count + fail_count + 1)) $total
                success_count=$((success_count + 1))
            else
                printf "  [失败] %-10s (JobID=%-8s) 无 vasprun.xml %s\n" \
                    "$dir" "$jid" "$(date '+%H:%M:%S')"
                fail_count=$((fail_count + 1))
            fi

            # 补充提交下一个任务
            if [ $next_idx -lt $total ]; then
                submit_one
            fi
        fi
    done
done

# ---- 6. 汇总 ----
echo ""
echo "===== 全部完成 ====="
echo "成功: $success_count / $total"
[ $fail_count -gt 0 ] && echo "失败: $fail_count"
echo ""
echo "提取能量:  python3 summary_tru.py   (或 extract_energies.py)"
