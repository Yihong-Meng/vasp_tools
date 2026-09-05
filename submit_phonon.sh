#!/bin/bash
# ============================================================
# 批量提交声子计算任务 (SLURM, 最多同时运行 4 个)
# 算完一个自动补一个，直到全部完成
# 用法: bash submit_phonon.sh
# 建议: 用 nohup 或 tmux 跑，防止 SSH 断开中断监控
#   nohup bash submit_phonon.sh > submit.log 2>&1 &
#   tail -f submit.log
# ============================================================

MAX_JOBS=4          # 同时运行的最大任务数
CHECK_INTERVAL=30   # 检查间隔(秒)

echo "===== 批量提交声子计算任务 (最多 $MAX_JOBS 个并行) ====="
echo "检查间隔: ${CHECK_INTERVAL}s"
echo ""

# ---- 1. 收集所有 dis-* 目录 ----
mapfile -t all_dirs < <(ls -d disp-* 2>/dev/null | sort)
if [ ${#all_dirs[@]} -eq 0 ]; then
    echo "[错误] 未找到 dis-* 目录，请先运行 setup_phonon.sh"
    exit 1
fi

total=${#all_dirs[@]}
echo "共 $total 个任务待计算"
echo ""

# ---- 2. 状态跟踪 ----
declare -A running    # jobid -> dir
next_idx=0
success_count=0
fail_count=0

submit_one() {
    # 没有更多任务可提交
    if [ $next_idx -ge $total ]; then return 1; fi

    local dir="${all_dirs[$next_idx]}"
    next_idx=$((next_idx + 1))

    # 检查 vasp.sh 是否存在
    if [ ! -f "$dir/vasp.sh" ]; then
        printf "  [跳过] %-10s (缺少 vasp.sh)\n" "$dir"
        fail_count=$((fail_count + 1))
        return 0
    fi

    # sbatch 提交
    local out jid
    out=$(cd "$dir" && sbatch vasp.sh 2>&1)
    jid=$(echo "$out" | grep -oE '[0-9]+$')

    if [ -n "$jid" ]; then
        running[$jid]="$dir"
        printf "  [提交] %-10s → JobID=%-8s  (%s)\n" "$dir" "$jid" "$(date '+%H:%M:%S')"
    else
        printf "  [失败] %-10s  %s\n" "$dir" "$out"
        fail_count=$((fail_count + 1))
    fi
}

# ---- 3. 提交初始批次 (最多 MAX_JOBS 个) ----
echo "--- 提交初始批次 ---"
while [ ${#running[@]} -lt $MAX_JOBS ] && [ $next_idx -lt $total ]; do
    submit_one
done
echo ""

# ---- 4. 监控 + 补充提交 ----
echo "--- 开始监控 (Ctrl+C 可中断，已提交的任务不受影响) ---"
echo ""

while [ ${#running[@]} -gt 0 ]; do
    sleep $CHECK_INTERVAL

    # 先快照当前 running 的 key (避免遍历时修改)
    keys=("${!running[@]}")

    for jid in "${keys[@]}"; do
        # 检查 job 是否还在 SLURM 队列中
        status=$(squeue -j "$jid" --noheader 2>/dev/null)
        rc=$?

        # squeue 命令本身出错 (如调度器暂时无响应)，跳过本次检查
        if [ $rc -ne 0 ]; then
            continue
        fi

        # 输出为空 = job 已离开队列 (完成或失败)
        if [ -z "$status" ]; then
            dir="${running[$jid]}"
            unset running[$jid]

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

# ---- 5. 汇总 ----
echo ""
echo "===== 全部完成 ====="
echo "成功: $success_count / $total"
[ $fail_count -gt 0 ] && echo "失败: $fail_count"
echo ""
echo "提取力常数:  phonopy -f dis-*/vasprun.xml"
echo "生成声子谱:  phonopy -p band.conf"
