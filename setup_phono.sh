#!/bin/bash
echo "===== 开始创建声子计算目录 ====="
missing=0
for f in INCAR POTCAR KPOINTS vasp.sh; do
    [ ! -f "$f" ] && { echo "  [错误] 缺少文件: $f"; missing=$((missing+1)); }
done
[ $missing -gt 0 ] && { echo "请确认以上文件存在。"; exit 1; }
files=$(ls POSCAR-[0-9]* 2>/dev/null | sort -t'-' -k2 -n)
[ -z "$files" ] && { echo "[错误] 未找到 POSCAR-[0-9]*"; exit 1; }
count=0
for poscar in $files; do
    num=$(echo "$poscar" | sed 's/POSCAR-//')
    dir="disp-${num}"
    mkdir -p "$dir"
    cp "$poscar" "$dir/POSCAR"
    ln -sf "$(pwd)/INCAR"   "$dir/INCAR"
    ln -sf "$(pwd)/POTCAR"  "$dir/POTCAR"
    ln -sf "$(pwd)/KPOINTS" "$dir/KPOINTS"
    cp "$(pwd)/vasp.sh" "$dir/"
    echo "  $dir  ←  $poscar"
    count=$((count + 1))
done
echo "===== 共创建 $count 个 dis-* 目录 ====="

