#!/bin/bash
# ============================================================
# vasp_ab.sh — VASP 晶格常数扫描模板（六角结构）
# ============================================================
# 功能：对六角结构（a=b）进行 a 轴扫描，b 轴按六角关系自动跟随，
#       c 轴固定（适用于 2D 层状材料）。
#
# 使用说明：
#   1. 检查并修改下方 "用户必须修改" 部分的路径。
#   2. 根据你的集群环境修改 #SBATCH 参数（分区、节点数、时间等）。
#   3. 运行：sbatch vasp_ab.sh
# ============================================================

# ==================== SLURM 作业参数 ====================
# ⚠️ 请根据你的集群配置修改以下参数
#SBATCH --partition=your_partition_name       # 分区名称
#SBATCH -N 1                                  # 节点数
#SBATCH --nodes=2                             # 节点数（与 -N 保留一个即可）
#SBATCH --ntasks-per-node=48                  # 每节点核心数
#SBATCH --time=99:00:00                       # 运行时间限制
#SBATCH -J job_name                           # 作业名称
#SBATCH -o %j.out                             # 标准输出文件
#SBATCH -e %j.err                             # 标准错误文件
#SBATCH --qos your_qos_name                   # QoS 名称（如需要）

# ==================== 环境加载（用户必须修改） ====================
# 方案一：使用 module load（推荐，适用于大多数超算中心）
# module load intel/your_version
# module load vasp/your_version

# 方案二：使用 Intel oneAPI 的 setvars.sh（需修改为你的实际路径）
# source /path/to/intel/oneapi/setvars.sh

# 方案三：直接加载各个组件（需修改为你的实际路径）
# source /path/to/intel/compiler/env/vars.sh intel64
# source /path/to/intel/mkl/env/vars.sh intel64
# source /path/to/intel/mpi/env/vars.sh intel64

# ============================================================

echo "Running on hosts "
echo "The job ID is: ${SLURM_JOB_ID}"
echo "The job's name is: ${SLURM_JOB_NAME}"
echo "Directory is: ${PWD}"
echo "This job runs on the following nodes: ${SLURM_JOB_NODELIST}"
echo "This job has allocated ${SLURM_JOB_CPUS_PER_NODE} cpu cores."

echo -n "start time  " > time
date >> time

# ==================== VASP 执行路径（用户必须修改） ====================
# 如果 vasp_std 已在系统 PATH 中，直接使用 vasp_std
VASP_CMD="vasp_std"

# 如果 vasp_std 在特定路径，取消注释并修改为实际路径
# VASP_CMD="/path/to/your/vasp_installation/bin/vasp_std"

# 如果使用环境变量（推荐），在 ~/.bashrc 中设置 VASP_BIN 后使用
# VASP_CMD="${VASP_BIN}"

# ============================================================

echo "Starting Time is `date`"

# ==================== 计算参数配置 ====================
# 六角结构：a=b，c 轴固定（2D 材料真空层）
C_FIX=24.0  # c 轴晶格常数（Å）

# a 轴扫描范围（可根据需要修改）
for a in 3.60 3.65 3.70 3.728 3.75 3.80 3.85 3.90
do

    # 用 awk 计算六角结构 b 轴分量
    # a 轴: ($a, 0, 0)
    # b 轴: (-a/2, a*sqrt(3)/2, 0)
    bx=$(awk -v a=$a 'BEGIN{printf "%.16f", -a/2.0}')
    by=$(awk -v a=$a 'BEGIN{printf "%.16f", a*sqrt(3.0)/2.0}')

    # 生成 POSCAR
    cat > POSCAR <<!
New structure
   1.00000000000000
     $a   0.0000000000000000    0.0000000000000000
     $bx   $by   0.0000000000000000
     0.0000000000000000    0.0000000000000000   $C_FIX
   XX   XX   XX
     2     1     1
Direct
  0.3333333429999996  0.6666666870000029  0.4650317275885359
  0.6666666570000004  0.3333333129999971  0.5303604201082806
  0.0000000000000000  0.0000000000000000  0.4006600477315120
  0.0000000000000000  0.0000000000000000  0.6039478045716719
!

    # 运行 VASP
    mpirun -np $SLURM_NPROCS $VASP_CMD >> log

    # 提取结果
    P=$(grep "press" OUTCAR | tail -1 | awk '{printf "%12.9f \n", $4 }')
    V=$(grep "volume" OUTCAR | tail -1 | awk '{printf "%12.9f \n", $5 }')
    E=$(grep "TOTEN" OUTCAR | tail -1 | awk '{printf "%12.9f \n", $5 }')

    # 备份文件
    cp CONTCAR POSCAR-$a
    cp OSZICAR OSZICAR-$a
    cp OUTCAR OUTCAR-$a
    echo $a $E >> summary.dat

done

echo "Ending Time is `date`"
echo -n "end   time  " >> time
date >> time