#!/bin/bash 
#SBATCH -J plink-ukb
#SBATCH -A lindgren.prj 
#SBATCH -o plink-output.out 
#SBATCH -e plink-error.err 
#SBATCH -c 24
#SBATCH --mem-per-cpu 8G
#SBATCH -p short
#SBATCH --array 1-21:1 
#SBATCH --requeue

module load Anaconda3/2022.05
module load OpenBLAS/0.3.8-GCC-9.2.0
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.8-GCC-9.2.0/lib/libopenblas.so
module load java/1.8.0_latest

CHR=${SLURM_ARRAY_TASK_ID}

python3 /well/lindgren/UKBIOBANK/dpalmer/gxe/gxe_hail.py