#!/bin/bash
#SBATCH --job-name=blur3D
#SBATCH --output blur3D_dyn.out
#SBATCH --error blur3D_dyn.err
#SBATCH -N 1
#SBATCH -n 128


module purge
#module load intel/oneapi/mpi/2021.4.0
module load slurm/16.05.8
module load gcc/9.4.0 

module list

hostname

echo "------------------"
echo
# echo "Job working directory: $SLURM_SUBMIT_DIR"
echo "Current working directory: $(pwd)"
echo

num=128
echo "Total processes: $num"
echo

echo "Job starting at $(date)"
echo

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=128

./blur3Do 1>out 2>log

echo
echo "Job finished at $(date)"
