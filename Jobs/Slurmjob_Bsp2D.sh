#!/bin/bash
#SBATCH --account=nn9999k 
#SBATCH --job-name=2D-Bsp
#SBATCH --partition=accel --gpus=1
#SBATCH --qos=devel
#SBATCH --time=00:01:00
#SBATCH --mem-per-cpu=2G
#SBATCH -o 2D-Bsp.out

#loading modules
module purge
module load NVHPC/21.2

#read the inputfile "input2D" from /IN
mypath=/cluster/home/hicham/BsplineInterpolation/IN/input2D

srun ./bspline_test-2d.exe<${mypath} 
