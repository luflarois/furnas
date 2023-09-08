#!/bin/bash
#SBATCH --job-name=PERA
#SBATCH --nodes=1
#SBATCH --partition=pesq1
#SBATCH --tasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem=160000M
#SBATCH --output=/home/lufla/bin/PERA.txt
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR
module load gnu8
export PATH=/home/public/CEMPA/opt/bin:$PATH
export LD_LIBRARY_PATH=/home/public/CEMPA/opt/lib:$LD_LIBRARY_PATH
export TMPDIR=./tmp
rm -rf tmp/*
ulimit -s unlimited
MPI_PARAMS="-iface ib0 -bind-to core -map-by core"
export OMP_NUM_THREADS=1
time mpirun /home/lufla/era5_stats/build/gfortran_0EA155BAE4FFA839/app/main &> PERA.out

exit
