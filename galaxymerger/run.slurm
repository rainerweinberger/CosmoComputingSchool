#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=debug
#SBATCH -t 02:00:00
#SBATCH --job-name=my_galaxymerger #[name of simulation]
##SBATCH --mail-user=[if you want to get notified, set your email hear and remove #]
##SBATCH --mail-type=BEGIN
##SBATCH --mail-type=FAIL
##SBATCH --mail-type=END

module load autotools
module load prun
module load gnu9
module load openmpi4
module load gsl
module load hdf5
module load fftw
module load hypre

wdir=$SLURM_SUBMIT_DIR
job_id=$SLURM_JOBID

mkdir ./logs
prun ./gadget4/Gadget4 param.txt 1> $wdir/logs/job.$job_id.out 2> $wdir/logs/job.$job_id.err
