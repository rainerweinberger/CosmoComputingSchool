rm -rf ./galaxy_merger_sim

# the paths assume when typing ls, a directory CosmoComputingSchool shows up 

# create simulation directory, check out code, copy initial conditions and input files needed by the simulation
mkdir ./galaxy_merger_sim
cd ./galaxy_merger_sim
git clone https://github.com/rainerweinberger/gadget4.git gadget4
cp ../CosmoComputingSchool/galaxymerger/1_1_merger_ics_30_15_45_0_rmin10_start320_lr_nogas.dat ./
cp ../CosmoComputingSchool/galaxymerger/param.txt ./
cp ../CosmoComputingSchool/galaxymerger/Config.sh ./gadget4/


# prepare compilation and compile code
cd gadget4
module load autotools
module load prun
module load gnu9
module load openmpi4
module load gsl
module load hdf5
module load fftw
module load hypre

export SYSTYPE="Newton21"
echo $SYSTYPE
make
cd ..

# interactive job on compute nodes (this part can be replaced by jobscript)
srun -p debug --job-name "interactive" --nodes=1 --tasks-per-node=12 --time=02:00:00 --pty bash
module load autotools
module load prun
module load gnu9
module load openmpi4
module load gsl
module load hdf5
module load fftw
module load hypre

# run
mpiexec -np 12 ./gadget4/Gadget4 param.txt

# run analysis on compute node
source  "/opt/aconda3/etc/profile.d/conda.sh"
conda activate py37
python ../CosmoComputingSchool/galaxymerger/analyze.py ./output/

# back to login node
exit
