
# the paths assume when typing ls, a directory CosmoComputingSchool shows up 

# create simulation directory, check out code, copy initial conditions and input files needed by the simulation
mkdir ./cosmo_box_sim
cd ./cosmo_box_sim
git clone https://github.com/rainerweinberger/arepo.git arepo 
cp ../CosmoComputingSchool/cosmobox/ics.hdf5 ./
cp ../CosmoComputingSchool/cosmobox/param.txt ./
cp ../CosmoComputingSchool/cosmobox/output_list.txt ./
cp ../CosmoComputingSchool/cosmobox/Config.sh ./arepo/


# prepare compilation and compile code
cd arepo
module load autotools
module load prun
module load gnu9
module load openmpi4
module load gsl
module load hdf5
module load fftw
module load hypre

source  "/opt/aconda3/etc/profile.d/conda.sh"
conda activate py37

export SYSTYPE="Newton21"
echo $SYSTYPE
make
cd ..

# submit via slurm job
cp ../CosmoComputingSchool/cosmobox/run.slurm ./
sbatch run.slurm

# run analysis on login node
python ../CosmoComputingSchool/cosmobox/analyze.py ./output/

# back to login node
exit
