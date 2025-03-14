rm -rf ./cosmo_box_sim

# create simulation directory, check out code, copy initial conditions and input files needed by the simulation
mkdir ./cosmo_box_sim
cd ./cosmo_box_sim
git clone https://github.com/rainerweinberger/arepo.git arepo #clone https://gitlab.mpcdf.mpg.de/vrs/arepo.git
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
mpiexec -np 4 ./arepo/Arepo param.txt

# run analysis on compute node
python ../CosmoComputingSchool/cosmobox/analyze.py ./output/

# back to login node
exit
