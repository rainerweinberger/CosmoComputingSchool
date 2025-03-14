rm -rf ./cosmo_box_sim

mkdir ./cosmo_box_sim
cd ./cosmo_box_sim
git clone https://github.com/rainerweinberger/arepo.git arepo #clone https://gitlab.mpcdf.mpg.de/vrs/arepo.git
cp ../CosmoComputingSchool/cosmobox/ics.hdf5 ./
cp ../CosmoComputingSchool/cosmobox/param.txt ./
cp ../CosmoComputingSchool/cosmobox/output_list.txt ./
cp ../CosmoComputingSchool/cosmobox/Config.sh ./arepo/
cd arepo
export SYSTYPE="Newton21"
echo $SYSTYPE

module load autotools
module load prun
module load gnu9
module load openmpi4
module load gsl
module load hdf5
module load fftw
module load hypre

make
cd ..
mpiexec -np 4 ./arepo/Arepo param.txt
python ../CosmoComputingSchool/cosmobox/analyze.py ./output/
