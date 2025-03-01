rm -rf ./cosmo_box_sim

mkdir ./cosmo_box_sim
cd ./cosmo_box_sim
git clone git@bitbucket.org:RWeinberger/rw_arepo.git arepo #clone https://gitlab.mpcdf.mpg.de/vrs/arepo.git
cp ../CosmoComputingSchool/cosmobox/ics.hdf5 ./
cp ../CosmoComputingSchool/cosmobox/param.txt ./
cp ../CosmoComputingSchool/cosmobox/output_list.txt ./
cp ../CosmoComputingSchool/cosmobox/Config.sh ./arepo/
cd arepo
export SYSTYPE="macOShomebrew"
echo $SYSTYPE
make
cd ..
mpiexec -np 4 ./arepo/Arepo param.txt
# python analyze.py
