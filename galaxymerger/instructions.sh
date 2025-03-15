rm -rf ./galaxy_merger_sim

mkdir ./galaxy_merger_sim
cd ./galaxy_merger_sim
git clone https://github.com/rainerweinberger/gadget4.git gadget4
cp ../CosmoComputingSchool/galaxymerger/1_1_merger_ics_30_15_45_0_rmin10_start320_lr_nogas.dat ./
cp ../CosmoComputingSchool/galaxymerger/param.txt ./
cp ../CosmoComputingSchool/galaxymerger/Config.sh ./gadget4/
cd gadget4
export SYSTYPE="macOShomebrew"
echo $SYSTYPE
make
cd ..
mpiexec -np 4 ./gadget4/Gadget4 param.txt
python ../CosmoComputingSchool/galaxymerger/analyze.py ./output/
