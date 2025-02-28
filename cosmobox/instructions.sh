 mkdir ./cosmo_box_sim
 cd ./cosmo_box_sim
 git clone https://gitlab.mpcdf.mpg.de/vrs/arepo.git
 cp ../cosmocomputingschool/cosmobox/ics.hdf5 ./
 cp ../cosmocomputingschool/cosmobox/param.txt ./
 cp ../cosmocomputingschool/cosmobox/Config.sh ./arepo/
 cd arepo
 echo $SYSTYPE
 make