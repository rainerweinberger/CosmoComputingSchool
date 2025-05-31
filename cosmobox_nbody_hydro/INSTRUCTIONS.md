
# COSMOBOX

In this example, we are simulating stucture formation in the Universe. The initial conditions are created with using the matter power spectrum at redshift z=127 and the Zel'dovich approximation to calculate the initial displacement and velocities of otherwise uniformly distributed simulation particles.

## Step-by-step instructions

1. create directory, e.g. ``` mkdir ./cosmo_box_sim ```
2. change into this directory ``` cd ./cosmo_box_sim ```
3. download the Arepo code via git ``` git clone https://gitlab.mpcdf.mpg.de/vrs/arepo.git ```
4. copy files ics.hdf5, Config.sh to simulation directory 
``` 
cp ../CosmoComputingSchool/cosmobox/ics.hdf5 ./
cp ../CosmoComputingSchool/cosmobox/param.txt ./
```
5. copy file Config.sh to code directory ``` cp ../CosmoComputingSchool/cosmobox/Config.sh ./arepo/ ```


6. change to arepo directory ``` cd arepo ```
7. check SYSTYPE ``` echo $SYSTYPE ```
8. compile Arepo ``` make ```
9. 


