#!/bin/bash 

#Set number of coils per period 
coils_per_p=10
#set number of  periods to plot 
periods=5
tot_coils=$((coils_per_p * periods))
thickness=0.2 # dr and dz (dr=dz, coils assume approx square cross section)

#make sure all geo files are removed 
rm 'S_coil'* 
#Ensure all coordinate files cleared 
rm coords/''*
#Generate coil coordinates
python3 3D_coil_gen.py --regcoil_path '../examples/regcoil_out.W7X.nc' --nescin_path '../examples/nescin.out' --coils_per_p $coils_per_p --p $periods --c_t $thickness
#Call code that makes coil step files here
python3 Coil_ind_splines_pythonAPI.py --coils $tot_coils #$tot_coils








