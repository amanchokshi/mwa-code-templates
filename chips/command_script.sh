#! /bin/bash
#############################################################################
#Slurm jobscript for running a single ObsID.  Second level program for 
#running firstpass on Oscar. First level program is 
#batch_firstpass.sh
#############################################################################


/fred/oz048/achokshi/mwa_dipole/full_sim/chips/run_CHIPS.py \
    --obs_list=/fred/oz048/achokshi/mwa_dipole/obslist/1125953248.txt \
    --data_dir=/fred/oz048/achokshi/mwa_dipole/full_sim/oskar/ \
    --uvfits_dir='/' \
    --uvfits_tag='1125953248_chips-t08_f0.080_band' \
    --output_tag=AC_1125953248_full \
    --band=high \
    --obs_range=0,1 \
    --no_delete \
    --field=1 \
    --timeres=8.0 \
    --base_freq=167.035e+6 \
    --freqres=80000
