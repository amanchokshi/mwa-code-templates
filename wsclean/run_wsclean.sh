#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60000
#SBATCH --account=oz048
#SBATCH --time=00:15:00

module load gcc/6.4.0 openmpi/3.0.0
module load fftw/3.3.7
module load gsl/2.4
module load cfitsio/3.420
module load boost/1.67.0-python-2.7.14
module load hdf5/1.10.1
module load openblas/0.2.20
module load cuda/9.0.176

cd /fred/oz048/achokshi/wsclean_III/clean_output


#time /fred/oz048/MWA/CODE/bin/wsclean -name uv_model_test \
#    -size 2048 2048 -auto-threshold 0.5 -auto-mask 3 -multiscale \
#    -niter 0 -mgain 0.85 -weight briggs 0 \
#    -small-inversion -pol I -channels-out 2 -j 12 \
#    -join-channels -fit-spectral-pol 1 -scale 0.005 \
#    /fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/CHIPS_input/1061315448/*.ms

#time /fred/oz048/MWA/CODE/bin/wsclean -name rephased_model_test \
#    -size 2048 2048 -auto-threshold 0.5 -auto-mask 3 -multiscale \
#    -niter 0 -mgain 0.85 -weight briggs 0 \
#    -small-inversion -pol I -channels-out 2 -j 12 \
#    -join-channels -fit-spectral-pol 1 -scale 0.005 \
#    /fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/CHIPS_input/1061315448/rephased_uv_model*.ms

time /fred/oz048/MWA/CODE/bin/wsclean -name 1125953008_wsclean \
    -size 2048 2048 -auto-threshold 0.5 -auto-mask 3 -multiscale \
    -niter 0 -mgain 0.85 -weight briggs 0 \
    -small-inversion -pol I -channels-out 24 -j 12 \
    -join-channels -fit-spectral-pol 1 -scale 0.005 \
    /fred/oz048/achokshi/oskar_III/1125953008/three_point_source_chips-t08_f0.080_band0*.ms




#time /fred/oz048/MWA/CODE/bin/wsclean -name phase1_VLA-ForA+GLEAM \
#    -size 4096 4096 -auto-threshold 0.5 -auto-mask 3 -multiscale \
#    -niter 1000000 -mgain 0.85 -save-source-list -weight briggs 0 \
#    -small-inversion -pol I -channels-out 8 -j 24 \
#    -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data \
#    -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
#    -join-channels -fit-spectral-pol 1 -scale 0.004 \
#    /fred/oz048/jline/ForA_OSKAR/data/1102864528/*.ms \
#    /fred/oz048/jline/ForA_OSKAR/data/1102865128/*.ms \
#    /fred/oz048/jline/ForA_OSKAR/data/1102865728/*.ms

