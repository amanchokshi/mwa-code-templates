from subprocess import call
import argparse
from numpy import *
from os import getcwd
from os.path import exists
from sys import exit

parser = argparse.ArgumentParser(description='Generates an sbatch script that takes 24 RTS uvdump*.uvfits files,\
all concats them into one single. Can be used to convert multiple observations with one controlling script. \
Uses fixmwams to add keywords to concatenated ms so wsclean can correct for primary beam')

parser.add_argument('--uvfits_location', default=False,help='If only converting one observation, \
put the full path to the uvdump*uvfits files. Use this OR a combination of --base_dir, --obs_list and --uvfits_loc_tag')

parser.add_argument('--base_dir',help='If converting multiple obs IDs, look in this base directory for \
the observations specfied in --obs_list, and then within each obs, look for uvfits files \
inside --uvfits_loc_tag (so searches for /base_dir/obs/uvfits_loc_tag/uvdump*.uvfits)')

parser.add_argument('--obs_list',help='List of observations to cycle through - use with --base_dir and --uvfits_loc_tag')
parser.add_argument('--uvfits_loc_tag', help='Name of the dir inside the obs ID dir where the uvdump files live')
# parser.add_argument('--output_tag', help='Name for concatinated measurement set, default="concat_obs.ms"',default='concat_obs.ms')
# parser.add_argument('--metafits', help='Metafits file for fixmwams to use')
parser.add_argument('--metafits_loc', default='/fred/oz048/MWA/data/',
      help='Location of metafits files - default /fred/oz048/MWA/data/*obs*/*obs*__metafits_ppds.fits')
# parser.add_argument('--no_delete_coarse', help='By default, delete',default=True,action='store_false')
parser.add_argument('--file_prepend', help='Prepend of file name - defaults to the RTS standare uvdump_',default='uvdump_')

args = parser.parse_args()

# output_tag = args.output_tag
# if output_tag[-3:] == '.ms':
#     pass
# else:
#     output_tag += '.ms'

prepend = args.file_prepend

def test_24_files(base_dir,obs,uvfits_loc_tag):
    for band in range(1,25):
        filename = '%s/%s/%s/%s%02d.uvfits' %(base_dir,obs,uvfits_loc_tag,prepend,band)
        if exists(filename):
            all_good = True
        else:
            exit('Could not find file %s\nCheck your inputs for --base_dir,--obs_list,and --uvfits_loc_tag\nExiting' %filename)

    return '%s/%s/%s/' %(base_dir,obs,uvfits_loc_tag)

def test_24_files_smaller(prepend):
    for band in range(1,25):
        filename = '%s%02d.uvfits' %(prepend,band)
        if exists(filename):
            all_good = True
        else:
            exit('Could not find file %s\nCheck your inputs for --uvfits_location. \nExiting' %filename)


try:
    obs = [line for line in open('%s' %args.obs_list,'r').read().split('\n') if line != '' and '#' not in line]
except:
    exit('Cannot open --obs_list', args.obs_list, 'exiting now')


convert_locs = []

if args.uvfits_location:
    print('====================================================================')
    print('User has defined --uvfits_location. Only converting ONE observation. \
Do not use --uvfits_location if processing multiple observations (see --help)')
    print('====================================================================')
    test_24_files_smaller(args.uvfits_location)
    convert_locs.append(args.uvfits_location)

else:
    for ob in obs:
        convert_loc = test_24_files(args.base_dir,ob,args.uvfits_loc_tag)
        convert_locs.append(convert_loc)
        # convert_scripts.append(convert_scipt)


control = open('run_convert-rts2ms.sh','w+')

control.write("#!/bin/bash -l\n")
control.write("#SBATCH --nodes=1\n")
control.write("#SBATCH --cpus-per-task=1\n")
control.write("#SBATCH --mem=20000\n")
control.write("#SBATCH --account=oz048\n")
control.write('#SBATCH --array=1-24\n')

##allow 25 minutes for each obs being converted
hours = floor(len(convert_locs) * (10.0 / 60.0))
mins = (len(convert_locs)*10) % 60

control.write("#SBATCH --time=%02d:%02d:00\n" %(int(hours),int(mins)))
control.write('module load gcc/6.4.0 openmpi/3.0.0\n')
control.write('module load fftw/3.3.7\n')
control.write('module load gsl/2.4\n')
control.write('module load cfitsio/3.420\n')
control.write('module load boost/1.67.0-python-2.7.14\n')
control.write('module load hdf5/1.10.1\n')
control.write('module load openblas/0.2.20\n')

control.write('BANDNUM=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")\n')


control.write("cd %s\n" %getcwd())


for ind,convert_loc in enumerate(convert_locs):
    if args.uvfits_location:

        cmd = 'time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa '
        cmd += '--nologger -c /fred/oz048/jline/ForA_OSKAR/convert_uv2ms.py '
        cmd += '%s $BANDNUM  \n' %(convert_loc)
        control.write(cmd)

        control.write('/fred/oz048/MWA/CODE/bin/fixmwams %s$BANDNUM.ms %s/%s/%s_metafits_ppds.fits\n' %(convert_loc,args.metafits_loc,obs[ind],obs[ind]))
        cmd = 'time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa '
        cmd += '--nologger -c /fred/oz048/jline/ForA_OSKAR/fix_uv2ms.py '
        cmd += '%s $BANDNUM \n' %(convert_loc)
        control.write(cmd)


        #control.write('time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa \
        #               --nologger -c /fred/oz048/jline/ForA_OSKAR/convert_uv2ms.py \
        #               %s $BANDNUM  \n' %(convert_loc))
        #control.write('/fred/oz048/MWA/CODE/bin/fixmwams %s$BANDNUM.ms %s/%s/%s_metafits_ppds.fits\n' %(convert_loc,args.metafits_loc,obs[ind],obs[ind]))
        #control.write('time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa \
        #               --nologger -c /fred/oz048/jline/ForA_OSKAR/fix_uv2ms.py \
        #               %s $BANDNUM  \n' %(convert_loc))
    #
    else:
        cmd = 'time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa '
        cmd += '--nologger -c /fred/oz048/jline/ForA_OSKAR/convert_uv2ms.py '
        cmd += '%s/%s $BANDNUM  \n' %(convert_loc,prepend)
        control.write(cmd)

        control.write('/fred/oz048/MWA/CODE/bin/fixmwams %s/%s$BANDNUM.ms %s/%s/%s_metafits_ppds.fits\n' %(convert_loc,prepend,args.metafits_loc,obs[ind],obs[ind]))
        cmd = 'time /fred/oz048/MWA/CODE/casa-release-5.3.0-143.el7/bin/casa '
        cmd += '--nologger -c /fred/oz048/jline/ForA_OSKAR/fix_uv2ms.py '
        cmd += '%s/%s $BANDNUM \n' %(convert_loc,prepend)
        control.write(cmd)

control.close()

print('Wrote controlling slurm sci:pt run_convert-rts2ms.sh . Finsihed!')
