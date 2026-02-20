#! /bin/bash
#
#SBATCH --job-name=tsrun_256vecs
#SBATCH --account=accountname
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=normal_q
#SBATCH --constraint=amd
#SBATCH --mem=747G

#
module reset
module load MATLAB/R2024b

## Start MATLAB and call the script
matlab -batch mbf_tsrun

exit 0
