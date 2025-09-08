#! /bin/bash
#
#SBATCH --job-name=csweepICsParll
#SBATCH --account=compsci
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal_q
#SBATCH --constraint=amd
#SBATCH --mem=747G

#
module reset
module load MATLAB/R2024b

## Start MATLAB and call the script
matlab -batch mbf_csweepICsParll

exit 0