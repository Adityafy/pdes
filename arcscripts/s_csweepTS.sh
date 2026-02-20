#! /bin/bash
#
#SBATCH --job-name=csweepTS
#SBATCH --account=accountname
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=normal_q
#SBATCH --constraint=amd
#SBATCH --mem=256G

#
module reset
module load MATLAB/R2024b

## Start MATLAB and call the script
matlab -batch mbf_csweepTS

exit 0