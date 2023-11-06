#!/bin/bash --login
#SBATCH --nodes=1 # number of nodes
#SBATCH --mem=30G # memory pool for all cores

#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR
#SBATCH --ntasks-per-node=1
#SBATCH --partition=spot-vhmem
#SBATCH --mail-type=ALL
#SBATCH --signal=SIGUSR1@90
#SBATCH --time=UNLIMITED

module load matlab/r2022a

# cd DMOEA

# don't use -singleCompThread in parallel computing

matlab -nodisplay  -r "try main 10 10; catch; end; quit"
