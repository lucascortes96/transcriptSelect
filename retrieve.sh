#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=00:20:00   # walltime
#SBATCH --ntasks=2   # number of tasks
#SBATCH --cpus-per-task=2   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=5G   # memory per node
#SBATCH -J "5primeretrieve"   # job name
#SBATCH -o "5primeretrieve.txt"   # job output file
#SBATCH --mail-user=lucascortes@ebi.ac.uk   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

perl retrieve5prime.pl --diff_limit 10
python3 parseOutput.py --number 10 