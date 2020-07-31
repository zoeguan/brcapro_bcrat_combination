#!/bin/bash 

#SBATCH -J penmod # A single job name for the array 
#SBATCH -n 12 # cores
#SBATCH -N 1 # nodes
#SBATCH -p serial_requeue # Partition 
#SBATCH -t 0-0:30 # Running time 
#SBATCH --mem 5000 # Memory request 
#SBATCH -o pm.out # Standard output 
#SBATCH -e pm.err # Standard error
#SBATCH --mail-type=END 
#SBATCH --mail-user=zguan@g.harvard.edu

module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER

seed=${SLURM_ARRAY_TASK_ID}
argString="--args "$seed

R CMD BATCH --quiet --no-restore --no-save "$argString" run_penmod_sim.R
