#!/bin/bash
#SBATCH --job-name=matlab_job      # Job name
#SBATCH --partition=batch-hsw           # Partition name
#SBATCH --time=10:00:00            # Runtime in HH:MM:SS
#SBATCH --ntasks=1                 # Number of tasks (processes)
#SBATCH --cpus-per-task=11         # Number of CPUs per task
#SBATCH --mem=16G                  # Total memory for the job
#SBATCH --output=/scratch/elec/t41011-enhance/LSM/Triton/logs/matlab_job__%j.out # Standard output log
#SBATCH --error=/scratch/elec/t41011-enhance/LSM/Triton/logs/matlab_job__%j.err  # Standard error log

# Ensure the logs directory exists
mkdir -p /scratch/elec/t41011-enhance/LSM/Triton/logs

# Load MATLAB module, adjust the version as necessary
module load matlab/r2021a

# Run the MATLAB script with a command that integrates all steps, ensure to exit MATLAB after execution
matlab -nodisplay -nosplash -nodesktop -r "ExamplePrescript; SpokenDigitsLSM; exit;"
