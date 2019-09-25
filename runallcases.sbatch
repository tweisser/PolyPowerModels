#! /bin/bash
#SBATCH --job-name=allcases_dense
#SBATCH --nodes=1
#SBATCH --exclusive=user
#SBATCH --time=2-0:00:00
#SBATCH --no-requeue
#SBATCH --output=out_job
#SBATCH --error=err_job
#SBATCH --qos=long
set -x                          # Output commands
set -e                          # Abort on errors
echo "Checking:"
pwd
hostname
date
env | sort > ENVIRONMENT
echo "Starting:"
# ---------------------------------------------
julia --project=. run_acopf.jl
# ---------------------------------------------
echo "Stopping:"
date
echo "Done."
