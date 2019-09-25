#! /bin/bash
#SBATCH --job-name=VLS_A21000
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
julia --project=. test/test.jl
# ---------------------------------------------
echo "Stopping:"
date
echo "Done."
