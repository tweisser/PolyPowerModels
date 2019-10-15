#! /bin/bash
#SBATCH --job-name=30fsr_api
#SBATCH --nodes=1
#SBATCH --exclusive=user
#SBATCH --time=2-0:00:00
#SBATCH --no-requeue
#SBATCH --output=sbatch_out/out_deg4_30fsr_api
#SBATCH --error=sbatch_out/err_deg4_30fsr_api
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
julia --project=. run_deg4.jl
# ---------------------------------------------
echo "Stopping:"
date
echo "Done."
