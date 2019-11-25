#!/bin/sh
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
##SBATCH --qos=1day
#
# Request it to run this for HH:MM:SS with ?G per core
#
##SBATCH --time=23:59:00
#SBATCH --time=05:59:00

source /scicore/home/neher/neher/miniconda3/etc/profile.d/conda.sh
conda activate

echo $@

~/miniconda3/bin/python $@
