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

source ~/.bashrc

echo $@

~/miniconda3/bin/python $@
