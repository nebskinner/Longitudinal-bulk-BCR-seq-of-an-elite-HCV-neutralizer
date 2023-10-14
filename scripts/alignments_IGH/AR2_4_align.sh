#!/bin/bash
#SBATCH --partition=himem
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=75GB
#SBATCH --time=7-00:00:00

set -e

module purge

ml load GCC/9.3.0
ml load OpenMPI/4.0.3
ml load R/4.2.2

R CMD BATCH /home/gdskinnerlab/nes002/serial_timepoints_C110/scripts/AR2_4_align.R > AR24_align.Rout
