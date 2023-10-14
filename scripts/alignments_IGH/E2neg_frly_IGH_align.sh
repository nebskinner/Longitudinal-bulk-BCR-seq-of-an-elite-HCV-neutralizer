#!/bin/bash
#SBATCH --partition=himem
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00

set -e

module purge

ml load GCC/9.3.0
ml load OpenMPI/4.0.3
ml load R/4.2.2

R CMD BATCH /home/gdskinnerlab/nes002/serial_timepoints_C110/scripts/E2neg_frly_IGH_align.R > E2neg_frly_IGH_align.Rout
