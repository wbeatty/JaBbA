#!/bin/bash
#SBATCH -J MMRF_SV
#SBATCH -t 7-00:00:00
#SBATCH -p tier3q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake_driver

which snakemake
snakemake --version

snakemake --unlock

snakemake --sdm conda env-modules --profile workflow/profiles/slurm/ 
