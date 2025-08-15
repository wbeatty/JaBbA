#!/bin/bash
#SBATCH -J MMRF_SV
#SBATCH -t 10-00:00:00
#SBATCH -p tier3q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000G

export FRAG_COUNTER="/gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library/fragCounter/extdata/frag"

module purge

source ~/miniconda3/etc/profile.d/conda.sh

conda activate snakemake_driver

snakemake --unlock

snakemake --sdm env-modules conda --profile workflow/profiles/slurm/
