Snakemake based pipeline for analyzing short-read DNA sequencing data based on JaBba 

How to Use (on Randi):

First Time Use:
conda env create -f envs/snakemake_base.yaml

Everytime:
conda activate snakemake_base

snakemake --unlock

nohup snakemake --sdm env-modules --profile workflow/profiles/slurm/ > snakemake.out 2>&1 &