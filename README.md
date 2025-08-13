Snakemake based pipeline for analyzing short-read DNA sequencing data based on JaBba 

How to Use (on Randi):

conda env create -f envs/snakemake-base.yaml

conda activate snakemake_base

snakemake --unlock (if needed)

snakemake --sdm conda env-modules --profile workflow/profiles/slurm/ 

nohup snakemake --sdm conda env-modules --profile workflow/profiles/slurm/ --until "fragcounter" > snakemake.out 2>&1 &