# Drazer Lab JaBbA
Snakemake based pipeline for analyzing short-read DNA sequencing data with JaBbA on the University of Chicago's RANDI 

## Usage

### If this is the first time using the pipeline in your scratch space:

1. Navigate to your scratch space and clone the repo (you may be asked to sign in):
```
git clone https://github.com/wbeatty/JaBbA.git
```

2. Navigate into the created folder:
```
cd JaBbA
```

3. Create base environment:
```
conda env create -f workflow/envs/snakemake_base.yaml
```

### How to Load Sample Pairs:

1. Activate R:
```
module purge && module load gcc/12.1.0 && module load R/4.2.1
```

2. Ensure Sample Sheet is in correct format:
> File Name, Sample ID, Case ID, Tissue Type (these are the only required columns, Case ID column can be left blank if NA)

3. Copy sample sheet from wherever it is located to cluster and into resources/samples directory

4. Run pair generation (insert sheet name, e.g: Test_samples.xlsx):
```
Rscript resources/samples/BAMSampleSheet.R <sample sheet name>
```

5. Samples are automatically loaded and ready to be executed in the pipeline (can double check pairs with `nano resources/samples/sample_pairs.py`)

### Running the Pipeline:

1. Before running, purge modules and activate the conda environment:
```
module purge && conda activate snakemake_base
```

2. Ensure BAM directory is set (located in config/config.yaml):

`nano config/config.yaml`, update BAM_DIR with path to files

3. Run the pipeline
```
sbatch run_analysis.sh
```

You may want to periodically check on the job and your scratch space storage over the next few days. You can run: `squeue -u <username>`
This would allow you to see how long the job has been running—if it had started running and is not currently pending. To get more information on the pipeline itself, in your scratch space, run: `cat slurm-<jobID>.out` (remove the <> for both commands)

You can Control+F to search for “finished job” and then you can see the rough percentage and most recent steps completed by flipping through that match search. Keep in mind that it is not linear with time so the percents won’t compute linearly (i.e., 1 day for 10% does not mean it will take 10 days to run)


### Questions or Issues
- Feel free to reach out to me via email at wbeatty@umich.edu or wbeatty@protonmail.com
- Currently there are still some issues with the fragcounter function which may cause errors in output

### Citations

> **Most large structural variants in cancer genomes can be detected without long reads.**
> Choo, ZN., Behr, J.M., Deshpande, A. et al.
>
> _Nat Genet_ 2023 Nov 09. doi: [https://doi.org/10.1038/s41588-023-01540-6](https://doi.org/10.1038/s41588-023-01540-6)

> **Distinct Classes of Complex Structural Variation Uncovered across Thousands of Cancer Genome Graphs**
> Hadi, Kevin et al.
>
> Cell, Volume 183, Issue 1. doi: [https://doi.org/10.1016/j.cell.2020.08.006](https://doi.org/10.1016/j.cell.2020.08.006)