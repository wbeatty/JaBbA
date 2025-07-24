def get_bam_paths(wildcards):
    for pair in sample_pairs:
        if pair["sample_name"] == wildcards.sample:
            return {
                "tumor": f"{BAM_DIR}/{pair['tumor_bam']}",
                "normal": f"{BAM_DIR}/{pair['normal_bam']}"
            }
    raise ValueError(f"No BAM paths found for sample {wildcards.sample}")