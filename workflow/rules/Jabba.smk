rule samsort_tumor:
    input:
        lambda wildcards: get_bam_paths(wildcards)["tumor"]
    output:
        temp("resources/bamtofastq/{sample}_tumor_sorted.bam")
    threads: 2
    resources:
        mem_mb = 4000
    benchmark:
        "workflow/report/benchmarks/{sample}/samsort_tumor.txt"
    log:
        "workflow/report/logs/{sample}/samsort_tumor.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -n -@ {threads} -o {output} {input}"

rule samsort_normal:
    input:
        lambda wildcards: get_bam_paths(wildcards)["normal"]
    output:
        temp("resources/bamtofastq/{sample}_normal_sorted.bam")
    threads: 2
    resources:
        mem_mb = 4000
    benchmark:
        "workflow/report/benchmarks/{sample}/samsort_normal.txt"
    log:
        "workflow/report/logs/{sample}/samsort_normal.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -n -@ {threads} -o {output} {input}"

rule bamtofastq:
    input:
        tumor = "resources/bamtofastq/{sample}_tumor_sorted.bam",
        normal = "resources/bamtofastq/{sample}_normal_sorted.bam"
    output:
        tumor_r1 = temp("resources/bamtofastq/fastq/{sample}_tumor_1.fastq"),
        tumor_r2 = temp("resources/bamtofastq/fastq/{sample}_tumor_2.fastq"),
        normal_r1 = temp("resources/bamtofastq/fastq/{sample}_normal_1.fastq"),
        normal_r2 = temp("resources/bamtofastq/fastq/{sample}_normal_2.fastq")
    resources:
        mem_mb = 16000
    benchmark:
        "workflow/report/benchmarks/{sample}/bamtofastq.txt"
    log:
        "workflow/report/logs/{sample}/bamtofastq.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools bamtofastq -i {input.tumor} -fq {output.tumor_r1} -fq2 {output.tumor_r2}
        bedtools bamtofastq -i {input.normal} -fq {output.normal_r1} -fq2 {output.normal_r2}
    if [ -f {output.tumor_r1} ] && [ -f {output.tumor_r2} ] && \
       [ -f {output.normal_r1} ] && [ -f {output.normal_r2} ]; then
        echo "Cleaning input files for {wildcards.sample}..."
        rm -f {input.tumor} {input.normal}
    fi
        """

rule bwa:
    input:
        tumor_r1 = rules.bamtofastq.output.tumor_r1,
        tumor_r2 = rules.bamtofastq.output.tumor_r2,
        normal_r1 = rules.bamtofastq.output.normal_r1,
        normal_r2 = rules.bamtofastq.output.normal_r2
    output:
        tumor_bam = temp("resources/bwa/bams/{sample}_tumor_unprocessed.bam"),
        normal_bam = temp("resources/bwa/bams/{sample}_normal_unprocessed.bam")
    threads: 6
    resources:
        mem_mb = 32000
    benchmark:
        "workflow/report/benchmarks/{sample}/bwa.txt"
    log:
        "workflow/report/logs/{sample}/bwa.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa-mem2 mem -t {threads} resources/bwa_mem2_index/hg19 {input.tumor_r2} | samtools sort -@ {threads} -o {output.tumor_bam}
        bwa-mem2 mem -t {threads} resources/bwa_mem2_index/hg19 {input.normal_r1} {input.normal_r2} | samtools sort -@ {threads} -o {output.normal_bam}
        if [ -f {output.tumor_bam} ] && [ -f {output.normal_bam} ]; then
            echo "Cleaning input files for {wildcards.sample}..."
            rm -f {input.tumor_r1} {input.tumor_r2} {input.normal_r1} {input.normal_r2}
        fi
        """
rule markdup:
    input:
        tumor = rules.bwa.output.tumor_bam,
        normal = rules.bwa.output.normal_bam
    output:
        tumor_clean = "resources/markdup/bams/{sample}_tumor.bam",
        normal_clean = "resources/markdup/bams/{sample}_normal.bam",
        tumor_index = "resources/markdup/bams/{sample}_tumor.bam.bai",
        normal_index = "resources/markdup/bams/{sample}_normal.bam.bai",
        tumor_metrics = "resources/markdup/metrics/{sample}_tumor.metrics",
        normal_metrics = "resources/markdup/metrics/{sample}_normal.metrics"
    threads: 2
    resources:
        mem_mb = 48000
    benchmark:
        "workflow/report/benchmarks/{sample}/markdup.txt"
    log:
        "workflow/report/logs/{sample}/markdup.log"
    conda:
        "../envs/markdup.yaml"
    shell:
        """
(
java -Xmx22G -Dsamjdk.snappy.disable=true -jar /apps/software/openjdk-17.0.2/picard/3.0.0/picard.jar MarkDuplicates \
    INPUT={input.tumor} OUTPUT={output.tumor_clean} METRICS_FILE={output.tumor_metrics} \
    CREATE_INDEX=false MAX_RECORDS_IN_RAM=20000000
) &
(
java -Xmx22G -Dsamjdk.snappy.disable=true -jar /apps/software/openjdk-17.0.2/picard/3.0.0/picard.jar MarkDuplicates \
    INPUT={input.normal} OUTPUT={output.normal_clean} METRICS_FILE={output.normal_metrics} \
    CREATE_INDEX=false MAX_RECORDS_IN_RAM=20000000
) &
wait
samtools index -@ {threads} {output.tumor_clean}
samtools index -@ {threads} {output.normal_clean}
if [ -f {output.tumor_clean} ] && [ -f {output.normal_clean} ] && \
   [ -f {output.tumor_index} ] && [ -f {output.normal_index} ]; then
    echo "Cleaning input files for {wildcards.sample}..."
    rm -f {input.tumor} {input.normal}
fi
        """
rule svaba:
    input:
        tumor = rules.markdup.output.tumor_clean,
        normal = rules.markdup.output.normal_clean,
        tumor_index = rules.markdup.output.tumor_index,
        normal_index = rules.markdup.output.normal_index
    output:
        vcf = "resources/svaba/{sample}/MMRF.svaba.unfiltered.somatic.sv.vcf",
        log = "resources/svaba/{sample}/MMRF.log"
    threads: 4
    resources:
        mem_mb = 10000
    benchmark:
        "workflow/report/benchmarks/{sample}/svaba.txt"
    log:
        "workflow/report/logs/{sample}/svaba.log"
    conda:
        "../envs/svaba.yaml"
    shell:
        """
        mkdir -p resources/svaba/{wildcards.sample}
        svaba run -t {input.tumor} -n {input.normal} -p {threads} \
            -a resources/svaba/{wildcards.sample}/MMRF -G resources/bwa_mem2_index/hg19.fa -v 0 > {output.log} 2>&1
        """
rule vembrane_filter_blacklist:
    input:
        vcf = rules.svaba.output.vcf
    output:
        blacklist = "resources/svaba/{sample}/junctions.nonautosomal.filtered.bnd.vcf"
    benchmark:
        "workflow/report/benchmarks/{sample}/vembrane_filter_blacklist.txt"
    log:
        "workflow/report/logs/{sample}/vembrane_filter_blacklist.log"
    conda:
        "../envs/vembrane_env.yaml"
    shell:
        """
        vembrane filter 'INFO["SVTYPE"] == "BND" and \
        (CHROM in {{"X", "Y", "MT", "GL000192.1", "GL000225.1", "GL000194.1", "GL000193.1", "GL000200.1", "GL000222.1", "GL000212.1", "GL000195.1", "GL000223.1", "GL000224.1", "GL000219.1", "GL000205.1", "GL000215.1", "GL000216.1", "GL000217.1", "GL000199.1", "GL000211.1", "GL000213.1", "GL000220.1", "GL000218.1", "GL000209.1", "GL000221.1", "GL000214.1", "GL000228.1", "GL000227.1", "GL000191.1", "GL000208.1", "GL000198.1", "GL000204.1", "GL000233.1", "GL000237.1", "GL000230.1", "GL000242.1", "GL000243.1", "GL000241.1", "GL000236.1", "GL000240.1", "GL000206.1", "GL000232.1", "GL000234.1", "GL000202.1", "GL000238.1", "GL000244.1", "GL000248.1", "GL000196.1", "GL000249.1", "GL000246.1", "GL000203.1", "GL000197.1", "GL000245.1", "GL000247.1", "GL000201.1", "GL000235.1", "GL000239.1", "GL000210.1", "GL000231.1", "GL000229.1", "GL000226.1", "GL000207.1"}} or \
        bool(re.search(r"(X|Y|MT|GL000\\d+\\.1):", ALT)))' \
        {input.vcf} > {output.blacklist}
        """
rule fragcounter:
    input:
        bam = rules.markdup.output.tumor_clean,
        bai = rules.markdup.output.tumor_index
    output:
        "resources/fragcounter/{sample}/cov.rds"
    resources:
        mem_mb = 6000
    benchmark:
        "workflow/report/benchmarks/{sample}/fragcounter.txt"
    log:
        "workflow/report/logs/{sample}/fragcounter.log"
    envmodules:
        "gcc/12.1.0",
        "R/4.2.1"
    shell:
        """
        mkdir -p resources/fragcounter/{wildcards.sample}
        /gpfs/data/icelake-apps/software/gcc-12.1.0/R/4.2.1/lib64/R/library/fragCounter/extdata/frag \
            -b {input.bam} \
            -d /gpfs/data/drazer-lab/Quinn_JaBbA/gc_map \
            -o resources/fragcounter/{wildcards.sample}
        """
rule dryclean:
    input:
        "resources/fragcounter/{sample}/cov.rds"
    output:
        "resources/dryclean/{sample}/dryclean_cov.rds"
    resources:
        mem_mb = 24000
    benchmark:
        "workflow/report/benchmarks/{sample}/dryclean.txt"
    log:
        "workflow/report/logs/{sample}/dryclean.log"
    envmodules:
        "gcc/12.1.0",
        "R/4.2.1"
    shell:
        """
        mkdir -p resources/dryclean/{wildcards.sample}
        Rscript workflow/scripts/Dryclean.R {input} resources/dryclean/{wildcards.sample}/dryclean_cov.rds
        """
rule jabba:
    input:
        coverage = "resources/dryclean/{sample}/dryclean_cov.rds",
        vcf = "resources/svaba/{sample}/MMRF.svaba.unfiltered.somatic.sv.vcf",
        blacklist = "resources/svaba/{sample}/junctions.nonautosomal.filtered.bnd.vcf"
    output:
        "results/jabba/{sample}/jabba.rds"
    resources:
        mem_mb = 32000
    benchmark:
        "workflow/report/benchmarks/{sample}/jabba.txt"
    log:
        "workflow/report/logs/{sample}/jabba.log"
    envmodules:
        "gcc/12.1.0",
        "go/1.20.1",
        "singularity/3.8.7"
    shell:
        """
    # Only delete markdup outputs if both dryclean and svaba outputs exist
if [ -f resources/dryclean/{wildcards.sample}/dryclean_cov.rds ] && \
   [ -f resources/svaba/{wildcards.sample}/MMRF.svaba.unfiltered.somatic.sv.vcf ] && \
   [ -f resources/svaba/{wildcards.sample}/junctions.nonautosomal.filtered.bnd.vcf ]; then
    echo "Pre-cleaning markdup outputs for {wildcards.sample}..."
    rm -f resources/markdup/bams/{wildcards.sample}_*.bam \
          resources/markdup/bams/{wildcards.sample}_*.bai \
          resources/markdup/metrics/{wildcards.sample}_*.metrics
fi
        mkdir -p results/jabba/{wildcards.sample}
        singularity exec --bind /scratch/wbeatty/:/data /apps/imgs/jabba-v2.sif jba \
            resources/svaba/{wildcards.sample}/MMRF.svaba.unfiltered.somatic.sv.vcf \
            resources/dryclean/{wildcards.sample}/dryclean_cov.rds \
            --field=foreground \
            --ppmethod=sequenza \
            -v \
            --blacklist.coverage=resources/non_autosomes.sorted.bed \
            --blacklist.junctions=resources/{input.blacklist} \
            -o results/jabba/{wildcards.sample}
    if [ -f results/jabba/{wildcards.sample}/jabba.rds ]; then
        echo "Post-cleaning svaba, fragcounter, and dryclean for {wildcards.sample}..."
        rm -rf resources/svaba/{wildcards.sample} \
               resources/fragcounter/{wildcards.sample} \
               resources/dryclean/{wildcards.sample}
    fi
        """
rule cleanup:
    input:
        jabba_rds = "results/jabba/{sample}/jabba.rds"
    output:
        touch("results/jabba/{sample}/cleanup.done")
    log:
        "workflow/report/logs/{sample}/cleanup.log"
    shell:
        """
        echo "Final cleanup pass for {wildcards.sample}..."
        # Remove any lingering intermediate files if they still exist
        rm -f resources/bamtofastq/{wildcards.sample}_tumor_sorted.bam \
              resources/bamtofastq/{wildcards.sample}_normal_sorted.bam \
              resources/bamtofastq/fastq/{wildcards.sample}_*.fastq \
              resources/bwa/bams/{wildcards.sample}_*.bam \
              resources/markdup/bams/{wildcards.sample}_*.bam \
              resources/markdup/bams/{wildcards.sample}_*.bai \
              resources/markdup/metrics/{wildcards.sample}_*.metrics
        rm -rf resources/svaba/{wildcards.sample} \
               resources/fragcounter/{wildcards.sample} \
               resources/dryclean/{wildcards.sample}
	touch results/jabba/{wildcards.sample}/cleanup.done
        """

