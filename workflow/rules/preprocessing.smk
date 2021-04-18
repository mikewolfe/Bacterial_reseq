## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

rule clean_preprocessing:
    shell:
        "rm -rf results/preprocessing/"

rule run_preprocessing:
    input:
        expand("results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz", sample = samples(pep))

def get_raw_fastas_per_sample(sample, pep):
    in1 = lookup_sample_metadata(sample, "infile_path", pep) + \
            lookup_sample_metadata(sample, 'filenameR1', pep)
    in2 = lookup_sample_metadata(sample, "infile_path", pep) + \
            lookup_sample_metadata(sample, 'filenameR2', pep)
    return [in1, in2]
    
rule cutadapt_pe:
    message: "Running cutadapt on {wildcards.sample}"
    input:
        lambda wildcards: get_raw_fastas_per_sample(wildcards.sample, pep)
    output:
        out1=temp("results/preprocessing/cutadapt/{sample}_cut_R1.fastq.gz"),
        out2=temp("results/preprocessing/cutadapt/{sample}_cut_R2.fastq.gz")
    threads: 5
    resources:
        mem_mb=10000
    params:
        cut_param_string = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "cutadapt_pe", "cut_param_string"], wildcards.sample, \
        default = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
    log:
        stdout="results/preprocessing/logs/cutadapt/{sample}_cutadapt.log",
        stderr="results/preprocessing/logs/cutadapt/{sample}_cutadapt.err"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "cutadapt {params.cut_param_string} "
        "--cores={threads} "
        "-o {output.out1} -p {output.out2} {input:q} > {log.stdout} 2> {log.stderr}"

rule trimmomatic_pe:
    message: "Running trimmomatic on {wildcards.sample}"
    input:
        in1="results/preprocessing/cutadapt/{sample}_cut_R1.fastq.gz",
        in2="results/preprocessing/cutadapt/{sample}_cut_R2.fastq.gz"
    output:
        out1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        out2=temp("results/preprocessing/trimmomatic/{sample}_trim_unpaired_R1.fastq.gz"),
        out3="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        out4=temp("results/preprocessing/trimmomatic/{sample}_trim_unpaired_R2.fastq.gz")
    threads: 1
    resources:
        mem_mb=10000
    params:
        trim_param_string = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "trimmomatic_pe", "trim_param_string"], wildcards.sample, default = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15")
    log:
        stdout="results/preprocessing/logs/trimmomatic/{sample}_trim.log",
        stderr="results/preprocessing/logs/trimmomatic/{sample}_trim.err"
    conda:
        "../envs/preprocessing.yaml"
    shell:
       "trimmomatic PE -phred33 {input.in1} {input.in2} {output.out1} {output.out2} "
       "{output.out3} {output.out4} {params.trim_param_string} > {log.stdout} 2> {log.stderr}"
