## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

rule clean_preprocessing:
    shell:
        "rm -rf results/preprocessing/"

def trimmed_files(pep):
    out = []
    for sample in samples(pep):
        if determine_single_end(sample, pep):
            out.append("results/preprocessing/trimmomatic/%s_trim_R0.fastq.gz"%sample)
        else:
            out.append("results/preprocessing/trimmomatic/%s_trim_paired_R1.fastq.gz"%sample)
            out.append("results/preprocessing/trimmomatic/%s_trim_paired_R2.fastq.gz"%sample)
    return out

rule run_preprocessing:
    input:
        trimmed_files(pep)

rule combine_fastq:
    input:
        lambda wildcards: determine_fastqs_to_combine(wildcards.sample, wildcards.read_num, pep)
    output:
        temp("results/preprocessing/combine_fastq/{sample}_{read_num}_combined.fastq.gz")
    log:
        stderr="results/preprocessing/logs/combine_fastq/{sample}_{read_num}_combined.err"
    threads: 1
    shell:
        "zcat {input} | gzip > {output} 2> {log.stderr}"


rule cutadapt_se:
    input:
        in1=lambda wildcards: match_fastq_to_sample(wildcards.sample, 'R0', pep),
    output:
        out1=temp("results/preprocessing/cutadapt/{sample}_cut_R0.fastq.gz")
    threads: 5
    resources:
        mem_mb=10000
    params:
        cut_param_string = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "cutadapt_se", "cut_param_string"], wildcards.sample, \
        default = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
    log:
        stdout="results/preprocessing/logs/cutadapt_se/{sample}_cutadapt_se.log",
        stderr="results/preprocessing/logs/cutadapt_se/{sample}_cutadapt_se.err"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "cutadapt {params.cut_param_string}  "
        "--cores={threads} "
        "-o {output.out1} {input.in1:q}  > {log.stdout} 2> {log.stderr}"
    
rule cutadapt_pe:
    message: "Running cutadapt on {wildcards.sample}"
    input:
        in1=lambda wildcards: match_fastq_to_sample(wildcards.sample, 'R1', pep),
        in2=lambda wildcards: match_fastq_to_sample(wildcards.sample, 'R2', pep)
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
        "-o {output.out1} -p {output.out2} {input.in1:q} {input.in2:q}> {log.stdout} 2> {log.stderr}"

rule trimmomatic_se:
    message: "Running trimmomatic on {wildcards.sample}"
    input:
        in1="results/preprocessing/cutadapt/{sample}_cut_R0.fastq.gz"
    output:
        out1="results/preprocessing/trimmomatic/{sample}_trim_R0.fastq.gz"
    threads: 1
    resources:
        mem_mb=10000
    params:
        trim_param_string = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "trimmomatic_se", "trim_param_string"], wildcards.sample, default = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15")
    conda:
        "../envs/preprocessing.yaml"
    log:
        stdout="results/preprocessing/logs/trimmomatic_se/{sample}_trim.log",
        stderr="results/preprocessing/logs/trimmomatic_se/{sample}_trim.err"
    shell:
       "trimmomatic SE -phred33 {input.in1} {output.out1} "
       "{params.trim_param_string} > {log.stdout} 2> {log.stderr}"

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
