## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# Rule to remove everything
rule clean_quality_control:
    shell:
        "rm -fr results/quality_control/"

# Rule to create everything
rule run_quality_control:
    input:
        "results/quality_control/read_qc.done",
    output:
        "results/quality_control/multiqc_report.html"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc results/ -o results/quality_control/ --config workflow/envs/multiqc_config.yaml"


def fastqc_files_output(pep, type_string = "processed"):
    out = []
    for sample in samples(pep):
        if determine_single_end(sample, pep):
            out.append("results/quality_control/fastqc_%s/%s_%s_R0_fastqc.html"%(type_string, sample, type_string))
        else:
            out.append("results/quality_control/fastqc_%s/%s_%s_R1_fastqc.html"%(type_string, sample, type_string))
            out.append("results/quality_control/fastqc_%s/%s_%s_R2_fastqc.html"%(type_string,sample, type_string))
    return out

# General read QC

rule read_qc:
    input:
        fastqc_files_output(pep, "processed"),
        fastqc_files_output(pep, "raw")
    output:
        touch("results/quality_control/read_qc.done")



rule fastqc_raw:
    message: "Running fastqc on raw reads for {wildcards.sample}"
    input:
       lambda wildcards: match_fastq_to_sample(wildcards.sample, wildcards.pair, pep) 
    output:
        "results/quality_control/fastqc_raw/{sample}_raw_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    log:
        stdout="results/quality_control/logs/fastqc_raw/{sample}_{pair}_raw.log",
        stderr="results/quality_control/logs/fastqc_raw/{sample}_{pair}_raw.err"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "zcat {input:q} | fastqc stdin:{wildcards.sample}_raw_{wildcards.pair} "
        "-o results/quality_control/fastqc_raw > {log.stdout} 2> {log.stderr}"

def fastqc_processed_input(sample, pair, pep):
    se = determine_single_end(sample, pep)
    if se and pair == "R0":
        out = "results/preprocessing/trimmomatic/%s_trim_R0.fastq.gz"%(sample)
    elif se:
        raise ValueError("Single end reads do not have an R1 or R2 file")
    else:
        out = "results/preprocessing/trimmomatic/%s_trim_paired_%s.fastq.gz"%(sample, pair)
    return out

rule fastqc_processed:
    message: "Running fastqc on processed reads for {wildcards.sample} {wildcards.pair}"
    input:
        lambda wildcards: fastqc_processed_input(wildcards.sample, wildcards.pair, pep)
    output:
        "results/quality_control/fastqc_processed/{sample}_processed_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    conda:
        "../envs/quality_control.yaml"
    log:
        stdout="results/quality_control/logs/fastqc_processed/{sample}_processed_{pair}.log",
        stderr="results/quality_control/logs/fastqc_processed/{sample}_processed_{pair}.err"
    shell:
        "zcat {input:q} | fastqc stdin:{wildcards.sample}_processed_{wildcards.pair} "
        "-o results/quality_control/fastqc_processed > {log.stdout} 2> {log.stderr}"

