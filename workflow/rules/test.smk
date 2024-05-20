rule pull_genbank_test:
    message: "Download genbank for {wildcards.accession}"
    output:
        "results/test/references/{accession}.gbk"
    log:
        stdout="results/test/logs/pull_genbank/{accession}.log",
        stderr="results/test/logs/pull_genbank/{accession}.err"
    threads: 1
    conda:
        "../envs/test.yaml"
    shell:
        "ncbi-acc-download {wildcards.accession} --out {output} > {log.stdout} "
        "2> {log.stderr}"

def pull_reference_seqs(sample, pep):
    all_refs = lookup_sample_metadata(sample, "reference_files", pep)
    all_refs = all_refs.rstrip().split(";")
    fasta_names = []
    path_to_fastas = "results/test/process_genbank/"
    for ref in all_refs:
        ref_name = ".".join(ref.split(".")[0:-1])
        fasta_names.append(path_to_fastas + ref_name + ".fna")
    return fasta_names
        

rule combine_fastas_test:
    message: "Generating fasta for sample {wildcards.sample}"
    input:
        input_files = lambda wildcards: pull_reference_seqs(wildcards.sample, pep),
        mask_regions = "test/{sample}.bed"
    output:
        "results/test/combine_fasta/{sample}.fa"
    log:
        stdout="results/test/logs/combine_fasta/{sample}.log",
        stderr="results/test/logs/combine_fasta/{sample}.err"
    threads: 1
    conda:
        "../envs/test.yaml"
    shell:
         "python3 workflow/scripts/combine_fasta.py "
         "results/test/combine_fasta/{wildcards.sample} "
         "{input.input_files} "
         "--masked_regions {input.mask_regions} "
         "--masked_delete "
         "> {log.stdout} 2> {log.stderr}"


rule process_genbank_test:
    message: "Processing genbank for {wildcards.accession}"
    input:
        "results/test/references/{accession}.gbk"
    output:
        "results/test/process_genbank/{accession}.fna"
    log:
        stdout="results/test/logs/process_genbank/{accession}.log",
        stderr="results/test/logs/process_genbank/{accession}.err"
    threads: 1
    conda:
        "../envs/test.yaml"
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} "
         "--outfmt fna "
         " > {output} 2> {log.stderr}"

rule sample_fastqs:
    message: "Sampling for test sample {wildcards.sample}"
    input:
        "results/test/combine_fasta/{sample}.fa"
    output:
        "results/test/input_fastqs/{sample}_R1.fastq.gz",
        "results/test/input_fastqs/{sample}_R2.fastq.gz"
    log:
        stdout="results/test/logs/sample_fastqs/{sample}.log",
        stderr="results/test/logs/sample_fastqs/{sample}.err"
    params:
        sample_param_string = lambda wildcards: lookup_in_config_persample(config,
        pep, ["test", "sample_fastqs", "sample_param_string"], wildcards.sample,
        default = "--rng_seed 42")
    conda:
        "../envs/test.yaml"
    shell:
        "python3 workflow/scripts/FastqSim.py {input} "
        "results/test/input_fastqs/{wildcards.sample} "
        "{params.sample_param_string} "
        "> {log.stdout} 2> {log.stderr}"


rule sample_fastqs_se:
    message: "Sampling for test sample {wildcards.sample}"
    input:
        "results/test/combine_fasta/{sample}.fa"
    output:
        "results/test/input_fastqs/{sample}_R0.fastq.gz",
    log:
        stdout="results/test/logs/sample_fastqs_se/{sample}.log",
        stderr="results/test/logs/sample_fastqs_se/{sample}.err"
    params:
        sample_param_string = lambda wildcards: lookup_in_config_persample(config,
        pep, ["test", "sample_fastqs", "sample_se_param_string"], wildcards.sample,
        default = "--rng_seed 42")
    conda:
        "../envs/test.yaml"
    shell:
        "python3 workflow/scripts/FastqSim.py {input} "
        "results/test/input_fastqs/{wildcards.sample} "
        "{params.sample_param_string} "
        "--single_end "
        "> {log.stdout} 2> {log.stderr}"
        
