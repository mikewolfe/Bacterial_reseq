rule unicycler_illumina:
    message: "Running unicycler with short reads only on {wildcards.sample}"
    input:
        fastq_R1 = "results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        fastq_R2 = "results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz"
    output:
        "results/assembly/unicycler_illumina/{sample}/assembly.fasta"
    threads: 16
    log:
        stdout = "results/assembly/logs/unicycler_illumina/{sample}.log",
        stderr = "results/assembly/logs/unicycler_illumina/{sample}.err"
    conda:
        "../envs/assembly.yaml"
    params:
        out_dir = lambda wildcards: "results/assembly/unicycler_illumina/%s"%(wildcards.sample),
        unicycler_illumina_params = lambda wildcards: lookup_in_config_persample(config, pep, ["assembly", "unicycler", "unicycler_ilummina_params"], wildcards.sample, " ")
    shell:
        "unicycler -1 {input.fastq_R1} -2 {input.fastq_R2} -o {params.out_dir} --threads {threads} "
        "{params.unicycler_illumina_params} > {log.stdout} 2> {log.stderr} "

rule quast:
    message: "checking genome assembly {wildcards.sample} against reference"
    input:
        assembly = "results/assembly/unicycler_illumina/{sample}/assembly.fasta",
        genome = "results/alignment/combine_fasta/{sample}/{sample}.fa",
        features = "results/alignment/combine_bed/{sample}/{sample}.bed"
    output:
        "results/assembly/quast/{sample}/report.txt"
    threads: 16
    log:
        stdout = "results/assembly/logs/quast/{sample}.log",
        stderr = "results/assembly/logs/quast/{sample}.err"
    conda:
        "../envs/assembly.yaml"
    params:
        out_dir = lambda wildcards: "results/assembly/quast/%s"%(wildcards.sample),
        quast_params = lambda wildcards: lookup_in_config_persample(config, pep, ["assembly", "quast", "quast_params"], wildcards.sample, " ")
    shell:
        "quast.py {input.assembly} -r {input.genome} -g {input.features} "
        "{params.quast_params} --threads {threads} "
        "-o {params.out_dir} > {log.stdout} 2> {log.stderr}"
