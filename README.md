[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.bitbucket.io) [![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
# Pipeline for processing bacterial re-sequencing data

The goal of this project is to create a reproducible analysis pipeline
for bacterial re-sequencing data that can be easily deployed in any
computational environment. 

# What the pipeline does 
Starting from raw fastq files this pipeline does the following:

- preprocessing to remove adapters and low quality sequences using
  [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/) and
  [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)
- variant calling against a reference genome using
  [`breseq`](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing)

In addition this pipeline uses [`multiqc`](https://multiqc.info/) to
compile the following quality control metrics into an interactive html
report:
- Read QC using
  [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  both before and after preprocessing

# Installing the pipeline

This pipeline has three dependencies:
- The package manager
  [`miniconda`](https://docs.conda.io/en/latest/miniconda.html)
- The workflow management system
  [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html)
- The API for working with Portable Encaspulated Projects (PEP)
  [`peppy`](http://peppy.databio.org/en/latest/)

`miniconda` can be installed following the installation instructions
for your system
[here](https://docs.conda.io/en/latest/miniconda.html).

Once `miniconda` is installed, both `snakemake` and `peppy` can be
installed in their own environment easily using: 

``` 
conda create -n Bacterial_reseq snakemake>=8 peppy 
conda activate Bacterial_reseq
```

Now you can pull the pipeline from GitHub using: 
``` 
git clone https://github.com/mikewolfe/Bacterial_reseq 
```

And you can change into the newly cloned `Bacterial_reseq` directory
and test your installation with: 
``` 
snakemake --use-conda --cores 10
```

This will create a small test data set consisting of example fastqs
sampled uniformly from genomes where rpoA rpoB and rpoC have been
deleted and then run the pipeline on those test datasets.

The first time you run the pipeline it will need to create dedicated
`conda` environments for each module which will take some time.
Afterwards, it will run in about 20-30 minutes, most of which is time
needed to create the test data. For more information on using `conda`
with `snakemake` including how to set things up to run offline check
out the documentation
[here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management).

If everything runs smoothly you can then clean up and remove the
results from the test data using:
```
snakemake clean_all --cores 1
```

## If installing on a job-managed computational cluster (like HT-condor)

If you are using a computational cluster that requires job management
software, you will want to install the required job management software with
your environment as well.  For example, if you are using an
`htcondor`-managed server you would instead create your environment
like so: 

``` 
conda create -n Bacterial_reseq snakemake>=8.0 peppy htcondor>=8.9.5 snakemake-executor-plugin-cluster-generic
conda activate Bacterial_reseq 
```

And then run the pipeline using the job management software
with an environment-specific
[profile](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles).

For running with `htcondor` I have included an example profile in
[htcondor_profile/](htcondor_profile/). This is based on the
work-in-progress v8.0 profile
[here](https://github.com/Snakemake-Profiles/htcondor).

To get this setup do the following:
- Create the directory `~/.config/snakemake/htcondor`
- Copy the files in [htcondor_profile/](htcondor_profile/) over to
  `~/.config/snakemake/htcondor`
- Create the directory `~/.condor_jobs`

Then you should be able to run the pipeline with `snakemake` version
8.0 and automatic submission of jobs through `htcondor` using:

``` 
snakemake --use-conda --cores 10 --profile htcondor
```


## If installing on an ARM-based Mac (i.e. 2020 M1 Chips and beyond)

Many of the bioinformatic packages in this pipeline still do not have ARM-based
binaries available for download through `conda`. Luckily, there is an easy
workaround.

Apple has `Rosetta` software that will translate instruction sets from
programs compiled for intel to the AMD processor. First you need to install
`Rosetta` to your computer. Instructions are
[here](https://support.apple.com/en-us/102527).

Next you need to set `conda` to look for intel-based programs rather than AMD
based programs. This is a simple as adding the following line to your `.condarc`
in your home directory: `subdir: osx-64`. Your overall `.condarc` should now
look like this:

```
channel_priority: 'strict'
channels:
  - conda-forge
  - bioconda
subdir: osx-64
```

Now install the running environment as before:

``` 
conda create -n Bacterial_reseq snakemake>=5.24.2 peppy htcondor>=8.9.5 
conda activate Bacterial_reseq 
```

And run the pipeline as described below.

# Running the pipeline

This pipeline uses `snakemake` to manage the workflow and familiarity
with `snakemake` will help with getting the most out of the pipeline.
Fortunately, `snakemake` has an excellent tutorial that can be found
[here](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
if you are unfamiliar with this workflow management system.

## Input

This pipeline takes as input a [Portable Encaspulated Project
(PEP)](http://pep.databio.org/en/latest/) which is essentially
a `.csv` of samples and metadata together with a `.yaml` file allowing
for extensions to the existing metadata.

The following required fields are needed for this pipeline
- `sample_name` - a unique identifier for each sample
- `filenameR1` - the base file name for read1 of a set of paired fastq
  files
- `filenameR2` - the base file name for read2 of a set of paired fastq
  files
- `infile_path` - the path to where the files for a given sample live
- `reference_files` - a semicolon separated list of files for the
  reference genome. Can be genbank or fasta format
- `reference_path` - the path to where the reference files for a given
  sample live.

An example of a sample sheet is included at
[pep/test_samples.csv](pep/test_samples.csv).

Additionally, the sample sheet can be augmented with a required
`config.yaml` file. In the included test example this is used to
replace the `file_path` field with a specific location. This example
can be found at [pep/config.yaml](pep/config.yaml).

The `pep/config.yaml` should be edited to point towards your
`pep/samples.csv` file. If you create your `samples.csv` file with
excel be sure to save it as `comma separated values` not any of the
other encodings for `.csv`.

## Configuration

The pipeline itself, including parameters controlling specific tools,
is controlled by a `.yaml` file in
[config/config.yaml](config/config.yaml). The included `config.yaml`
has all possible options specified with comments describing what those
options control.

## Rules

The pipeline is organized into modules each of which runs a specific
task needed for analysis.

- [workflow/rules/preprocessing.smk](workflow/rules/preprocessing.smk)
  includes rules for trimming raw reads for adapters and quality
- [workflow/rules/variant_calling.smk](workflow/rules/variant_calling.smk)
  includes rules for calling changes from the reference files
- [workflow/rules/quality_control.smk](workflow/rules/quality_control.smk)
  includes rules for performing summarizing quality control on the
  reads themselves
- [workflow/rules/assembly.smk](workflow/rules/assembly.smk)
  includes rules for doing de novo assembly with `unicycler` and report
  generation with `quast`
- [workflow/rules/test.smk](workflow/rules/test.smk)
  includes rules for creating the small test set.

Each of these rules can be run individually using:
```
snakemake run_module_name --use-conda --cores 10
```

For example:
```
snakemake run_preprocessing --use-conda --cores 10
```

Additionally to remove the output of a given module run:
```
snakemake clean_module_name --use-conda --cores 1
```

For example:
```
snakemake clean_preprocessing --use-conda --cores 1
```

Many of later modules are dependent on the earlier modules and running
a later module will run the required rules in an earlier module
automatically.

# Issues with the pipeline

If you run into any issues with the pipeline and would like help
please submit it to the Issues page. Please include your
`config/config.yaml` file, your `pep/config.yaml` file, your
`pep/samples.yaml` file, and the output from `snakemake` that includes
your error.

# Version history

Currently at version 0.0.2

See the [Changelog](CHANGELOG.md) for version history and upcoming
features.
