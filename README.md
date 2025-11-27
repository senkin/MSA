# Mutational Signature Attribution pipeline

![logo](MSA.png)

Mutational signature attribution analysis, featuring automated optimisation study with simulated data.

## Introduction
The purpose of this README is to provide a guide for the quick start of using MSA. An extensive Wiki page detailing the usage of this tool can be found [here](https://gitlab.com/s.senkin/MSA/-/wikis/home).

## Running with Nextflow
The best way to run the code is by using [Nextflow](https://www.nextflow.io/).
Once you have installed Nextflow, run the test job locally or on your favourite cluster:

```bash
nextflow run https://gitlab.com/s.senkin/MSA -profile conda,test
```

**Important**: The pipeline has been updated to Nextflow DSL2 and requires Nextflow version 21.04.0 or later. The version can be changed using the following command:

```bash
export NXF_VER=25.04.8
```

It is recommended to use [docker](https://www.docker.com/) or [singularity](https://sylabs.io/singularity/) profiles as these normally provide greater stability and reproducibility than [conda](https://conda.io).

The pipeline should run and produce all the results automatically. You can also retrieve the code ([see below](#getting-started)) in order to adjust all the inputs and parameters. In the [nextflow.config](nextflow.config) file various parameters can be specified.

## Running on SigProfiler output

MSA natively supports [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) and [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) outputs.

The simplest way to run is as follows:

```bash
nextflow run https://gitlab.com/s.senkin/MSA -profile conda \
    --dataset SP_test \
    --SP_extractor_output_path /full/path/to/SP_extractor_output/
```

If [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) output is provided, it will take priority over the [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) one for input mutation matrices:

```bash
nextflow run https://gitlab.com/s.senkin/MSA -profile conda \
    --dataset SP_test \
    --SP_matrix_generator_output_path /full/path/to/SP_matrix_generator_output/ \
    --SP_extractor_output_path /full/path/to/SP_extractor_output/
```

## Running with specific input files

MSA now supports direct specification of individual mutation tables and signature files:

```bash
nextflow run https://gitlab.com/s.senkin/MSA -profile conda \
    --dataset my_data \
    --input_mutation_table /path/to/mutations.txt \
    --signatures_file /path/to/signatures.txt \
    --mutation_types SBS \
    --SBS_context 96
```

**Note**: When using specific files, you must specify exactly ONE mutation type and (for SBS) ONE context.

## Options

All parameters are described in the dedicated [wiki page](https://gitlab.com/s.senkin/MSA/-/wikis/Parameters-description-table). Most general parameters are listed below.

### Input Priority

The pipeline processes inputs in the following priority order:
1. **Specific files** (highest): `--input_mutation_table` and `--signatures_file`
2. **SigProfiler outputs**: `--SP_extractor_output_path` and `--SP_matrix_generator_output_path`
3. **Default directories** (lowest): `--input_tables` and `--signature_tables`

### General parameters

| Parameter | Default value | Description |
|-----------|---------------|-------------|
| --help | null | Print usage and optional parameters |
| --input_mutation_table | null | Specific mutation table file to convert and use (overrides all other inputs) |
| --signatures_file | null | Specific signature file to convert and use (overrides all other inputs) |
| --SP_matrix_generator_output_path | null | Use SigProfilerMatrixGenerator output from specified **full path** |
| --SP_extractor_output_path | null | Use SigProfilerExtractor output from specified **full path** to attribute signatures |
| --dataset | SIM_test | Set the name of the dataset |
| --input_tables | $baseDir/input_mutation_tables | **Full path** to input mutation tables directory |
| --signature_tables | $baseDir/signature_tables | **Full path** to input signature tables directory |
| --signature_prefix | sigProfiler | Prefix of signature files (e.g. sigProfiler, sigRandom) |
| --output_path | . | Output path for plots and tables |
| --temp_path | $output_path/temp | Temporary path for converted inputs |
| --mutation_types | ['SBS'] | Mutation types to analyse. Only ONE can be specified from command line; use list in script for multiple |
| --number_of_samples | -1 | Number of samples to analyse (-1 means all available) |
| --SBS_context | 96 | SBS context to use (96, 192, 288, 1536, or 4608) |
| --COSMIC_signatures | false | If true, use COSMIC signatures from SigProfiler output; otherwise use de-novo |
<!-- | --signatures_to_prioritise | null | set a list of signatures to prioritise in penalty calculation (only these will be used to calculate optimal penalty) |
| --signatures_to_deprioritise | null | set a list of signatures to deprioritise (these will be excluded in penalty calculation) |
| --no_CI_for_penalties | false | do not use confidence intervals for optimal penalties calculation | -->

### Output structure

The pipeline organizes outputs as follows:

```
output_path/
├── output_tables/              # Final NNLS results
│   └── {dataset}/
│       ├── output_*_mutations_table.csv
│       ├── output_*_weights_table.csv
│       ├── signatures_prevalences_*.csv
│       └── bootstrap_output/
├── output_tables_unoptimised/  # Unoptimized results
├── plots/                      # All generated plots
└── temp/                       # Temporary conversions
    ├── input_tables/
    └── signature_tables/
```

## Running manually

### Getting started

Retrieve the code:
```bash
git clone https://gitlab.com/s.senkin/MSA.git
cd MSA
```

Run the fully automated pipeline with optimisation on the test dataset:
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda,test \
    --output_path test
```

Alternatively, to run the pipeline without optimisation, using fixed penalties (DSL1 only, please downgrade Nextflow with `export NXF_VER=22.10.4`):
```bash
nextflow run run_analysis.nf -profile conda \
    --dataset SIM_test \
    --weak_threshold 0.02 \
    --output_path test
```

### Using resume

The pipeline partially supports Nextflow's `-resume` functionality to restart from cached processes:

```bash
nextflow run run_auto_optimised_analysis.nf -profile conda -resume \
    --dataset my_data \
    --output_path results
```

### Setting up dependencies

If you cannot run *nextflow*, you can still run some basic analysis manually (scripts in the *./bin* folder).
Dependencies: *pandas*, *numpy*, *scipy*, *matplotlib* and *seaborn*. 

Set up the virtual environment using [conda](https://conda.io):

```bash
conda env create -f environment.yml
```

This only needs to be done once. Afterwards, activate the environment:

```bash
conda activate msa
```

Alternatively, you can use *docker* with the provided *Dockerfile*, or use ready-made images ([docker](https://hub.docker.com/r/ssenkin/msa/tags) or [singularity](https://cloud.sylabs.io/library/ssenkin/default/msa)).

### Simulating data

* [input_mutation_tables/SIM](input_mutation_tables/SIM) folder contains a set produced with existing (PCAWG or COSMIC) signatures using normal distributions of mutational burdens.

* To reproduce the simulated set of samples with reshuffled *SBS1/5/22/40* PCAWG signatures:
```bash
python bin/simulate_data.py -t SBS -c 96 -n 100 -s signature_tables
```

* [input_mutation_tables/SIMrand](input_mutation_tables/SIMrand) folder contains simulated samples where each contains contributions from **5** randomly selected signatures out of **100** Poisson-generated signatures. To reproduce:
```bash
python bin/generate_random_signatures.py -t SBS -c 96 -n 100
python bin/simulate_data.py -r -t SBS -c 96 -n 100 -d SIMrand
```

Both scripts support normal distribution noise using `-z` option with standard deviation set by `-Z` option (**2** by default). Additional flags: use `-h` option.

## Example workflows

### Basic analysis with repository test data
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda
```

### Analysis with SigProfiler extractor output
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda \
    --dataset my_cohort \
    --SP_extractor_output_path /data/sigprofiler_output \
    --output_path results
```

### Analysis with specific mutation table
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda \
    --dataset my_sample \
    --input_mutation_table data/my_mutations_SBS96.txt \
    --mutation_types SBS \
    --SBS_context 96 \
    --output_path results
```

### Multiple mutation types (must be set in script)
Edit `run_auto_optimised_analysis.nf`:
```groovy
params.mutation_types = ['SBS', 'DBS', 'ID']
```
Then run:
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda \
    --dataset my_data \
    --output_path results
```

## Changes from DSL1

Key improvements in the DSL2 version:
- Modern Nextflow syntax (DSL2)
- Improved caching and resume functionality
- Support for specific file inputs
- Cleaner module organization
- All temporary files in configurable temp directory
- Better error handling and validation
- Unified output directory structure

## Citation

If you use this pipeline, please cite the original publication and the MSA repository.
