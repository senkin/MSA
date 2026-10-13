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

**Important**: The pipeline uses Nextflow DSL2 and is tested with Nextflow 26.04.6 (`nextflowVersion = '>=25.04.0'` in the manifest). Pin the version for reproducible runs:

```bash
export NXF_VER=26.04.6
```

**Strict syntax**: the pipeline follows Nextflow's strict syntax (`nextflow lint run_auto_optimised_analysis.nf modules nextflow.config` reports no errors), so it runs with the strict parser, the default from Nextflow 26.04, as well as with the legacy one (`NXF_SYNTAX_PARSER=v1`, the default of earlier releases); no parser setting is needed.
With Nextflow 25.04 or 25.10, keep their default (legacy) parser: their strict parser (`NXF_SYNTAX_PARSER=v2`) rejects `-profile` together with an extra `-c` file, e.g. `-c gpu.config` ("Unknown configuration profile", in any pipeline; fixed in Nextflow 26.04).

It is recommended to use [docker](https://www.docker.com/) or [singularity](https://sylabs.io/singularity/) profiles as these normally provide greater stability and reproducibility than [conda](https://conda.io).

**On a cluster**, choose your executor and queue in your own config file, and see [Running on a cluster](#running-on-a-cluster) for the resources each task needs.

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
| --help | false | Print usage and optional parameters |
| --input_mutation_table | null | Specific mutation table file to convert and use (overrides all other inputs) |
| --signatures_file | null | Specific signature file to convert and use (overrides all other inputs) |
| --SP_matrix_generator_output_path | null | Use SigProfilerMatrixGenerator output from the specified path |
| --SP_extractor_output_path | null | Use SigProfilerExtractor output from the specified path to attribute signatures |
| --dataset | SIM_test | Set the name of the dataset |
| --input_tables | $projectDir/input_mutation_tables | Input mutation tables directory: the tables of a dataset are read from `<input_tables>/<dataset>/` |
| --signature_tables | $projectDir/signature_tables | Input signature tables directory |
| --signature_prefix | sigProfiler | Prefix of signature files (e.g. sigProfiler, sigRandom) |
| --output_path | . | Output path for plots and tables (absolute, or relative to the launch directory) |
| --cleanup_temp | true | `true`: only results are published to the output path; `false`: intermediate files too (see [Output structure](#output-structure)) |
| --mutation_types | ['SBS'] | Mutation types to analyse: a list in the config, or comma-separated on the command line (e.g. `--mutation_types SBS,DBS,ID`) |
| --number_of_samples | -1 | Number of samples to analyse (-1 means all available) |
| --SBS_context | 96 | SBS context to use (96, 192, 288, 1536, or 4608) |
| --COSMIC_signatures | false | If true, use COSMIC signatures from SigProfiler output; otherwise use de-novo |
| --seed | null | RNG seed for the simulations and bootstrap resampling. If unset, a seed is derived from the Nextflow session ID; the seed used is always recorded in `pipeline_info/` (see [Reproducibility](#reproducibility)) |
| --bootstrap_chunk_size | 0 | Bootstrap iterations per job. 0: one per job on CPU, all in one job on GPU. Results do not depend on it |
| --bootstrap_method | binomial | Bootstrap resampling: `binomial`, `multinomial`, `multinomial_weight`, `poisson`, `classic` (mutation channels resampled with replacement) or `bootstrap_residuals` (residuals of the central fit resampled across channels and added back to it) |
| --optimisation_strategy | removal | `removal`: start from all signatures and remove the weak ones; `addition`: start from SBS1 and SBS5 (those in the signature table, otherwise its first signature) and add the strong ones; `add-remove`: from the same start, alternate addition and removal. See [Penalty scan](#penalty-scan) |
| --weak_thresholds | auto | Weak penalties scanned in the optimisation: `auto`, or an explicit list (config list, or comma-separated on the command line, e.g. `0,0.0001,0.0002,0.0003`). See [Penalty scan](#penalty-scan) |
| --strong_thresholds | [0] | Strong penalties scanned (similarity increase needed to add a signature), used by `addition` and `add-remove`; same formats as `--weak_thresholds`, except `auto` |
| --auto_weak_thresholds_number | 12 | Number of non-zero penalties chosen by `auto` (0 is always included) |
| --auto_weak_thresholds_samples | 300 | Simulated samples whose removal trajectories choose the `auto` penalties (-1: all) |
| --threshold_scan | trajectory | `trajectory`: one greedy removal per sample, replayed at every penalty (identical results); `separate`: one removal per penalty, as before |
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
│       ├── CIs_*.csv, signatures_prevalences_*.csv, pruned_attribution_*.csv
│       ├── truth_studies/      # simulated data only (with a truth table)
│       ├── *.json              # bootstrap distributions (*)
│       └── bootstrap_output/   # bootstrap iterations (*)
├── outputs_optimisation/
│   ├── SIM_{dataset}/{type}/   # penalty scan: thresholds, metrics, optimal penalties
│   └── SIM_{dataset}_{context}_NNLS_{weak}_{strong}/  # per-penalty attributions and bootstrap (*)
├── optimisation_plots/
├── plots/                      # All generated plots
├── pipeline_info/              # run_info_<run name>.txt: seed, code revision, versions, parameters
└── temp/                       # (*) converted inputs, unoptimised attributions, simulated data
    ├── input_tables/
    ├── signature_tables/
    ├── output_tables_unoptimised/
    └── output_tables/
```
(*) intermediate files, published only with `--cleanup_temp false`. Each task reads its inputs
from its own work directory (staged through channels) and the pipeline never reads back from the
output path, so results do not depend on it being shared with the tasks, and `-resume` reuses
everything cached in the work directory.

## Running on a cluster

Keep site-specific settings in your own config file, passed with `-c`, for example for SLURM:

```groovy
// cluster.config
process {
    executor = 'slurm'                                       // or sge, lsf, pbspro, ...
    queue = 'your-partition'
    resourceLimits = [cpus: 8, memory: 64.GB, time: 24.h]   // your cluster's largest job
}
executor {
    queueSize = 500                                          // jobs queued at once
}
```

and run with e.g. `nextflow run https://gitlab.com/s.senkin/MSA -profile apptainer -c cluster.config --dataset ...`.

The numbers below come from a typical run: one mutation type, 1000 samples, the 65 SBS signatures of the
default reference set and default settings (11,000 simulated samples, 100 bootstrap iterations in the
optimisation and 1000 in the final attribution).

* **Tasks.** Each bootstrap iteration is a task of its own (`--bootstrap_chunk_size 0`, the default), so
  such a run has about 1,130 tasks: 100 for the optimisation bootstrap, 1000 for the final one and about
  30 others. Mutation types run in parallel. The `test` profile has only 60 tasks for three mutation
  types, as it uses 10 iterations in chunks of 5.
* **CPUs.** Each task requests one CPU and uses one. The attribution, the tables and the plots are
  single-threaded, and the maths libraries are limited to one thread (`OMP_NUM_THREADS`,
  `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` in [nextflow.config](nextflow.config)); without that,
  NumPy's OpenBLAS starts a thread per core of the node. More CPUs per task would stay idle: their
  threads do not speed up these small problems.
* **Memory.** The pipeline requests memory per process, and retries a task killed for exceeding it
  (exit status 130-145, e.g. 137) with twice, then three times the memory:

  | Processes | Requested | Peak in the typical run | Grows with |
  |---|---|---|---|
  | `make_bootstrap_tables` (final attribution) | 8 GB | 3.8 GB | samples × iterations × signatures |
  | `make_optimisation_bootstrap_tables` | 4 GB | 2.2 GB | simulated samples × iterations × signatures |
  | `plot_bootstrap_attributions` | 4 GB | 2.8 GB | samples × iterations × signatures |
  | GPU jobs (`--use_GPU true`) | 8 GB | (not measured) | |
  | all other processes | 2 GB | 0.15-0.45 GB | samples |

  For a larger cohort, raise the first rows in your config, e.g. for 5000 samples and 1000 final
  iterations: `process { withName: 'make_bootstrap_tables' { memory = 24.GB } }`. Use `withName` or
  `withLabel` (`memory_high`, `memory_medium`) for these: the pipeline's own setting for them takes
  precedence over a plain `process.memory` in your config.
* **Time.** On one 2.8 GHz Xeon core, the attribution takes 0.8-0.9 s per sample and iteration with
  65 signatures, and roughly scales with the square of the number of signatures (a SigProfiler
  decomposition of 10-20 signatures is 10-40 times faster). In the typical run that is 2.5-3 h for
  `run_optimized_NNLS` and for each of the 100 optimisation bootstrap tasks (11,000 simulated samples),
  and about 15 min for each final bootstrap task. The input spectra plots take about 1 h, alongside the
  rest; the bootstrap plots, at the end, about 3 h (`--plot_bootstrap_attributions false` skips them).
  The other tasks take minutes. Set time limits accordingly, e.g.
  `process { withName: 'run_optimized_NNLS|run_bootstrap_NNLS' { time = 8.h } }`.
* **Retries.** When a retried task then succeeds, Nextflow's summary line still counts the failed
  attempt (`[FAILED] completed=... failed=1`), but the run's exit status is 0.
* **On one machine.** The local executor refuses a task that requests more memory than is free
  ("Process requirement exceeds available memory"). Cap the requests with e.g.
  `-c local.config` holding `process.resourceLimits = [memory: 8.GB]`; the `test` profile caps them at 2 GB.

## GPU acceleration (experimental)

The greedy NNLS signature-optimisation hot path can be offloaded to the GPU via
[cuML](https://github.com/rapidsai/cuml)'s batched NNLS solver. This is opt-in and does not
affect the default CPU behaviour.

1. **Build the GPU environment** (once per cluster + GPU architecture, as a build job on a
   node with an NVIDIA driver — not the login node). This builds the source-only cuML into a
   conda env and tops it up with MSA's dependencies:

   ```bash
   ./setup_gpu_env.sh                              # CUDA 13.3 (default)
   ./setup_gpu_env.sh -c all_cuda-129_arch-x86_64 -e all_cuda-129_arch-x86_64   # if the GPU driver only supports CUDA 12.x
   ```

   Match the CUDA version to the driver on your GPU nodes (`nvidia-smi` shows the maximum
   supported CUDA). The script prints the full env path to use in the next step.

   cuML is built from a pinned commit of the `fea-nnls` branch (`CUML_COMMIT` in the script;
   `--commit` or `--unpinned` to change it), and the build is recorded in the env:
   `etc/msa_cuml_build.json` (commit, options, versions — copied into every GPU result as
   `gpu_provenance.json` / the task log) and `etc/msa_gpu_env.lock.txt` (every conda package).
   The cuML build env uses RAPIDS nightly packages, which are removed after a while, so keep
   that lock: `./setup_gpu_env.sh --lock msa_gpu_env.lock.txt` rebuilds the same env elsewhere.

   Before relying on a new env, validate it on a GPU node against the CPU (fp64) reference.
   This checks GPU fp64 against CPU, repeatability, independence from the batch size, and the
   default `mixed` precision against fp64 (a small fraction of attributions may legitimately
   differ: greedy decisions that sit within rounding error of the threshold):

   ```bash
   tests/gpu_acceptance.sh --gpu-python /path/to/gpu-env/bin/python --cpu-python /path/to/msa-env/bin/python
   ```

2. **Enable it** in [nextflow.config](nextflow.config) (or on the command line). `gpu_conda_env`
   must be the **full path** to the built env — a bare name is treated by Nextflow's `conda`
   directive as a package to install:

   ```groovy
   use_GPU       = true
   gpu_conda_env = '/home/you/miniforge3/envs/all_cuda-133_arch-x86_64'
   ```

3. **Route the NNLS jobs to a GPU queue.** Cluster-specific queue/GPU-request settings are
   deliberately not committed; put them in your own config (e.g. `gpu.config`) and pass it with
   `-c gpu.config` only for GPU runs. See [gpu.config.example](gpu.config.example) for a template
   scoped to the `gpu_nnls`-labelled processes.

Then run with a conda-enabled profile:

```bash
nextflow run run_auto_optimised_analysis.nf -profile conda,test_GPU -c gpu.config
```

The penalty scan replays one batched removal trajectory per sample on the GPU as well; the
acceptance check above also verifies that this gives exactly the results of one run per penalty
on your GPU (`--threshold_scan separate` is the fallback).

By default each GPU bootstrap runs as a single job. To spread it over several GPUs, set
`--bootstrap_chunk_size` (e.g. `250` runs 1000 iterations as 4 GPU jobs; `test_GPU` does this).
Iteration `i` is always seeded from `(seed, i)`, so the split does not change the results.

## Penalty scan

The optimisation attributes the simulated data at a range of weak penalties (`weak_thresholds`)
and picks the one that best meets the chosen metric (`metric_to_prioritise`, `metric_threshold`).

* **One pass for all penalties.** In the greedy removal (default `optimisation_strategy = "removal"`)
  the penalty only decides when to stop; which signature is removed at each step does not depend
  on it. So the removal is run once per sample down to a single signature, recording the similarity
  drop at each step, and replayed at every penalty. Results are identical to one run per penalty
  (on CPU and in the batched GPU path), and one job covers all penalties, so scanning many penalties
  costs little extra computation. `--threshold_scan separate` restores one run per penalty.
* **Automatic range (`weak_thresholds = 'auto'`, default).** A sample's attributed signatures only
  change at the running maxima of its recorded similarity drops. Before the scan, these are
  collected from the trajectories of `auto_weak_thresholds_samples` simulated samples, and the
  penalties are 0 plus `auto_weak_thresholds_number` values log-spaced between their 5th and 95th
  percentiles (values below 1e-6, numerical noise, are ignored). The range therefore adapts to the
  cohort and signature set instead of a fixed `[0, 0.0001, 0.0002, 0.0003]`, which can be far below
  the optimum for some cohorts.
* **Explicit ranges** work as before, e.g. `weak_thresholds = [0, 0.0001, 0.0002, 0.0003]` in the
  config or `--weak_thresholds 0,0.0001,0.0002,0.0003`, and give the same results as before.
  `auto` requires the removal strategy; other strategies need an explicit list.
* **Addition strategies.** `optimisation_strategy = "addition"` starts from SBS1 and SBS5 (those of
  them in the signature table, otherwise its first signature) and adds, one at a time, the signature
  that increases the similarity most, while the increase exceeds the strong penalty. `"add-remove"`
  starts from the same signatures and alternates addition and removal until the similarity no longer
  changes. Scan strong penalties with `strong_thresholds` (e.g. `--strong_thresholds 0.001,0.005,0.01`);
  `addition` ignores the weak penalty, so a single weak value (e.g. `--weak_thresholds 0`) avoids
  repeating identical runs.

The penalties used are in `outputs_optimisation/SIM_<dataset>/<mutation type>/weak_thresholds.txt`.
Each penalty still adds its own bootstrap output files and summary-table job, so very long lists
increase disk usage and post-processing time. The per-penalty post-processing (final fit, metrics,
output files) runs on the CPU host even in GPU mode, so on a GPU it, rather than the NNLS, dominates
the optimisation jobs; set `bootstrap_chunk_size` to spread the optimisation bootstrap over several
jobs. If the selected penalty is the largest one scanned, a warning suggests extending the range.

## Reproducibility

* **Random seed.** One seed drives the simulations and the bootstrap resampling. Without
  `--seed`, a seed is derived from the Nextflow session ID (new for each fresh run, unchanged
  on `-resume`). The seed actually used is logged and written, with the code revision, Nextflow
  version, profile and all parameters, to `<output_path>/pipeline_info/run_info_<run name>.txt`;
  passing it back with `--seed` reproduces the run.
* **Software.** `-profile conda` installs the exact lock of [environment.yml](environment.yml) for
  your platform from [conda-lock/](conda-lock) (`--conda_env environment.yml` for an unpinned
  solve); after editing `environment.yml`, regenerate the locks with `conda-lock/update_locks.sh`.
  `-profile docker/singularity/apptainer` use the MSA image pinned by digest
  (`params.container_image`). Pin Nextflow with `export NXF_VER=26.04.6`.
* **Tests.** `python -m unittest discover -s tests -v` (CPU only, a few minutes) checks that the
  GPU code path gives the same attributions as the CPU path (including MSA's cuML wrapper, on a
  NumPy stand-in for the GPU), that simulations and bootstrap chunks are reproducible, every
  bootstrap method on every mutation type and SBS context, simulations from inputs with empty
  samples, the bootstrap summary tables and the penalty selection.
  To compare the attributions of any two runs, e.g. CPU vs GPU, use
  `python bin/compare_attributions.py <reference dir> <candidate dir>`.

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

Alternatively, to run the pipeline without optimisation, using fixed penalties (DSL1 only, now under `legacy/`; please downgrade Nextflow with `export NXF_VER=22.10.4`):
```bash
nextflow run legacy/run_analysis.nf -profile conda \
    --dataset SIM_test \
    --weak_threshold 0.02 \
    --output_path test
```

### Using resume

Nextflow's `-resume` restarts a run from its cached tasks: unchanged tasks are not run again, and
their results are published again (also if the output folder was deleted in between):

```bash
nextflow run run_auto_optimised_analysis.nf -profile conda -resume \
    --dataset my_data \
    --output_path results
```

### Setting up dependencies

If you cannot run *nextflow*, you can still run some basic analysis manually (scripts in the *./bin* folder).
Dependencies: *pandas*, *numpy*, *scipy*, *numba*, *matplotlib*, *seaborn* and *statsmodels*. 

Set up the virtual environment using [conda](https://conda.io):

```bash
conda env create -f environment.yml
```

For the exact, pinned environment use the lock file for your platform instead
(`conda create -n msa --file conda-lock/environment.linux-64.txt`).

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

* [input_mutation_tables/SIMrand](input_mutation_tables/SIMrand) folder contains simulated samples where each contains contributions from **5** randomly selected signatures out of **100** Poisson-generated signatures (the `sigRandom` tables in [signature_tables](signature_tables), generated with the legacy [legacy/generate_random_signatures.py](legacy/generate_random_signatures.py) script). To reproduce:
```bash
python legacy/generate_random_signatures.py -t SBS -c 96 -n 100
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

### Multiple mutation types
Mutation types run in parallel. Give them comma-separated on the command line, or as a list in
your config (`params.mutation_types = ['SBS', 'DBS', 'ID']`):
```bash
nextflow run run_auto_optimised_analysis.nf -profile conda \
    --dataset my_data \
    --mutation_types SBS,DBS,ID \
    --output_path results
```

## Changes from DSL1

Key improvements in the DSL2 version:
- Modern Nextflow syntax (DSL2)
- Improved caching and resume functionality
- Support for specific file inputs
- Cleaner module organization
- Intermediate files kept in the work directory (published to `temp/` with `--cleanup_temp false`)
- Better error handling and validation
- Unified output directory structure

## Citation

If you use this pipeline, please cite the original publication and the MSA repository.
