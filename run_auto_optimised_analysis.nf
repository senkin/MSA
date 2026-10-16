#!/usr/bin/env nextflow
// Copyright (C) 2025 Sergey Senkin

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Data flow: every process receives the files it reads through channels, staged in its work
// directory, and writes only there. Results are published to <output_path> for the user and
// never read back by the pipeline; intermediate files are published to <output_path> only with
// cleanup_temp = false. Mutation types run in parallel, as (dataset, mutation_type) keys.

include { convert_data as convert_signatures; convert_data as convert_inputs } from './modules/data_conversion'
include { run_unoptimized_NNLS; resolve_thresholds; run_optimized_NNLS; run_bootstrap_NNLS; run_final_NNLS; run_final_bootstrap_NNLS } from './modules/nnls'
include { simulate_data } from './modules/simulations'
include { make_optimisation_bootstrap_tables; make_bootstrap_tables } from './modules/bootstrap_tables'
include { calculate_optimal_penalties } from './modules/optimal_penalties'
include { plot_optimisation_plots } from './modules/optimisation_plots'
include { plot_signatures; plot_input_spectra; plot_bootstrap_attributions; plot_metrics; plot_fitted_spectra; plot_residuals } from './modules/plotting'
include { bootstrap_chunks } from './modules/common'

def print_banner() {
    log.info ''
    log.info '--------------------------------------------------------'
    log.info '              __   __  _____                             '
    log.info '             |  \\/  |/ ____|  /\\                      '
    log.info '             | \\  / | (___   /  \\                     '
    log.info '             | |\\/| |\\___ \\ / /\\ \\                 '
    log.info '             | |  | |____) / ____ \\                    '
    log.info '             |_|  |_|_____/_/    \\_\\                  '
    log.info '                                                        '
    log.info "          MUTATIONAL SIGNATURE ANALYSIS v${workflow.manifest.version}"
    log.info '--------------------------------------------------------'
    log.info 'Copyright (C) Sergey Senkin'
    log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE'
    log.info 'This is free software, and you are welcome to redistribute it'
    log.info 'under certain conditions; see LICENSE for details.'
    log.info '--------------------------------------------------------'
    log.info ''
}

def print_help() {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run run_auto_optimised_analysis.nf [options]"
    log.info ""
    log.info "INPUT OPTIONS (in order of priority):"
    log.info "  1. Specific files (highest priority):"
    log.info "     --input_mutation_table <path>  : Specific mutation table file"
    log.info "     --signatures_file <path>       : Specific signature file"
    log.info "     NOTE: When using specific files, you MUST specify:"
    log.info "           --mutation_types <type>  : Single mutation type (e.g., SBS, DBS, ID)"
    log.info "           --SBS_context <context>  : Single context for SBS (e.g., 96, 192, 288)"
    log.info ""
    log.info "  2. SigProfiler outputs:"
    log.info "     --SP_extractor_output_path <path>         : SigProfiler extractor output"
    log.info "     --SP_matrix_generator_output_path <path>  : SigProfiler matrix generator output"
    log.info ""
    log.info "  3. Default directories (lowest priority):"
    log.info "     --input_tables <path>      : Directory with mutation tables (in <path>/<dataset>/)"
    log.info "     --signature_tables <path>  : Directory with signature tables"
    log.info ""
    log.info "EXAMPLE USAGE:"
    log.info "  # Using specific files"
    log.info "  nextflow run run_auto_optimised_analysis.nf \\"
    log.info "    --input_mutation_table my_sbs.txt \\"
    log.info "    --signatures_file my_sigs.txt \\"
    log.info "    --mutation_types SBS \\"
    log.info "    --SBS_context 96 \\"
    log.info "    --dataset my_data"
    log.info ""
    log.info "  # Using default directories"
    log.info "  nextflow run run_auto_optimised_analysis.nf --dataset my_data"
    log.info ""
    log.info "Note: Specific file inputs override SigProfiler outputs, which override default directories."
    log.info "      Converted inputs are kept in <output_path>/temp with --cleanup_temp false."
    log.info ""
    log.info "All parameters: https://gitlab.com/s.senkin/MSA/-/wikis/Parameters-description-table"
    log.info ""
}

// Mutation types as a list: from the config, or from the command line (e.g. --mutation_types SBS,DBS)
def mutation_types_list() {
    return params.mutation_types instanceof String ? params.mutation_types.tokenize(',')*.trim() : params.mutation_types
}

def is_positive_integer(value) {
    return value.toString().isInteger() && value.toString().toInteger() > 0
}

// Boolean params are read as params.<name>.toBoolean() throughout, which takes booleans (config,
// params file, legacy parser) and strings: with Nextflow's strict parser (the default from 26.04),
// values given on the command line stay strings, and 'false' would be true in Groovy
def validate_parameters(mutation_types) {
    ['help', 'COSMIC_signatures', 'cleanup_temp', 'run_only_simulations', 'add_noise', 'run_only_optimisation', 'optimised',
     'use_GPU', 'no_CI_for_penalties', 'calculate_penalty_on_average', 'plot_optimisation_plots',
     'plot_bootstrap_attributions', 'plot_metrics', 'plot_signatures', 'plot_input_spectra', 'plot_fitted_spectra',
     'plot_residuals', 'show_poisson_errors', 'show_strands', 'show_nontranscribed_region'].each { name ->
        if (!(params[name].toString().trim().toLowerCase() in ['true', 'false'])) {
            error "ERROR: ${name} must be true or false (got '${params[name]}')"
        }
    }
    // Options the scripts check only when they run, possibly hours into a run (and noise_type not at
    // all: any other value silently adds no noise)
    def choices = [
        optimisation_strategy: ['removal', 'addition', 'add-remove'],
        bootstrap_method     : ['binomial', 'multinomial', 'multinomial_weight', 'poisson', 'classic', 'bootstrap_residuals'],
        metric_to_prioritise : ['specificity', 'sensitivity'],
        noise_type           : ['gaussian', 'Gaussian', 'normal', 'Normal', 'poisson', 'Poisson', 'negative_binomial', 'Negative_binomial'],
        gpu_precision        : ['mixed', 'fp32', 'fp64'],
    ]
    choices.each { name, allowed ->
        if (!(params[name].toString() in allowed)) {
            error "ERROR: ${name} must be one of ${allowed.join(', ')} (got '${params[name]}')"
        }
    }
    // Penalties: numbers ('auto' for the weak ones), from a config list or comma-separated
    ['weak_thresholds', 'strong_thresholds'].each { name ->
        def value = params[name]
        if (!(name == 'weak_thresholds' && value.toString().trim() == 'auto')) {
            def items = value instanceof Collection ? value : value.toString().tokenize(', ')
            if (!items || items.any { item -> !item.toString().trim().isNumber() }) {
                error "ERROR: ${name} must be ${name == 'weak_thresholds' ? "'auto' or " : ''}a list of numbers (got '${value}')"
            }
        }
    }
    // Threshold scan settings
    if (params.weak_thresholds.toString().trim() == 'auto' && params.optimisation_strategy != 'removal') {
        error "ERROR: weak_thresholds = 'auto' needs optimisation_strategy = 'removal' (got '${params.optimisation_strategy}'). Give an explicit list, e.g. --weak_thresholds 0,0.0001,0.0002,0.0003"
    }
    if (!(params.threshold_scan in ['trajectory', 'separate'])) {
        error "ERROR: threshold_scan must be 'trajectory' or 'separate' (got '${params.threshold_scan}')"
    }
    // The optimisation (unless only simulating) and the final attribution (unless only optimising)
    // need at least one bootstrap iteration
    def only_simulations = params.run_only_simulations.toBoolean()
    def only_optimisation = params.run_only_optimisation.toBoolean()
    def bootstrap_counts_needed = [number_of_bootstrapped_samples_in_optimisation: !only_simulations,
                                   number_of_bootstrapped_samples: !only_simulations && !only_optimisation]
    bootstrap_counts_needed.each { name, needed ->
        if (needed && !is_positive_integer(params[name])) {
            error "ERROR: ${name} must be a positive integer (got '${params[name]}')"
        }
    }
    // Specific file inputs
    if (params.input_mutation_table) {
        if (mutation_types.size() != 1) {
            error "ERROR: When using --input_mutation_table, please specify exactly ONE mutation type. Got: ${mutation_types}. Use: --mutation_types SBS (or DBS, ID, etc.)"
        }
        if (mutation_types[0] == 'SBS') {
            log.info "Using specific mutation table for ${mutation_types[0]} with context ${params.SBS_context}"
        } else {
            log.info "Using specific mutation table for ${mutation_types[0]}"
        }
    }
    if (params.signatures_file) {
        if (mutation_types.size() != 1) {
            error "ERROR: When using --signatures_file, please specify exactly ONE mutation type. Got: ${mutation_types}. Use: --mutation_types SBS (or DBS, ID, etc.)"
        }
        if (mutation_types[0] == 'SBS') {
            log.info "Using specific signature file for ${mutation_types[0]} with context ${params.SBS_context}"
        } else {
            log.info "Using specific signature file for ${mutation_types[0]}"
        }
    }
}

// Which backend runs the NNLS, and how many bootstrap iterations each task runs. Without
// use_GPU, the NNLS jobs of a GPU config still go to its GPU queue but run on the CPU.
def log_nnls_backend() {
    def use_gpu = params.use_GPU.toBoolean()
    def chunk_size = (params.bootstrap_chunk_size ?: 0) as int
    def per_task = chunk_size > 0 ? "up to ${chunk_size} bootstrap iterations per task"
        : (use_gpu ? 'all bootstrap iterations of an attribution in one task' : 'one bootstrap iteration per task')
    if (use_gpu) {
        log.info "NNLS backend: GPU (cuML, ${params.gpu_precision} precision), ${per_task}, env: ${params.gpu_conda_env ?: 'the default one'}"
        if (!params.gpu_conda_env && workflow.profile.tokenize(',').contains('conda')) {
            log.warn "use_GPU = true without gpu_conda_env: the GPU NNLS processes run in the default conda env, which has no cuML. Set gpu_conda_env to the full path of the env built by setup_gpu_env.sh"
        }
    } else {
        log.info "NNLS backend: CPU, ${per_task}"
    }
}

// Write what is needed to reproduce this run (resolved seed, code revision, versions,
// profile, all parameters) next to its results, in <output_path>/pipeline_info/.
def write_run_info(run_seed, seed_source) {
    def git_commit = 'n/a'
    try {
        def proc = ['git', '-C', "${workflow.projectDir}", 'rev-parse', 'HEAD'].execute()
        proc.waitFor()
        if (proc.exitValue() == 0) {
            git_commit = proc.text.trim()
        }
    } catch (Exception _ignored) { }
    def info = [
        'run name'         : workflow.runName,
        'session id'       : workflow.sessionId,
        'start'            : workflow.start,
        'command line'     : workflow.commandLine,
        'project dir'      : workflow.projectDir,
        'revision'         : workflow.revision ?: 'n/a',
        'commit id'        : workflow.commitId ?: git_commit,
        'MSA version'      : workflow.manifest.version,
        'nextflow version' : workflow.nextflow.version,
        'profile'          : workflow.profile,
        'container engine' : workflow.containerEngine ?: 'none',
        // the container directive is a closure (see nextflow.config), so workflow.container is unset
        'container'        : workflow.containerEngine ? params.container_image : 'none',
        'rng seed'         : run_seed,
        'rng seed source'  : seed_source,
    ]
    def lines = info.collect { k, v -> "${k.padRight(18)}: ${v}" }
    lines << '' << 'params:'
    lines.addAll(params.sort().collect { k, v -> "  ${k.padRight(46)}: ${v}" })
    def info_dir = file("${params.output_path}/pipeline_info")
    info_dir.mkdirs()
    file("${info_dir}/run_info_${workflow.runName}.txt").text = lines.join('\n') + '\n'
}

// <weak>_<strong> labels of an optimisation output folder SIM_<dataset>_<context>_NNLS_<weak>_<strong>
def pair_suffix(folder, dataset) {
    def prefix = "SIM_${dataset}_${params.SBS_context}_NNLS_".toString()
    if (!folder.startsWith(prefix)) {
        error "ERROR: unexpected optimisation output folder ${folder} (expected ${prefix}<weak>_<strong>)"
    }
    return folder.substring(prefix.length())
}

// Choice of the penalties on the simulated data: threshold scan, bootstrap of every threshold
// pair, truth studies and the optimal penalties of each mutation type
workflow OPTIMISATION {
    take:
    simulated           // (dataset, mutation_type, simulated tables and truth)
    signature_files
    seed

    main:
    resolve_thresholds(simulated, signature_files)
    def scan = simulated.join(resolve_thresholds.out, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
    run_optimized_NNLS(scan, signature_files)
    run_bootstrap_NNLS(scan, signature_files,
        bootstrap_chunks(params.number_of_bootstrapped_samples_in_optimisation, params.bootstrap_chunk_size), seed)

    // The central attribution and the bootstrap iterations of each threshold pair, identified by
    // its output folder SIM_<dataset>_<context>_NNLS_<weak>_<strong>
    def central = run_optimized_NNLS.out.tables
        .transpose()
        .map { dataset, mutation_type, f -> tuple(dataset, mutation_type, pair_suffix(f.parent.name, dataset), f) }
        .groupTuple(by: [0, 1, 2])
    def bootstrap = run_bootstrap_NNLS.out.tables
        .transpose()
        .map { dataset, mutation_type, f -> tuple(dataset, mutation_type, pair_suffix(f.parent.parent.name, dataset), f) }
        .groupTuple(by: [0, 1, 2])
    make_optimisation_bootstrap_tables(central.join(bootstrap, by: [0, 1, 2], failOnMismatch: true, failOnDuplicate: true),
        signature_files, params.number_of_bootstrapped_samples_in_optimisation)

    // Every threshold pair of a mutation type, for the penalty choice and its plots
    def pairs = make_optimisation_bootstrap_tables.out.bundle.groupTuple(by: [0, 1])
    calculate_optimal_penalties(scan.join(pairs, by: [0, 1], failOnMismatch: true, failOnDuplicate: true), signature_files)

    if (params.plot_optimisation_plots.toBoolean()) {
        def plot_inputs = scan
            .map { dataset, mutation_type, _simulated_tables, weak, strong -> tuple(dataset, mutation_type, weak, strong) }
            .join(calculate_optimal_penalties.out.analysis, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
            .join(pairs, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
        plot_optimisation_plots(plot_inputs, signature_files)
    }

    emit:
    calculate_optimal_penalties.out.penalties   // (dataset, mutation_type, optimal weak penalty, optimal strong penalty)
}

// Final attribution of the input data with the optimal penalties: central attribution,
// bootstrap, summary tables and plots
workflow FINAL_ATTRIBUTION {
    take:
    penalties           // (dataset, mutation_type, optimal weak penalty, optimal strong penalty)
    input_files
    signature_files
    seed

    main:
    run_final_NNLS(penalties, input_files, signature_files)
    run_final_bootstrap_NNLS(penalties, input_files, signature_files,
        bootstrap_chunks(params.number_of_bootstrapped_samples, params.bootstrap_chunk_size), seed)
    def bootstrap = run_final_bootstrap_NNLS.out.tables.transpose().groupTuple(by: [0, 1])
    make_bootstrap_tables(run_final_NNLS.out.tables.join(bootstrap, by: [0, 1], failOnMismatch: true, failOnDuplicate: true),
        input_files, signature_files, params.number_of_bootstrapped_samples)

    // Attribution tables, bootstrap summaries and truth table (simulated data); truth studies
    def results = run_final_NNLS.out.tables
        .join(make_bootstrap_tables.out.tables, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
        .map { dataset, mutation_type, central, cis, prevalences, per_sample, per_signature, stat_metrics, pruned, truth, truth_studies ->
            tuple(dataset, mutation_type, [central, cis, prevalences, per_sample, per_signature, stat_metrics, pruned, truth].flatten(),
                [truth_studies].flatten())
        }

    if (params.plot_bootstrap_attributions.toBoolean()) {
        plot_bootstrap_attributions(results.map { dataset, mutation_type, tables, _truth_studies -> tuple(dataset, mutation_type, tables) },
            input_files, signature_files)
    }
    if (params.plot_metrics.toBoolean()) {
        plot_metrics(results)
    }
    if (params.plot_fitted_spectra.toBoolean()) {
        plot_fitted_spectra(run_final_NNLS.out.tables)
    }
    if (params.plot_residuals.toBoolean()) {
        plot_residuals(run_final_NNLS.out.tables)
    }
}

workflow {
    print_banner()
    if (params.help.toBoolean()) {
        print_help()
        return
    }
    log.info "help:                               ${params.help}"
    // add parameter values to log output (.nextflow.log)
    log.info params.collect { k, v -> "${k.padRight(34)}: $v" }.join("\n")

    def mutation_types = mutation_types_list()
    validate_parameters(mutation_types)
    log_nnls_backend()

    // Resolve the run-wide RNG seed once; it seeds both the simulations and the bootstrap
    // resampling. Without --seed, a seed is derived from the Nextflow session ID: new for
    // every fresh run, unchanged on -resume (so cached tasks stay valid), and recorded in
    // pipeline_info/ so any run can be repeated exactly with --seed <value>.
    def run_seed = params.seed != null ? (params.seed as long)
        : ((workflow.sessionId.mostSignificantBits ^ workflow.sessionId.leastSignificantBits) & 0x7fffffffL)
    def seed_source = params.seed != null ? '--seed' : 'derived from session ID'
    log.info "RNG seed: ${run_seed} (${seed_source}; rerun with --seed ${run_seed} to reproduce)"
    write_run_info(run_seed, seed_source)

    def dataset = params.dataset

    // Signature tables. Priority: specific file > SigProfiler extractor output > signature folder
    def reference_signatures = channel.fromPath("${params.signature_tables}/*.csv", checkIfExists: true).collect()
    def signature_files
    if (params.signatures_file) {
        convert_signatures(dataset, file(params.signatures_file, checkIfExists: true), 'specific_signature_files', mutation_types, reference_signatures)
        signature_files = convert_signatures.out.signature_files
            .ifEmpty { error "ERROR: no signature tables converted from ${params.signatures_file}" }.collect()
    } else if (params.SP_extractor_output_path) {
        convert_signatures(dataset, file(params.SP_extractor_output_path, checkIfExists: true), 'signature_tables', mutation_types, reference_signatures)
        signature_files = convert_signatures.out.signature_files
            .ifEmpty { error "ERROR: no signature tables converted from ${params.SP_extractor_output_path}" }.collect()
    } else {
        signature_files = reference_signatures
    }

    // Input mutation tables. Priority: specific file > SigProfiler outputs > <input_tables>/<dataset>/
    def input_files
    if (params.input_mutation_table) {
        convert_inputs(dataset, file(params.input_mutation_table, checkIfExists: true), 'specific_mutation_files', mutation_types, reference_signatures)
        input_files = convert_inputs.out.input_files
            .ifEmpty { error "ERROR: no mutation tables converted from ${params.input_mutation_table}" }.collect()
    } else if (params.SP_extractor_output_path && !params.SP_matrix_generator_output_path) {
        convert_inputs(dataset, file(params.SP_extractor_output_path, checkIfExists: true), 'extractor_matrices', mutation_types, reference_signatures)
        input_files = convert_inputs.out.input_files
            .ifEmpty { error "ERROR: no mutation tables converted from ${params.SP_extractor_output_path}" }.collect()
    } else if (params.SP_matrix_generator_output_path) {
        convert_inputs(dataset, file(params.SP_matrix_generator_output_path, checkIfExists: true), 'matrix_generator_matrices', mutation_types, reference_signatures)
        input_files = convert_inputs.out.input_files
            .ifEmpty { error "ERROR: no mutation tables converted from ${params.SP_matrix_generator_output_path}" }.collect()
    } else {
        input_files = channel.fromPath("${params.input_tables}/${dataset}/*.csv", checkIfExists: true).collect()
    }

    // One (dataset, mutation type) key per mutation type: they run in parallel
    def keys = channel.fromList(mutation_types).map { mutation_type -> tuple(dataset, mutation_type) }

    if (params.plot_signatures.toBoolean()) {
        plot_signatures(keys, signature_files)
    }
    if (params.plot_input_spectra.toBoolean()) {
        plot_input_spectra(keys, input_files)
    }

    // Simulated data from the unoptimised attribution of the input data
    run_unoptimized_NNLS(keys, input_files, signature_files)
    simulate_data(run_unoptimized_NNLS.out.mutations_table, signature_files, run_seed)

    if (!params.run_only_simulations.toBoolean()) {
        OPTIMISATION(simulate_data.out.tables, signature_files, run_seed)
        if (!params.run_only_optimisation.toBoolean()) {
            FINAL_ATTRIBUTION(OPTIMISATION.out, input_files, signature_files, run_seed)
        }
    }
}
