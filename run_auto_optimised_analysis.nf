#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
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

log.info ''
log.info '--------------------------------------------------------'
log.info '              __   __  _____                             '
log.info '             |  \\/  |/ ____|  /\\                      '
log.info '             | \\  / | (___   /  \\                     '
log.info '             | |\\/| |\\___ \\ / /\\ \\                 '
log.info '             | |  | |____) / ____ \\                    '
log.info '             |_|  |_|_____/_/    \\_\\                  '
log.info '                                                        '
log.info '          MUTATIONAL SIGNATURE ANALYSIS v3.0            '
log.info '--------------------------------------------------------'
log.info 'Copyright (C) Sergey Senkin'
log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE'
log.info 'This is free software, and you are welcome to redistribute it'
log.info 'under certain conditions; see LICENSE for details.'
log.info '--------------------------------------------------------'
log.info ''

if (params.help) {
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
    log.info "     --input_tables <path>      : Directory with mutation tables"
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
    log.info "      All inputs are converted and stored in the temp directory for processing."
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

// add parameter values to log output (.nextflow.log)
log.info params.collect { k,v -> "${k.padRight(34)}: $v" }.join("\n")

// Normalize mutation_types to always be a list (handle command-line string input)
if (params.mutation_types instanceof String) {
    mutation_types = [params.mutation_types]
} else {
    mutation_types = params.mutation_types
}

// Validate the threshold scan settings
if (params.weak_thresholds.toString().trim() == 'auto' && params.optimisation_strategy != 'removal') {
    error "ERROR: weak_thresholds = 'auto' needs optimisation_strategy = 'removal' (got '${params.optimisation_strategy}'). Give an explicit list, e.g. --weak_thresholds 0,0.0001,0.0002,0.0003"
}
if (!(params.threshold_scan in ['trajectory', 'separate'])) {
    error "ERROR: threshold_scan must be 'trajectory' or 'separate' (got '${params.threshold_scan}')"
}
// The optimisation (unless only simulating) and the final attribution (unless only optimising)
// need at least one bootstrap iteration
def is_positive_integer(value) { value.toString().isInteger() && value.toString().toInteger() > 0 }
def bootstrap_counts_needed = [number_of_bootstrapped_samples_in_optimisation: !params.run_only_simulations,
                               number_of_bootstrapped_samples: !params.run_only_simulations && !params.run_only_optimisation]
bootstrap_counts_needed.each { name, needed ->
    if (needed && !is_positive_integer(params[name])) {
        error "ERROR: ${name} must be a positive integer (got '${params[name]}')"
    }
}

// Validate parameters when using specific file inputs
if (params.input_mutation_table || params.signatures_file) {
    if (params.input_mutation_table) {
        if (mutation_types.size() != 1) {
            error "ERROR: When using --input_mutation_table, please specify exactly ONE mutation type. Got: ${mutation_types}. Use: --mutation_types SBS (or DBS, ID, etc.)"
        }
        if (mutation_types[0] == 'SBS') {
            // For SBS, also check that only one context is specified
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


// Include modules
include { stage_default_signatures_workflow } from './modules/stage_default_inputs'
include { stage_default_inputs_workflow } from './modules/stage_default_inputs'
include { plot_spectra_workflow } from './modules/plotting'
include { ALL_FINAL_PLOTS_workflow } from './modules/plotting'
include { convert_data_workflow } from './modules/data_conversion'
include { NNLS_unoptimized_workflow as UnoptimizedNNLS } from './modules/nnls'
include { NNLS_optimized_workflow as OptimizedNNLS } from './modules/nnls'
include { NNLS_bootstrap_workflow as OptimizedNNLSforBootstrap } from './modules/nnls'
include { WEAK_THRESHOLDS_workflow } from './modules/nnls'
include { FINAL_NNLS_workflow } from './modules/nnls'
include { FINAL_NNLS_BOOTSTRAP_workflow } from './modules/nnls'
include { simulate_data_workflow } from './modules/simulations'
include { BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables'
include { FINAL_BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables'
include { OPTIMAL_PENALTIES_workflow } from './modules/optimal_penalties'
include { OPTIMISATION_PLOTS_workflow } from './modules/optimisation_plots'
include { cleanup_workflow } from './modules/cleanup'

// Write what is needed to reproduce this run (resolved seed, code revision, versions,
// profile, all parameters) next to its results, in <output_path>/pipeline_info/.
def write_run_info(run_seed, seed_source) {
    def git_commit = 'n/a'
    try {
        def proc = ['git', '-C', "${workflow.projectDir}", 'rev-parse', 'HEAD'].execute()
        proc.waitFor()
        if (proc.exitValue() == 0) git_commit = proc.text.trim()
    } catch (Exception ignored) { }
    def info = [
        'run name'         : workflow.runName,
        'session id'       : workflow.sessionId,
        'start'            : workflow.start,
        'command line'     : workflow.commandLine,
        'project dir'      : workflow.projectDir,
        'revision'         : workflow.revision ?: 'n/a',
        'commit id'        : workflow.commitId ?: git_commit,
        'nextflow version' : workflow.nextflow.version,
        'profile'          : workflow.profile,
        'container engine' : workflow.containerEngine ?: 'none',
        'container'        : workflow.container ?: 'none',
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

// Main workflow
workflow {
    // Resolve the run-wide RNG seed once; it seeds both the simulations and the bootstrap
    // resampling. Without --seed, a seed is derived from the Nextflow session ID: new for
    // every fresh run, unchanged on -resume (so cached tasks stay valid), and recorded in
    // pipeline_info/ so any run can be repeated exactly with --seed <value>.
    def run_seed = params.seed != null ? (params.seed as long)
        : ((workflow.sessionId.mostSignificantBits ^ workflow.sessionId.leastSignificantBits) & 0x7fffffffL)
    def seed_source = params.seed != null ? '--seed' : 'derived from session ID'
    log.info "RNG seed: ${run_seed} (${seed_source}; rerun with --seed ${run_seed} to reproduce)"
    write_run_info(run_seed, seed_source)

    // Set up input channels based on parameters
    // Priority: specific files > SigProfiler outputs > default directories
    
    // Handle specific file inputs (highest priority)
    if (params.signatures_file) {
        // Convert specific signature file
        convert_data_workflow(params.dataset, file(params.signatures_file).toAbsolutePath(), 'specific_signature_files', mutation_types)
        signature_files_channel = convert_data_workflow.out.signature_files.collect()
        
        if (params.plot_signatures) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.output_path}/temp/signature_tables", 'signatures', signature_files_channel)
            }
        }
    } else if (params.SP_extractor_output_path) {
        // Convert SigProfiler extractor output to temp location
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'signature_tables', mutation_types)
        signature_files_channel = convert_data_workflow.out.signature_files.collect()
        
        if (params.plot_signatures) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.output_path}/temp/signature_tables", 'signatures', signature_files_channel)
                }
        }
    } else {
        // Use default signature tables
        stage_default_signatures_workflow()
        signature_files_channel = stage_default_signatures_workflow.out.staged_signature_tables.collect()
        
        if (params.plot_signatures) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, params.signature_tables, 'signatures', signature_files_channel)
            }
        }
    }
    
    // Handle specific mutation table input (highest priority)
    if (params.input_mutation_table) {
        // Convert specific mutation table file
        convert_data_workflow(params.dataset, file(params.input_mutation_table).toAbsolutePath(), 'specific_mutation_files', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files.collect()
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.output_path}/temp/input_tables", 'mutation_spectra', input_files_channel)
            }
        }
    } else if (params.SP_extractor_output_path && !params.SP_matrix_generator_output_path) {
        // Convert extractor matrices to temp location
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'extractor_matrices', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files.collect()
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.output_path}/temp/input_tables", 'mutation_spectra', input_files_channel)
            }
        }
    } else if (params.SP_matrix_generator_output_path) {
        // Convert matrix generator output to temp location
        convert_data_workflow(params.dataset, params.SP_matrix_generator_output_path, 'matrix_generator_matrices', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files.collect()
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.output_path}/temp/input_tables", 'mutation_spectra', input_files_channel)
            }
        }
    } else {
        // Use default input tables - collect them into a channel
        stage_default_inputs_workflow(params.dataset)
        input_files_channel = stage_default_inputs_workflow.out.staged_input_tables.collect()
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, params.input_tables, 'mutation_spectra', input_files_channel)
            }
        }
    }
    def all_plot_outputs = []

    // Process each mutation type sequentially to avoid channel conflicts
    for (mutation_type in mutation_types) {
        // Run unoptimized NNLS
        UnoptimizedNNLS(
            params.dataset,
            mutation_type,
            input_files_channel,
            signature_files_channel
        )

        // Run simulation workflow
        simulate_data_workflow(
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            UnoptimizedNNLS.out.mutations_table,
            run_seed
        )
        
        // Thresholds to scan: 'auto' (from removal trajectories on the simulated data) or
        // the configured list, resolved into the labels every later step uses
        WEAK_THRESHOLDS_workflow(
            params.dataset,
            mutation_type,
            simulate_data_workflow.out.all_simulation_outputs,
            signature_files_channel
        )

        // Optimized NNLS on the simulated data: one job for all threshold pairs
        OptimizedNNLS(
            params.dataset,
            mutation_type,
            simulate_data_workflow.out.all_simulation_outputs,
            signature_files_channel,
            WEAK_THRESHOLDS_workflow.out.weak_file,
            WEAK_THRESHOLDS_workflow.out.strong_file
        )

        // Optimized NNLS with bootstrap: each chunk of iterations covers all threshold pairs
        OptimizedNNLSforBootstrap(
            params.dataset,
            mutation_type,
            simulate_data_workflow.out.all_simulation_outputs,
            signature_files_channel,
            WEAK_THRESHOLDS_workflow.out.weak_file,
            WEAK_THRESHOLDS_workflow.out.strong_file,
            params.number_of_bootstrapped_samples_in_optimisation,
            run_seed
        )

        // Generate bootstrap tables after ALL optimized runs complete
        // Wait for both OptimizedNNLS and OptimizedNNLSforBootstrap to finish
        BOOTSTRAP_TABLES_workflow(
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            OptimizedNNLS.out.dataset_mutation_pairs.collect().mix(OptimizedNNLSforBootstrap.out.bootstrap_done.collect()).collect(),
            WEAK_THRESHOLDS_workflow.out.weak_labels,
            WEAK_THRESHOLDS_workflow.out.strong_labels,
            params.number_of_bootstrapped_samples_in_optimisation
        )
        
        // Calculate optimal penalties
        OPTIMAL_PENALTIES_workflow(
            BOOTSTRAP_TABLES_workflow.out.bootstrap_tables,
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            WEAK_THRESHOLDS_workflow.out.weak_labels,
            WEAK_THRESHOLDS_workflow.out.strong_labels
        )
        
        // Generate optimization plots
        OPTIMISATION_PLOTS_workflow(
            OPTIMAL_PENALTIES_workflow.out.penalties_for_optimisation_plotting,
            WEAK_THRESHOLDS_workflow.out.weak_labels,
            WEAK_THRESHOLDS_workflow.out.strong_labels
        )
        
        // Run final NNLS with optimal penalties
        FINAL_NNLS_workflow(
            OPTIMAL_PENALTIES_workflow.out.penalties_for_central_NNLS_attribution,
            OPTIMAL_PENALTIES_workflow.out.optimal_weak_penalty,
            OPTIMAL_PENALTIES_workflow.out.optimal_strong_penalty,
            input_files_channel,
            signature_files_channel
        )
        
        // Run final bootstrap NNLS with optimal penalties
        FINAL_NNLS_BOOTSTRAP_workflow(
            OPTIMAL_PENALTIES_workflow.out.penalties_for_bootstrap_NNLS_attribution,
            OPTIMAL_PENALTIES_workflow.out.optimal_weak_penalty,
            OPTIMAL_PENALTIES_workflow.out.optimal_strong_penalty,
            input_files_channel,
            signature_files_channel,
            params.number_of_bootstrapped_samples,
            run_seed
        )

        // Generate final bootstrap tables after final bootstrap NNLS completes
        FINAL_BOOTSTRAP_TABLES_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_NNLS_BOOTSTRAP_workflow.out.bootstrap_done,
            params.number_of_bootstrapped_samples
        )

        // Generate final plots
        ALL_FINAL_PLOTS_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_BOOTSTRAP_TABLES_workflow.out.attributions_per_sample,
            FINAL_BOOTSTRAP_TABLES_workflow.out.signature_prevalences
        )
        all_plot_outputs.add(ALL_FINAL_PLOTS_workflow.out.bootstrap_plots)
        all_plot_outputs.add(ALL_FINAL_PLOTS_workflow.out.metrics_plots)
        all_plot_outputs.add(ALL_FINAL_PLOTS_workflow.out.fitted_plots)
        all_plot_outputs.add(ALL_FINAL_PLOTS_workflow.out.residuals_plots)
    }
    // Cleanup temporary files if specified
    if (params.cleanup_temp) {
        // Combine all accumulated channels
        all_plots_done = Channel.empty().mix(*all_plot_outputs).collect()
        cleanup_workflow(all_plots_done)
    }
}

