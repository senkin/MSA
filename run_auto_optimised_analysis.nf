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

// input data
params.signatures_file = null // optional direct path to signatures file, if null, default signature tables will be used from params.signature_tables
params.input_mutation_table = null // optional direct path to input mutation table file, if null, default input tables will be used from params.input_tables
params.SP_extractor_output_path = null // optional path to SigProfilerExtractor output, if null, default input will be used
params.SP_matrix_generator_output_path = null // optional path to SigProfilerMatrixGenerator output, if null, default input will be used
params.COSMIC_signatures = false // if set to true, COSMIC signatures are used form SigProfiler output, otherwise de-novo ones are used
params.dataset = 'SIM_test' // dataset name. Input matrices to be provided in params.input_tables/params.dataset, unless SigProfiler inputs are used
params.mutation_types = ['SBS'] // add or remove mutation types if needed
params.input_tables = "${workflow.projectDir}/input_mutation_tables"
params.SBS_context = 96 // 96, 192, 288, 1536 context matrices can be provided (SBS only)
params.number_of_samples = -1 // number of samples to analyse (-1 means all available)

// output paths
params.output_path = "${workflow.launchDir}" // specify full path to output folder, default is current working directory
params.tables_output_path = params.output_path + "/output_tables"
params.plots_output_path = params.output_path + "/plots"
params.temp_path = params.output_path + "/temp"

// signatures to use
params.signature_tables = "${workflow.projectDir}/signature_tables"
params.signature_prefix = "sigProfiler" // prefix of signature files to use (e.g. sigProfiler, sigRandom)

// simulations parameters
params.run_only_simulations = false // set to true if only simulations are needed, these will be produced in output_tables folder
params.number_of_simulated_samples = -1 // number of simulations to run (-1 means automatically apply a rounded factor of ten but not less than 1000)
params.add_noise = true // add noise in simulations (recommended)
params.noise_type = "gaussian" // set the type of noise in simulations: gaussian, poisson or negative_binomial (Gaussian by default)
params.noise_stdev = 10 // set standard deviation of gaussian noise, in percentage of sample mutation burden (10 percent by default)
params.zero_inflation_threshold = 0.05 // set the relative threshold below which all simulated signature activities are set to zero (0.01 by default)

// optimisation flag and parameters
params.run_only_optimisation = false // set to true if only optimisation is required, without final attributions
params.optimisation_NNLS_output_path = params.output_path + "/outputs_optimisation"
params.optimisation_plots_output_path = params.plots_output_path + "/optimisation_plots"
params.optimised = true // if set to false, optimisation will run but not be used in final attributions
params.optimisation_strategy = "removal" // optimisation strategy (removal, addition or add-remove)
params.weak_thresholds = [0, 0.0001, 0.0002, 0.0003] // range of L2 similarity decrease thresholds to be scanned, excluding weakest signatures - adjust if needed
params.strong_thresholds = [0] // range of L2 similarity increase thresholds to be scanned, including strongest signatures: only one is sufficient in default removal strategy
params.bootstrap_method = "binomial" // bootstrap flag and method (binomial, multinomial, residuals, classic, bootstrap_residuals)
params.number_of_bootstrapped_samples_in_optimisation = 100 // at least 100 is recommended
params.number_of_bootstrapped_samples = 1000 // bootstrap variations in final attribution, at least 1000 is recommended
params.confidence_level = 0.95 // specify the confidence level for CI calculation (default: 0.95)
params.use_absolute_attributions = false // use absolute mutation counts in final bootstrap outputs (relative by default)
params.metric_to_prioritise = "specificity" // set a metric to prioritise (default: specificity), requiring at least the specified threshold or closest alternative
params.metric_threshold = 0.95 // specify the minimum threshold of the prioritised metric
params.signatures_to_prioritise = [] // set a list of signatures to prioritise (empty list means all, by default)
params.no_CI_for_penalties = false // do not use confidence intervals for optimal penalties calculation
params.calculate_penalty_on_average = false // apply criteria based on signatures overall (on average, less conservative), rather than maximising prioritised metric for every signature (more conservative)

// plotting flags
params.plot_optimisation_plots = true
params.plot_bootstrap_attributions = true
params.plot_metrics = true
params.plot_signatures = true
params.plot_input_spectra = true
params.plot_fitted_spectra = false
params.plot_residuals = false
params.show_poisson_errors = false
params.show_strands = false // only works with higher contexts (192, 288)
params.show_nontranscribed_region = false // only wortks with higher contexts (288)

// if SIM in dataset name (synthetic data), use the following percentage range for measuring signature attirbution sensitivities
params.signature_attribution_thresholds = 0..20

params.help = null

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
include { FINAL_NNLS_workflow } from './modules/nnls'
include { FINAL_NNLS_BOOTSTRAP_workflow } from './modules/nnls'
include { simulate_data_workflow } from './modules/simulations'
include { BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables'
include { FINAL_BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables'
include { OPTIMAL_PENALTIES_workflow } from './modules/optimal_penalties'
include { OPTIMISATION_PLOTS_workflow } from './modules/optimisation_plots'

// Main workflow
workflow {
    // Set up input channels based on parameters
    // Priority: specific files > SigProfiler outputs > default directories
    
    // Handle specific file inputs (highest priority)
    if (params.signatures_file) {
        // Convert specific signature file
        convert_data_workflow(params.dataset, file(params.signatures_file).toAbsolutePath(), 'specific_signature_files', mutation_types)
        signature_files_channel = convert_data_workflow.out.signature_files
        
        if (params.plot_signatures) {
            convert_data_workflow.out.signature_files.collect().view { "Conversion complete, starting plots..." }
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/signature_tables", 'signatures')
            }
        }
    } else if (params.SP_extractor_output_path) {
        // Convert SigProfiler extractor output to temp location
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'signature_tables', mutation_types)
        signature_files_channel = convert_data_workflow.out.signature_files
        
        if (params.plot_signatures) {
            convert_data_workflow.out.signature_files.collect().view { "Conversion complete, starting plots..." }
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/signature_tables", 'signatures')
            }
        }
    } else {
        // Use default signature tables
        stage_default_signatures_workflow()
        signature_files_channel = stage_default_signatures_workflow.out.staged_signature_tables
        
        if (params.plot_signatures) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, params.signature_tables, 'signatures')
            }
        }
    }
    
    // Handle specific mutation table input (highest priority)
    if (params.input_mutation_table) {
        // Convert specific mutation table file
        convert_data_workflow(params.dataset, file(params.input_mutation_table).toAbsolutePath(), 'specific_mutation_files', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files
        
        if (params.plot_input_spectra) {
            convert_data_workflow.out.input_files.collect().view { "Conversion complete, starting plots..." }
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/input_tables", 'mutation_spectra')
            }
        }
    } else if (params.SP_extractor_output_path && !params.SP_matrix_generator_output_path) {
        // Convert extractor matrices to temp location
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'extractor_matrices', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files
        
        if (params.plot_input_spectra) {
            convert_data_workflow.out.input_files.collect().view { "Conversion complete, starting plots..." }
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/input_tables", 'mutation_spectra')
            }
        }
    } else if (params.SP_matrix_generator_output_path) {
        // Convert matrix generator output to temp location
        convert_data_workflow(params.dataset, params.SP_matrix_generator_output_path, 'matrix_generator_matrices', mutation_types)
        input_files_channel = convert_data_workflow.out.input_files
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/input_tables", 'mutation_spectra')
            }
        }
    } else {
        // Use default input tables - collect them into a channel
        stage_default_inputs_workflow(params.dataset)
        input_files_channel = stage_default_inputs_workflow.out.staged_input_tables
        
        if (params.plot_input_spectra) {
            for (mutation_type in mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, params.input_tables, 'mutation_spectra')
            }
        }
    }

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
            UnoptimizedNNLS.out.mutations_table
        )
        
        def optimized_outputs = []
        def bootstrap_outputs = []
        // Run optimized NNLS for each threshold combination
        for (weak_threshold in params.weak_thresholds) {
            for (strong_threshold in params.strong_thresholds) {
                // Standard optimized NNLS
                OptimizedNNLS(
                    params.dataset,
                    mutation_type,
                    simulate_data_workflow.out.all_simulation_outputs,
                    signature_files_channel,
                    weak_threshold,
                    strong_threshold
                )
                optimized_outputs.add(OptimizedNNLS.out.dataset_mutation_pairs)
                
                // Optimized NNLS with bootstrap
                OptimizedNNLSforBootstrap(
                    params.dataset,
                    mutation_type,
                    simulate_data_workflow.out.all_simulation_outputs,
                    signature_files_channel,
                    weak_threshold,
                    strong_threshold,
                    params.number_of_bootstrapped_samples_in_optimisation
                )
                bootstrap_outputs.add(OptimizedNNLSforBootstrap.out.bootstrap_indices)
            }
        }

        // Wait for ALL processes to complete
        all_optimized = Channel.empty().mix(*optimized_outputs).collect()
        all_bootstrap = Channel.empty().mix(*bootstrap_outputs).collect()
        
        // Generate bootstrap tables after ALL optimized runs complete
        // Wait for both OptimizedNNLS and OptimizedNNLSforBootstrap to finish
        BOOTSTRAP_TABLES_workflow(
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            all_optimized.mix(all_bootstrap).collect(),
            params.weak_thresholds,
            params.strong_thresholds,
            params.number_of_bootstrapped_samples_in_optimisation
        )
        
        // Calculate optimal penalties
        OPTIMAL_PENALTIES_workflow(
            BOOTSTRAP_TABLES_workflow.out.bootstrap_tables,
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            params.weak_thresholds,
            params.strong_thresholds
        )
        
        // Generate optimization plots
        OPTIMISATION_PLOTS_workflow(
            OPTIMAL_PENALTIES_workflow.out.penalties_for_optimisation_plotting,
            params.weak_thresholds,
            params.strong_thresholds
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
            params.number_of_bootstrapped_samples
        )

        // Generate final bootstrap tables after final bootstrap NNLS completes
        FINAL_BOOTSTRAP_TABLES_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_NNLS_BOOTSTRAP_workflow.out.bootstrap_indices,
            params.number_of_bootstrapped_samples
        )

        // Generate final plots
        ALL_FINAL_PLOTS_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_BOOTSTRAP_TABLES_workflow.out.attributions_per_sample,
            FINAL_BOOTSTRAP_TABLES_workflow.out.signature_prevalences
        )
    }
}

