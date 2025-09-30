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
params.input_tables = "$baseDir/input_mutation_tables"
params.SBS_context = 96 // 96, 192, 288, 1536 context matrices can be provided (SBS only)
params.number_of_samples = -1 // number of samples to analyse (-1 means all available)

// output paths
params.output_path = "."
params.tables_output_path = params.output_path + "/output_tables"
params.plots_output_path = params.output_path + "/plots"
params.temp_path = params.output_path + "/temp"

// signatures to use
params.signature_tables = "$baseDir/signature_tables"
params.signature_prefix = "sigProfiler" // prefix of signature files to use (e.g. sigProfiler, sigRandom)

// simulations parameters
params.run_only_simulations = false // set to true if only simulations are needed, these will be produced in $baseDir/output_tables folder
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
params.weak_thresholds = [0, 0.0001]//, '0.0002', '0.0003'] // range of L2 similarity decrease thresholds to be scanned, excluding weakest signatures - adjust if needed
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

// helper flags for scripts (automatic based on parameters)
optimised_flag = (params.optimised) ? "-x" : ''
abs_flag = (params.use_absolute_attributions) ? "-a" : ''
suffix = (params.use_absolute_attributions) ? "abs_mutations" : 'weights'
error_flag = (params.show_poisson_errors) ? "-e" : ''
strands_flag = (params.show_strands) ? "-b" : ''
nontranscribed_flag = (params.show_nontranscribed_region) ? "-n" : ''
COSMIC_flag = (params.COSMIC_signatures) ? "-C" : ''
signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix
noise_flag = (params.add_noise) ? "-z" : ''
no_CI_for_penalties_flag = (params.no_CI_for_penalties) ? "--no_CI" : ''
calculate_penalty_on_average_flag = (params.calculate_penalty_on_average) ? "--average" : ''
// override number of samples/variations for test run
test_run = (params.dataset == 'SIM_test') ? true : false
number_of_bootstrapped_samples_in_optimisation = (test_run) ? 10 : params.number_of_bootstrapped_samples_in_optimisation
number_of_bootstrapped_samples = (test_run) ? 10 : params.number_of_bootstrapped_samples
number_of_simulated_samples = (test_run) ? 10 : params.number_of_simulated_samples
weak_thresholds = (test_run) ? [0, 0.01, 0.02] : params.weak_thresholds
strong_thresholds = (test_run) ? [0] : params.strong_thresholds
// mutation_types = (test_run) ? ['SBS', 'DBS', 'ID'] : params.mutation_types
mutation_types = (test_run) ? ['DBS'] : params.mutation_types

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
    log.info "nextflow run run_analysis.nf"
    log.info ""
    log.info "Nextflow currently does not support list parameters,"
    log.info "so please specify the parameters directly in the script."
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

// add parameter values to log output (.nextflow.log)
log.info params.collect { k,v -> "${k.padRight(34)}: $v" }.join("\n")


// Include modules
include { plot_spectra_workflow } from './modules/plotting' addParams(
    strands_flag: strands_flag,
    nontranscribed_flag: nontranscribed_flag,
    signature_prefix: signature_prefix
)

include { BOOTSTRAP_ATTRIBUTIONS_PLOTS_workflow } from './modules/plotting'
include { METRICS_PLOTS_workflow } from './modules/plotting'
include { FITTED_SPECTRA_PLOTS_workflow } from './modules/plotting'
include { RESIDUALS_PLOTS_workflow } from './modules/plotting'
include { ALL_FINAL_PLOTS_workflow } from './modules/plotting'

include { convert_data_workflow } from './modules/data_conversion' addParams(
    strands_flag: strands_flag,
    nontranscribed_flag: nontranscribed_flag,
    error_flag: error_flag,
    COSMIC_flag: COSMIC_flag,
    signature_prefix: signature_prefix
)


include { convert_specific_files_workflow } from './modules/data_conversion' addParams(
    strands_flag: strands_flag,
    nontranscribed_flag: nontranscribed_flag,
    error_flag: error_flag,
    COSMIC_flag: COSMIC_flag,
    signature_prefix: signature_prefix
)

include { NNLS_unoptimized_workflow as UnoptimizedNNLS } from './modules/nnls' addParams(
    signature_prefix: signature_prefix
)
include { NNLS_optimized_workflow as OptimizedNNLS } from './modules/nnls' addParams(
    signature_prefix: signature_prefix
)
include { NNLS_bootstrap_workflow as OptimizedNNLSforBootstrap } from './modules/nnls' addParams(
    signature_prefix: signature_prefix
)
include { FINAL_NNLS_workflow } from './modules/nnls' addParams(
    signature_prefix: signature_prefix
)
include { FINAL_NNLS_BOOTSTRAP_workflow } from './modules/nnls' addParams(
    signature_prefix: signature_prefix
)
include { simulate_data_workflow } from './modules/simulations' addParams(
    signature_prefix: signature_prefix
)
include { BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables' addParams(
    signature_prefix: signature_prefix
)
include { FINAL_BOOTSTRAP_TABLES_workflow } from './modules/bootstrap_tables' addParams(
    signature_prefix: signature_prefix
)
include { OPTIMAL_PENALTIES_workflow } from './modules/optimal_penalties' addParams(
    signature_prefix: signature_prefix
)
include { OPTIMISATION_PLOTS_workflow } from './modules/optimisation_plots' addParams(
    signature_prefix: signature_prefix
)

// Main workflow
workflow {

    if (params.SP_extractor_output_path) {
        // Convert SigProfiler extractor output to temp location
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'signature_tables')
        signature_files_channel = convert_data_workflow.out.signature_files
        
        // Plot signatures if requested
        if (params.plot_signatures) {
            for (mutation_type in params.mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/signature_tables", 'signatures')
            }
        }
        
        if (!params.SP_matrix_generator_output_path) {
            // Convert extractor matrices to temp location
            convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'extractor_matrices')
            input_files_channel = convert_data_workflow.out.input_files
            
            if (params.plot_input_spectra) {
                for (mutation_type in params.mutation_types) {
                    plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/input_tables", 'mutation_spectra')
                }
            }
        }
    } else {
        // Use default signature tables - collect them into a channel
        signature_files_channel = Channel.fromPath("${params.signature_tables}/*.csv").collect()
        
        if (params.plot_signatures) {
            for (mutation_type in params.mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, params.signature_tables, 'signatures')
            }
        }
    }

    if (params.SP_matrix_generator_output_path) {
        // Convert matrix generator output to temp location
        convert_data_workflow(params.dataset, params.SP_matrix_generator_output_path, 'matrix_generator_matrices')
        input_files_channel = convert_data_workflow.out.input_files
        
        if (params.plot_input_spectra) {
            for (mutation_type in params.mutation_types) {
                plot_spectra_workflow(params.dataset, mutation_type, "${params.temp_path}/input_tables", 'mutation_spectra')
            }
        }
    } else if (!params.SP_extractor_output_path) {
        // Use default input tables - collect them into a channel
        input_files_channel = Channel.fromPath("${params.input_tables}/**/*.csv").collect()
        
        if (params.plot_input_spectra) {
            for (mutation_type in params.mutation_types) {
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
        for (weak_threshold in weak_thresholds) {
            for (strong_threshold in strong_thresholds) {
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
                    number_of_bootstrapped_samples_in_optimisation
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
            weak_thresholds,
            strong_thresholds,
            number_of_bootstrapped_samples_in_optimisation
        )
        
        // Calculate optimal penalties
        OPTIMAL_PENALTIES_workflow(
            BOOTSTRAP_TABLES_workflow.out.bootstrap_tables,
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            weak_thresholds,
            strong_thresholds
        )
        
        // Generate optimization plots
        OPTIMISATION_PLOTS_workflow(
            OPTIMAL_PENALTIES_workflow.out.penalties_for_optimisation_plotting,
            weak_thresholds,
            strong_thresholds
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
            number_of_bootstrapped_samples
        )

        // Generate final bootstrap tables after final bootstrap NNLS completes
        FINAL_BOOTSTRAP_TABLES_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_NNLS_BOOTSTRAP_workflow.out.bootstrap_indices,
            number_of_bootstrapped_samples,
            suffix
        )

        // Generate final plots
        ALL_FINAL_PLOTS_workflow(
            FINAL_NNLS_workflow.out.dataset_mutation_pairs,
            FINAL_BOOTSTRAP_TABLES_workflow.out.attributions_per_sample,
            FINAL_BOOTSTRAP_TABLES_workflow.out.signature_prevalences
        )
    }
}

