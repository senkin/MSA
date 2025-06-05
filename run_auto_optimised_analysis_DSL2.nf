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
params.SP_extractor_output_path = null // optional path to SigProfilerExtractor output
params.SP_matrix_generator_output_path = null // optional path to SigProfilerMatrixGenerator output
params.COSMIC_signatures = false // if set to true, COSMIC signatures are used form SigProfiler output, otherwise de-novo ones are used
params.dataset = 'SIM_test' // dataset name. Input matrices to be provided in params.input_tables/params.dataset, unless SigProfiler inputs are used
params.mutation_types = ['DBS'] // add or remove mutation types if needed
params.input_tables = "$baseDir/input_mutation_tables"
params.SBS_context = 96 // 96, 192, 288, 1536 context matrices can be provided (SBS only)
params.number_of_samples = -1 // number of samples to analyse (-1 means all available)

// output paths
params.output_path = "."
params.tables_output_path = params.output_path + "/output_tables"
params.plots_output_path = params.output_path + "/plots"

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
params.optimisation_NNLS_output_path = "$baseDir/outputs_optimisation"
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
prioritised_signatures_flag = (params.signatures_to_prioritise) ? "--signatures_to_prioritise " + params.signatures_to_prioritise.join(' ') : ''
// override number of samples/variations for test run
test_run = (params.dataset == 'SIM_test') ? true : false
number_of_bootstrapped_samples_in_optimisation = (test_run) ? 10 : params.number_of_bootstrapped_samples_in_optimisation
number_of_bootstrapped_samples = (test_run) ? 10 : params.number_of_bootstrapped_samples
number_of_simulated_samples = (test_run) ? 10 : params.number_of_simulated_samples
weak_thresholds = (test_run) ? ['0.0000', '0.0100', '0.0200'] : params.weak_thresholds
strong_thresholds = (test_run) ? ['0.0000'] : params.strong_thresholds
mutation_types = (test_run) ? ['SBS', 'DBS', 'ID'] : params.mutation_types

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

// Include the plotting module and pass helper flags
include { plot_spectra_workflow } from './modules/plot_spectra' addParams(
    strands_flag: strands_flag,
    nontranscribed_flag: nontranscribed_flag,
    signature_prefix: signature_prefix
)

// Include the SP -> MSA data conversion module and pass helper flags
include { convert_data_workflow } from './modules/data_conversion' addParams(
    strands_flag: strands_flag,
    nontranscribed_flag: nontranscribed_flag,
    error_flag: error_flag,
    COSMIC_flag: COSMIC_flag,
    signature_prefix: signature_prefix
)

// include NNLS workflows
include { NNLS_workflow as UnoptimizedNNLS } from './modules/nnls'
include { NNLS_workflow as OptimizedNNLS } from './modules/nnls' addParams(optimised: true)
include { NNLS_bootstrap_workflow as OptimizedNNLSforBootstrap } from './modules/nnls' addParams(optimised: true)
include { simulate_data_workflow } from './modules/simulations'

// Main workflow
workflow {
    // Create placeholder channels if SP_extractor_output_path and SP_matrix_generator_output_path are null
    if (params.SP_extractor_output_path == null && params.SP_matrix_generator_output_path == null) {
        // Create a dummy input directory
        def defaultInputDir = file("${params.signature_tables}")
        defaultInputDir.mkdirs()

        // Emit the input directory as a channel
        Channel.fromPath(defaultInputDir).set { converted_SP_to_MSA_for_spectra }
        Channel.fromPath(defaultInputDir).set { converted_SP_to_MSA_for_unoptimised_NNLS }
        Channel.fromPath(defaultInputDir).set { signatures_for_spectra }
        Channel.fromPath(defaultInputDir).set { signatures_for_unoptimised_NNLS }
    }

    def converted_SP_to_MSA_for_unoptimised_NNLS = null
    def converted_SP_to_MSA_for_spectra = null
    def signatures_for_spectra = null
    def signatures_for_unoptimised_NNLS = null

    if (params.SP_extractor_output_path) {
        convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'signature_tables')
        signatures_for_spectra = convert_data_workflow.out.signatures_for_spectra
        signatures_for_unoptimised_NNLS = convert_data_workflow.out.signatures_for_unoptimised_NNLS
        if (params.plot_signatures) {
            params.mutation_types.each { mutation_type ->
                plot_spectra_workflow(params.dataset, mutation_type, signatures_for_spectra, 'signatures')
        }
        }
        if (!params.SP_matrix_generator_output_path) {
            convert_data_workflow(params.dataset, params.SP_extractor_output_path, 'extractor_matrices')
            converted_SP_to_MSA_for_spectra = convert_data_workflow.out.converted_SP_to_MSA_for_spectra
            converted_SP_to_MSA_for_unoptimised_NNLS = convert_data_workflow.out.converted_SP_to_MSA_for_unoptimised_NNLS
            if (params.plot_input_spectra) {
                params.mutation_types.each { mutation_type ->
                plot_spectra_workflow(params.dataset, mutation_type, converted_SP_to_MSA_for_spectra, 'mutation_spectra')
                }
            }
        }
    } else if (params.plot_signatures) {
        params.mutation_types.each { mutation_type ->
            plot_spectra_workflow(params.dataset, mutation_type, signatures_for_spectra, 'signatures')
      }
    }

    if (params.SP_matrix_generator_output_path) {
        convert_data_workflow(params.dataset, params.SP_matrix_generator_output_path, 'matrix_generator_matrices')
        converted_SP_to_MSA_for_spectra = convert_data_workflow.out.converted_SP_to_MSA_for_spectra
        converted_SP_to_MSA_for_unoptimised_NNLS = convert_data_workflow.out.converted_SP_to_MSA_for_unoptimised_NNLS
        if (params.plot_input_spectra) {
            params.mutation_types.each { mutation_type ->
                plot_spectra_workflow(params.dataset, mutation_type, convert_data_workflow.out.converted_SP_to_MSA_for_spectra, 'mutation_spectra')
            }
        }
    } else if ((params.plot_input_spectra) && (!params.SP_extractor_output_path)) {
        params.mutation_types.each { mutation_type ->
            plot_spectra_workflow(params.dataset, mutation_type, converted_SP_to_MSA_for_spectra, 'mutation_spectra')
        }
    }
    // Initial unoptimized NNLS runs
    params.mutation_types.each { mutation_type ->
        // Run unoptimized NNLS
        UnoptimizedNNLS(
            params.dataset,
            mutation_type,
            converted_SP_to_MSA_for_unoptimised_NNLS,
            signatures_for_unoptimised_NNLS,
            0.0,      // weak_threshold (not used for unoptimized)
            0.0,      // strong_threshold (not used for unoptimized)
            1         // bootstrap_samples
        )
        
        // Run simulation workflow
        simulate_data_workflow(
            UnoptimizedNNLS.out.dataset_mutation_pairs,
            UnoptimizedNNLS.out.mutations_table
        )
        
        // Run optimized NNLS for each threshold combination
        params.weak_thresholds.each { weak_threshold ->
            params.strong_thresholds.each { strong_threshold ->
                // Standard optimized NNLS
                OptimizedNNLS(
                    params.dataset,
                    mutation_type,
                    simulate_data_workflow.out.all_simulation_outputs,
                    signatures_for_unoptimised_NNLS,
                    weak_threshold,
                    strong_threshold,
                    1
                )
                
                // Optimized NNLS with bootstrap
                OptimizedNNLSforBootstrap(
                    params.dataset,
                    mutation_type,
                    simulate_data_workflow.out.all_simulation_outputs,
                    signatures_for_unoptimised_NNLS,
                    weak_threshold,
                    strong_threshold,
                    number_of_bootstrapped_samples_in_optimisation
                )
            }
        }
    }
}

