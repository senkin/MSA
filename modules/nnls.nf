// modules/nnls.nf
//
// Every process stages the files it reads in the layout run_NNLS.py expects:
//   input_tables/<dataset>/WGS_<dataset>.<table>.csv     input mutation tables
//   simulated_tables/SIM_<dataset>/WGS_SIM_<dataset>...  simulated tables (and their truth)
//   signature_tables/<prefix>_<type>_signatures.csv      signature tables
// and writes its outputs in its own work directory only.

include { signature_prefix; gpu_flag; table_name; list_args } from './common'

// Unoptimised attribution of the input data, whose signature activities seed the simulations
process run_unoptimized_NNLS {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/temp/output_tables_unoptimised", mode: 'copy', overwrite: true, enabled: !params.cleanup_temp.toBoolean()

    input:
    tuple val(dataset), val(mutation_type)
    path "input_tables/${dataset}/*"
    path signature_files, stageAs: 'signature_tables/*'

    output:
    tuple val(dataset), val(mutation_type), path("${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv"), emit: mutations_table
    path "${dataset}/output_${dataset}_${mutation_type}_{weights_table,stat_info,fitted_values,residuals}.csv", emit: other_tables

    script:
    """
    run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
        -p ${signature_prefix()} -i input_tables -s signature_tables -o "./" \\
        -n ${params.number_of_samples}
    """
}

// Resolve the weak/strong thresholds scanned in the optimisation into the labels used in
// output names (written by run_NNLS.py, so every step names folders the same way).
// weak_thresholds = 'auto' picks them from removal trajectories of a subsample of the
// simulated data; an explicit list is only formatted.
process resolve_thresholds {
    tag "${mutation_type}/${dataset}"
    // a closure, since the path uses the inputs (as a plain string it fails with the strict parser)
    publishDir { "${params.output_path}/outputs_optimisation/SIM_${dataset}/${mutation_type}" }, mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("simulated_tables/SIM_${dataset}/*")
    path signature_files, stageAs: 'signature_tables/*'

    output:
    tuple val(dataset), val(mutation_type), path('weak_thresholds.txt'), path('strong_thresholds.txt')

    script:
    """
    run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
        -p ${signature_prefix()} -x --optimisation_strategy ${params.optimisation_strategy} \\
        -W ${list_args(params.weak_thresholds)} -S ${list_args(params.strong_thresholds)} \\
        --auto_thresholds ${params.auto_weak_thresholds_number} -n ${params.auto_weak_thresholds_samples} \\
        -i simulated_tables -s signature_tables \\
        -o "./" --resolve_thresholds_only --thresholds_dir ./
    """
}

// Optimised attribution of the simulated data: one task for all threshold pairs, one output
// folder SIM_<dataset>_<context>_NNLS_<weak>_<strong> per pair
process run_optimized_NNLS {
    tag "${mutation_type}/${dataset}"
    label 'gpu_nnls'
    publishDir "${params.output_path}/outputs_optimisation", mode: 'copy', overwrite: true, enabled: !params.cleanup_temp.toBoolean()

    input:
    tuple val(dataset), val(mutation_type), path("simulated_tables/SIM_${dataset}/*"), path(weak_thresholds_file), path(strong_thresholds_file)
    path signature_files, stageAs: 'signature_tables/*'

    output:
    tuple val(dataset), val(mutation_type), path("SIM_${dataset}_${params.SBS_context}_NNLS_*/*.csv"), emit: tables

    script:
    def truth = "simulated_tables/SIM_${dataset}/WGS_SIM_${dataset}.${table_name(mutation_type)}.weights.csv"
    """
    # Run NNLS analysis for all threshold pairs
    run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
        -p ${signature_prefix()} --optimisation_strategy ${params.optimisation_strategy} \\
        -W \$(cat ${weak_thresholds_file}) -S \$(cat ${strong_thresholds_file}) --threshold_scan ${params.threshold_scan} \\
        -i simulated_tables -s signature_tables \\
        -o "./" -x --add_suffix ${gpu_flag()}

    # Every threshold pair's folder also holds the simulated truth, which the tables step scores against
    for folder in SIM_${dataset}_${params.SBS_context}_NNLS_*/; do
        cp ${truth} "\$folder"
    done
    """
}

// Bootstrap of the optimised attribution: each task runs one chunk of iterations in-process
// via run_NNLS.py --n_bootstrap, for all threshold pairs, writing one indexed output set per
// iteration and pair (indices start_index .. start_index + n_per_task - 1) into each pair's
// bootstrap_output/ folder. Chunks come from bootstrap_chunks().
process run_bootstrap_NNLS {
    tag "${mutation_type}/${dataset}/${chunk[0]}"
    label 'gpu_nnls'
    publishDir "${params.output_path}/outputs_optimisation", enabled: !params.cleanup_temp.toBoolean()

    input:
    tuple val(dataset), val(mutation_type), path("simulated_tables/SIM_${dataset}/*"), path(weak_thresholds_file), path(strong_thresholds_file)
    path signature_files, stageAs: 'signature_tables/*'
    each chunk              // [start_index, n_per_task]
    val seed

    output:
    tuple val(dataset), val(mutation_type), path("SIM_${dataset}_${params.SBS_context}_NNLS_*/bootstrap_output/*.csv"), emit: tables

    script:
    def start_index = chunk[0]
    def n_per_task = chunk[1]
    """
    run_NNLS.py -B --n_bootstrap ${n_per_task} --bootstrap_start_index ${start_index} \\
        -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \\
        --optimisation_strategy ${params.optimisation_strategy} \\
        --bootstrap_method ${params.bootstrap_method} \\
        -W \$(cat ${weak_thresholds_file}) -S \$(cat ${strong_thresholds_file}) --threshold_scan ${params.threshold_scan} \\
        --add_suffix -p ${signature_prefix()} -i simulated_tables \\
        -s signature_tables -o "./" ${gpu_flag()} --seed ${seed}
    """
}

// Final attribution of the input data with the optimal penalties
process run_final_NNLS {
    tag "${mutation_type}/${dataset}"
    label 'gpu_nnls'
    publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path(weak_penalty), path(strong_penalty)
    path "input_tables/${dataset}/*"
    path signature_files, stageAs: 'signature_tables/*'

    output:
    tuple val(dataset), val(mutation_type), path("${dataset}/output_${dataset}_${mutation_type}_*.csv"), emit: tables
    path "${dataset}/gpu_provenance.json", optional: true, emit: gpu_provenance

    script:
    def optimised_flag = params.optimised.toBoolean() ? "-x" : ""
    """
    run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
        --optimisation_strategy ${params.optimisation_strategy} \\
        -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
        -p ${signature_prefix()} -i input_tables -s signature_tables -o "./" ${gpu_flag()}
    """
}

// Bootstrap of the final attribution, in chunks of iterations (see run_bootstrap_NNLS)
process run_final_bootstrap_NNLS {
    tag "${mutation_type}/${dataset}/${chunk[0]}"
    label 'gpu_nnls'
    publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true, enabled: !params.cleanup_temp.toBoolean()

    input:
    tuple val(dataset), val(mutation_type), path(weak_penalty), path(strong_penalty)
    path "input_tables/${dataset}/*"
    path signature_files, stageAs: 'signature_tables/*'
    each chunk                   // [start_index, n_per_task]
    val seed

    output:
    tuple val(dataset), val(mutation_type), path("${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*.csv"), emit: tables

    script:
    def start_index = chunk[0]
    def n_per_task = chunk[1]
    def optimised_flag = params.optimised.toBoolean() ? "-x" : ""
    """
    run_NNLS.py -B --n_bootstrap ${n_per_task} --bootstrap_start_index ${start_index} \\
        -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
        --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \\
        -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
        -p ${signature_prefix()} -i input_tables -s signature_tables -o "./" ${gpu_flag()} --seed ${seed}
    """
}
