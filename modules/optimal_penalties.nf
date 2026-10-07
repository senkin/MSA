// modules/optimal_penalties.nf

include { signature_prefix; list_args } from './common'

// Optimal penalties of a mutation type from the truth studies of every threshold pair
process calculate_optimal_penalties {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/outputs_optimisation", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("simulated_tables/SIM_${dataset}/*"), path(weak_thresholds_file), path(strong_thresholds_file),
        path(pairs, stageAs: 'attributions/*')     // bundle/<pair> of every threshold pair
    path signature_files, stageAs: 'signature_tables/*'

    output:
    tuple val(dataset), val(mutation_type), path("SIM_${dataset}/${mutation_type}/optimal_weak_penalty"), path("SIM_${dataset}/${mutation_type}/optimal_strong_penalty"), emit: penalties
    tuple val(dataset), val(mutation_type), path("SIM_${dataset}/${mutation_type}/*.*"), emit: analysis   // tables (csv, json)

    script:
    def no_CI_for_penalties_flag = params.no_CI_for_penalties.toString().toBoolean() ? "--no_CI" : ""
    def calculate_penalty_on_average_flag = params.calculate_penalty_on_average.toString().toBoolean() ? "--average" : ""
    def prioritised_signatures_flag = params.signatures_to_prioritise ? "--signatures_to_prioritise " + list_args(params.signatures_to_prioritise) : ''
    def deprioritised_signatures_flag = params.signatures_to_deprioritise ? "--signatures_to_deprioritise " + list_args(params.signatures_to_deprioritise) : ''
    """
    calculate_optimal_penalties.py -d SIM_${dataset} -t ${mutation_type} \\
        -I simulated_tables/SIM_${dataset} \\
        -i attributions -o "./" \\
        -c ${params.SBS_context} \\
        ${no_CI_for_penalties_flag} \\
        ${calculate_penalty_on_average_flag} \\
        ${prioritised_signatures_flag} \\
        ${deprioritised_signatures_flag} \\
        -M ${params.metric_to_prioritise} \\
        -T ${params.metric_threshold} \\
        -W \$(cat ${weak_thresholds_file}) \\
        -S \$(cat ${strong_thresholds_file}) \\
        --signature_path signature_tables \\
        -p ${signature_prefix()}
    """
}
