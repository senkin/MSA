// modules/simulations.nf

include { signature_prefix } from './common'

// Simulated data (SIM_<dataset>) from the unoptimised signature activities of the input data:
// mutation tables and their truth (signature weights) table
process simulate_data {
    tag "${dataset}_${mutation_type}"
    // intermediate files, kept in <output_path>/temp with cleanup_temp = false or when only simulating
    publishDir "${params.output_path}/temp/output_tables", mode: 'copy', overwrite: true,
        enabled: !params.cleanup_temp.toString().toBoolean() || params.run_only_simulations.toString().toBoolean()

    input:
    tuple val(dataset), val(mutation_type), path(mutations_table)
    path signature_files, stageAs: 'signature_tables/*'
    val seed

    output:
    tuple val(dataset), val(mutation_type), path("SIM_${dataset}/WGS_SIM_${dataset}.*.csv"), emit: tables

    script:
    def noise_flag = params.add_noise.toString().toBoolean() ? "-z" : ""
    """
    simulate_data.py \\
        -d SIM_${dataset} \\
        -t ${mutation_type} \\
        -c ${params.SBS_context} \\
        -n ${params.number_of_simulated_samples} \\
        -p ${signature_prefix()} \\
        -B ${noise_flag} \\
        --noise_type ${params.noise_type} \\
        -Z ${params.noise_stdev} \\
        --zero_inflation_threshold ${params.zero_inflation_threshold} \\
        -i ${mutations_table} \\
        -s signature_tables \\
        -o "./" \\
        --seed ${seed}
    """
}
