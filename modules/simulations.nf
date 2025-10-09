// modules/simulations.nf

// parameters for helper flags
// params.number_of_simulated_samples = -1
// params.noise_type = "gaussian"
// params.noise_stdev = 10
// params.zero_inflation_threshold = 0.05
// params.SBS_context = 96
def noise_flag = params.add_noise ? "-z" : ""
def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow simulate_data_workflow {
    take:
    dataset_mutation_pairs // tuple of (dataset, mutation_type)
    mutations_table // path to mutations table
    // signature_files // path to signature files

    main:
    process simulate_data {
        tag "${dataset}_${mutation_type}"
        publishDir "${params.temp_path}/output_tables", mode: 'copy', overwrite: true

        input:
        tuple val(dataset), val(mutation_type)
        path mutations_table
        // path signature_files

        output:
        path "SIM_${dataset}/WGS_SIM_${dataset}.*.csv", optional: true
        tuple val(dataset), val(mutation_type), emit: simulation_outputs
        tuple val(dataset), val(mutation_type), emit: simulation_outputs_for_bootstrap

        script:
        """
        python ${workflow.projectDir}/bin/simulate_data.py \
            -d SIM_${dataset} \
            -t ${mutation_type} \
            -c ${params.SBS_context} \
            -n ${params.number_of_simulated_samples} \
            -p ${signature_prefix} \
            -B ${noise_flag} \
            --noise_type ${params.noise_type} \
            -Z ${params.noise_stdev} \
            --zero_inflation_threshold ${params.zero_inflation_threshold} \
            -i ${mutations_table} \
            -s ${params.temp_path}/signature_tables \
            -o "./"
        """
    }

    // Execute simulation for each dataset/mutation_type pair
    simulate_data(
        dataset_mutation_pairs,
        mutations_table
    )

    emit:
    simulation_outputs = simulate_data.out.simulation_outputs  // (dataset, mutation_type) channel
    simulation_outputs_for_bootstrap = simulate_data.out.simulation_outputs_for_bootstrap  // (dataset, mutation_type) channel
    all_simulation_outputs = simulate_data.out[0]  // All output files
}