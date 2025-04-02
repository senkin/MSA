// modules/nnls.nf

// parameters for helper flags
params.optimised = false
params.bootstrap_samples = 1

workflow NNLS_workflow {
    take:
    dataset
    mutation_type
    input_files
    signature_files
    weak_thresholds
    strong_thresholds
    bootstrap_samples

    main:
    // Create channels from input parameters
    ch_weak_thresholds = weak_thresholds instanceof Channel ? weak_thresholds : Channel.from(weak_thresholds)
    ch_strong_thresholds = strong_thresholds instanceof Channel ? strong_thresholds : Channel.from(strong_thresholds)
    ch_bootstrap_samples = bootstrap_samples instanceof Channel ? bootstrap_samples : Channel.from(1..bootstrap_samples)

    // Unified NNLS process
    process run_NNLS {
        tag "${mutation_type}/${dataset}"
        
        publishDir (
            params.optimised ? "${params.optimisation_NNLS_output_path}" : "$baseDir/output_tables_unoptimised",
            mode: 'copy',
            overwrite: true,
            saveAs: { filename ->
                if (params.optimised) {
                    "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/$filename"
                } else {
                    "${dataset}/$filename"
                }
            }
        )

        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        val weak_threshold
        val strong_threshold
        val bootstrap_sample

        output:
        path "./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv", emit: mutations_table
        path "./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv", emit: weights_table
        path "./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv", emit: stat_info
        path "./${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv", optional: true, emit: fitted_values
        path "./${dataset}/output_${dataset}_${mutation_type}_residuals.csv", optional: true, emit: residuals

        script:
        """
        python $baseDir/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \
                      -p ${params.signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./" \
                      -n ${params.number_of_samples} \
                      ${params.optimised ? "-x --optimisation_strategy ${params.optimisation_strategy} -W ${weak_threshold} -S ${strong_threshold}" : ""} \
                      ${bootstrap_sample > 1 ? "-B --bootstrap_method ${params.bootstrap_method}" : ""}
        """
    }

    // Prepare input channels
    ch_dataset = Channel.from(dataset)
    ch_mutation_type = Channel.from(mutation_type)
    ch_input_files = input_files
    ch_signature_files = signature_files

    if (params.optimised) {
        ch_weak_thresholds
            .combine(ch_strong_thresholds)
            .combine(ch_bootstrap_samples)
            .set { optimization_runs }
        
        run_NNLS(
            ch_dataset,
            ch_mutation_type,
            ch_input_files,
            ch_signature_files,
            optimization_runs.map { it[0] },
            optimization_runs.map { it[1] },
            optimization_runs.map { it[2] }
        )
    } else {
        run_NNLS(
            ch_dataset,
            ch_mutation_type,
            ch_input_files,
            ch_signature_files,
            Channel.from('0.0000'),
            Channel.from('0.0000'),
            Channel.from(1)
        )
    }

    emit:
    mutations_table = run_NNLS.out.mutations_table
    weights_table = run_NNLS.out.weights_table
    stat_info = run_NNLS.out.stat_info
    fitted_values = run_NNLS.out.fitted_values
    residuals = run_NNLS.out.residuals
}