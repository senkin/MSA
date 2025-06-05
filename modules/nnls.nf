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
    weak_threshold
    strong_threshold
    bootstrap_samples
    
    main:
    // Unified NNLS process
    process run_NNLS {
        tag "${mutation_type}/${dataset}"
        publishDir (
            params.optimised ? "${params.optimisation_NNLS_output_path}" : "$baseDir/output_tables_unoptimised",
            mode: 'copy',
            overwrite: true,
            saveAs: { filename ->
                if (params.optimised) {
                    "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/${filename}"
                } else {
                    "${dataset}/${filename}"
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
        val bootstrap_samples
        
        output:
        path "./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv", emit: mutations_table
        path "./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv", emit: weights_table
        path "./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv", emit: stat_info
        path "./${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv", optional: true, emit: fitted_values
        path "./${dataset}/output_${dataset}_${mutation_type}_residuals.csv", optional: true, emit: residuals
        tuple val(dataset), val(mutation_type), emit: dataset_mutation_pairs
        
        script:
        """
        python $baseDir/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${params.signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./" \\
            -n ${params.number_of_samples} \\
            ${params.optimised ? "-x --optimisation_strategy ${params.optimisation_strategy} -W ${weak_threshold} -S ${strong_threshold}" : ""} \\
            ${bootstrap_samples > 1 ? "-B --bootstrap_method ${params.bootstrap_method}" : ""}
        """
    }
    
    // Run the process directly with the input values
    run_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files,
        weak_threshold,
        strong_threshold,
        bootstrap_samples
    )
    
    emit:
    dataset_mutation_pairs = run_NNLS.out.dataset_mutation_pairs
    mutations_table = run_NNLS.out.mutations_table
    weights_table = run_NNLS.out.weights_table
    stat_info = run_NNLS.out.stat_info
    fitted_values = run_NNLS.out.fitted_values
    residuals = run_NNLS.out.residuals
}

// Bootstrap-specific workflow that runs multiple iterations
workflow NNLS_bootstrap_workflow {
    take:
    dataset
    mutation_type
    input_files
    signature_files
    weak_threshold
    strong_threshold
    num_bootstrap_samples
    
    main:
    // Bootstrap NNLS process - runs once per bootstrap sample
    process run_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${bootstrap_index}"
        publishDir "${params.optimisation_NNLS_output_path}"
        
        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        val weak_threshold
        val strong_threshold
        val bootstrap_index
        
        output:
        path "./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_mutations_table.csv", emit: mutations_table
        path "./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_weights_table.csv", emit: weights_table
        path "./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_stat_info.csv", emit: stat_info
        val bootstrap_index, emit: bootstrap_indices
        
        script:
        """
        python $baseDir/bin/run_NNLS.py -B -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            --bootstrap_method ${params.bootstrap_method} \\
            -W ${weak_threshold} -S ${strong_threshold} --add_suffix \\
            -p ${params.signature_prefix} -i $baseDir/output_tables -s ${params.signature_tables} -o "./"
        
        mkdir -p SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output
        
        mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_mutations_table.csv \\
           SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_mutations_table.csv
        
        mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_weights_table.csv \\
           SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_weights_table.csv
        
        mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_stat_info.csv \\
           SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${bootstrap_index}_stat_info.csv
        """
    }
    
    // Create a channel with bootstrap indices (1 to num_bootstrap_samples)
    bootstrap_indices_ch = Channel.from(1..num_bootstrap_samples)
    
    // Run bootstrap process for each index
    run_bootstrap_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files,
        weak_threshold,
        strong_threshold,
        bootstrap_indices_ch
    )
    
    emit:
    mutations_tables = run_bootstrap_NNLS.out.mutations_table
    weights_tables = run_bootstrap_NNLS.out.weights_table
    stat_infos = run_bootstrap_NNLS.out.stat_info
    bootstrap_indices = run_bootstrap_NNLS.out.bootstrap_indices
}