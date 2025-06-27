// modules/nnls.nf

// Unoptimized NNLS workflow
workflow NNLS_unoptimized_workflow {
    take:
    dataset
    mutation_type
    input_files
    signature_files
    
    main:
    // Unoptimized NNLS process
    process run_unoptimized_NNLS {
        tag "${mutation_type}/${dataset}"
        publishDir "$baseDir/output_tables_unoptimised", mode: 'copy', overwrite: true, saveAs: { filename -> "${dataset}/${filename}" }
        
        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        
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
            -n ${params.number_of_samples}
        """
    }
    
    // Run the unoptimized process
    run_unoptimized_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files
    )
    
    emit:
    dataset_mutation_pairs = run_unoptimized_NNLS.out.dataset_mutation_pairs
    mutations_table = run_unoptimized_NNLS.out.mutations_table
    weights_table = run_unoptimized_NNLS.out.weights_table
    stat_info = run_unoptimized_NNLS.out.stat_info
    fitted_values = run_unoptimized_NNLS.out.fitted_values
    residuals = run_unoptimized_NNLS.out.residuals
}

// Optimized NNLS workflow (single run per threshold combination)
workflow NNLS_optimized_workflow {
    take:
    dataset
    mutation_type
    input_files
    signature_files
    weak_threshold
    strong_threshold
    
    main:
    // Optimized NNLS process
    process run_optimized_NNLS {
        tag "${mutation_type}/${dataset}/${weak_threshold}/${strong_threshold}"
        publishDir "${params.optimisation_NNLS_output_path}", mode: 'copy', overwrite: true
        
        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        val weak_threshold
        val strong_threshold
        
        output:
        path "SIM_${dataset}_${params.SBS_context}_NNLS_*/output_SIM_${dataset}_${mutation_type}_mutations_table.csv", emit: mutations_table
        path "SIM_${dataset}_${params.SBS_context}_NNLS_*/output_SIM_${dataset}_${mutation_type}_weights_table.csv", emit: weights_table
        path "SIM_${dataset}_${params.SBS_context}_NNLS_*/output_SIM_${dataset}_${mutation_type}_stat_info.csv", emit: stat_info
        path "SIM_${dataset}_${params.SBS_context}_NNLS_*/output_SIM_${dataset}_${mutation_type}_fitted_values.csv", optional: true, emit: fitted_values
        path "SIM_${dataset}_${params.SBS_context}_NNLS_*/output_SIM_${dataset}_${mutation_type}_residuals.csv", optional: true, emit: residuals
        tuple val(dataset), val(mutation_type), emit: dataset_mutation_pairs
        
        when:
        !params.run_only_simulations
        
        script:
        """
        # Create output directory
        mkdir -p ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}
        
        # Copy the appropriate simulated dataset based on mutation type
        if [[ ${mutation_type} == "SBS" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi
        
        # Run NNLS analysis
        python $baseDir/bin/run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${params.signature_prefix} --optimisation_strategy ${params.optimisation_strategy} \\
            -W ${weak_threshold} -S ${strong_threshold} \\
            -i $baseDir/output_tables -s ${params.signature_tables} \\
            -o "./" -x --add_suffix
        """
    }
    
    // Run the optimized process
    run_optimized_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files,
        weak_threshold,
        strong_threshold
    )
    
    emit:
    dataset_mutation_pairs = run_optimized_NNLS.out.dataset_mutation_pairs
    mutations_table = run_optimized_NNLS.out.mutations_table
    weights_table = run_optimized_NNLS.out.weights_table
    stat_info = run_optimized_NNLS.out.stat_info
    fitted_values = run_optimized_NNLS.out.fitted_values
    residuals = run_optimized_NNLS.out.residuals
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
        
        when:
        !params.run_only_simulations
        
        script:
        """
        # Create output directory
        mkdir -p ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}
        
        # Copy the appropriate simulated dataset based on mutation type
        if [[ ${mutation_type} == "SBS" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi
        
        # Run bootstrap NNLS analysis
        python $baseDir/bin/run_NNLS.py -B -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            --bootstrap_method ${params.bootstrap_method} \\
            -W ${weak_threshold} -S ${strong_threshold} --add_suffix \\
            -p ${params.signature_prefix} -i $baseDir/output_tables -s ${params.signature_tables} -o "./"
        
        # Create bootstrap output directory and move files
        mkdir -p SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output
        
        # Move the output files with bootstrap index
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