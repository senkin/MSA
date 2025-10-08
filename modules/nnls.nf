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
        publishDir "${params.temp_path}/output_tables_unoptimised", mode: 'copy', overwrite: true, saveAs: { filename -> "${dataset}/${filename}" }
        
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
        mkdir -p ${params.temp_path}
        cp -a ${params.signature_tables} ${params.temp_path}
        python ${workflow.projectDir}/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${params.signature_prefix} -i ${params.temp_path}/input_tables -s ${params.temp_path}/signature_tables -o "./" \\
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
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi
        
        # Run NNLS analysis
        python ${workflow.projectDir}/bin/run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${params.signature_prefix} --optimisation_strategy ${params.optimisation_strategy} \\
            -W ${weak_threshold} -S ${strong_threshold} \\
            -i ${params.temp_path}/output_tables -s ${params.temp_path}/signature_tables \\
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
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp ${params.temp_path}/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi
        
        # Run bootstrap NNLS analysis
        python ${workflow.projectDir}/bin/run_NNLS.py -B -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            --bootstrap_method ${params.bootstrap_method} \\
            -W ${weak_threshold} -S ${strong_threshold} --add_suffix \\
            -p ${params.signature_prefix} -i ${params.temp_path}/output_tables -s ${params.temp_path}/signature_tables -o "./"
        
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

// Final NNLS workflow using optimal penalties (reuses existing process logic)
workflow FINAL_NNLS_workflow {
    take:
    penalties_for_attribution  // tuple of (dataset, mutation_type)
    optimal_weak_penalty_files // optimal weak penalty files  
    optimal_strong_penalty_files // optimal strong penalty files
    input_files
    signature_files
    
    main:
    // Simple final NNLS process using penalty files directly
    process run_final_NNLS {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.tables_output_path}", mode: 'copy', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        path weak_penalty
        path strong_penalty
        path input_files
        path signature_files
        
        output:
        path "./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv", emit: mutations_table
        path "./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv", emit: weights_table
        path "./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv", emit: stat_info
        path "./${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv", emit: fitted_values
        path "./${dataset}/output_${dataset}_${mutation_type}_residuals.csv", emit: residuals
        tuple val(dataset), val(mutation_type), emit: dataset_mutation_pairs
        
        when:
        !params.run_only_optimisation
        
        script:
        def optimised_flag = params.optimised ? "-x" : ""
        """
        python ${workflow.projectDir}/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
            -p ${params.signature_prefix} -i ${params.temp_path}/input_tables -s ${params.temp_path}/signature_tables -o "./"
        
        # Copy residuals and fitted values back to input tables
        cp ${dataset}/output_${dataset}_${mutation_type}_residuals.csv ${params.temp_path}/input_tables/${dataset}/
        cp ${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv ${params.temp_path}/input_tables/${dataset}/
        """
    }
    
    run_final_NNLS(
        penalties_for_attribution,
        optimal_weak_penalty_files,
        optimal_strong_penalty_files,
        input_files,
        signature_files
    )
    
    emit:
    dataset_mutation_pairs = run_final_NNLS.out.dataset_mutation_pairs
    mutations_table = run_final_NNLS.out.mutations_table
    weights_table = run_final_NNLS.out.weights_table
    stat_info = run_final_NNLS.out.stat_info
    fitted_values = run_final_NNLS.out.fitted_values
    residuals = run_final_NNLS.out.residuals
}

// Final bootstrap workflow using optimal penalties
workflow FINAL_NNLS_BOOTSTRAP_workflow {
    take:
    penalties_for_attribution  // tuple of (dataset, mutation_type)
    optimal_weak_penalty_files // optimal weak penalty files
    optimal_strong_penalty_files // optimal strong penalty files
    input_files
    signature_files
    num_bootstrap_samples
    
    main:
    // Simple bootstrap process using penalty files directly
    process run_final_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${bootstrap_index}"
        publishDir "${params.tables_output_path}", mode: 'copy', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        path weak_penalty
        path strong_penalty
        path input_files
        path signature_files
        each bootstrap_index
        
        output:
        path "./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_mutations_table.csv", emit: mutations_table
        path "./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_stat_info.csv", emit: stat_info
        path "./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_weights_table.csv", emit: weights_table
        val bootstrap_index, emit: bootstrap_indices
        
        when:
        !params.run_only_optimisation
        
        script:
        def optimised_flag = params.optimised ? "-x" : ""
        """
        python ${workflow.projectDir}/bin/run_NNLS.py -B -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
            --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \\
            -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
            -p ${params.signature_prefix} -i ${params.temp_path}/input_tables -s ${params.temp_path}/signature_tables -o "./"
        
        mkdir -p ${dataset}/bootstrap_output
        mv ${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_mutations_table.csv
        mv ${dataset}/output_${dataset}_${mutation_type}_weights_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_weights_table.csv
        mv ${dataset}/output_${dataset}_${mutation_type}_stat_info.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${bootstrap_index}_stat_info.csv
        """
    }
    
    run_final_bootstrap_NNLS(
        penalties_for_attribution,
        optimal_weak_penalty_files,
        optimal_strong_penalty_files,
        input_files,
        signature_files,
        Channel.from(1..num_bootstrap_samples)
    )
    
    emit:
    mutations_tables = run_final_bootstrap_NNLS.out.mutations_table
    stat_infos = run_final_bootstrap_NNLS.out.stat_info
    weights_tables = run_final_bootstrap_NNLS.out.weights_table
    bootstrap_indices = run_final_bootstrap_NNLS.out.bootstrap_indices
}