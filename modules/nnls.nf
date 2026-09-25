// modules/nnls.nf

// In GPU mode, run_NNLS.py runs the greedy signature optimisation transposed
// across samples (one large masked NNLS batch per greedy step), offloaded to the
// cuML batched NNLS solver via batched_solve_and_score().
def gpu_flag = params.use_GPU ? "--use_gpu --gpu_precision ${params.gpu_precision} --gpu_batch_size ${params.gpu_batch_size}" : ''

// Split bootstrap iterations 1..n into [start_index, n_iterations] chunks, one chunk per task.
// run_NNLS.py seeds every iteration from (seed, index), so results do not depend on the split.
// chunk_size <= 0 keeps the mode default: one iteration per task on CPU (maximum scheduler
// width), all iterations in a single task on GPU (one reused CUDA context).
def bootstrap_chunks(n, chunk_size) {
    int size = (chunk_size ?: 0) > 0 ? chunk_size as int : (params.use_GPU ? n as int : 1)
    return (1..n).step(size).collect { start -> [start, Math.min(size, (n as int) - start + 1)] }
}

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
        publishDir "${params.output_path}/temp/output_tables_unoptimised", mode: 'copy', overwrite: true, saveAs: { filename -> "${dataset}/${filename}" }
        
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
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        mkdir -p ${params.output_path}/temp
        cp -a ${params.signature_tables} ${params.output_path}/temp
        python ${workflow.projectDir}/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${signature_prefix} -i ${params.output_path}/temp/input_tables -s ${params.output_path}/temp/signature_tables -o "./" \\
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
        label 'gpu_nnls'
        publishDir "${params.output_path}/outputs_optimisation", mode: 'copy', overwrite: true
        
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
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        # Create output directory
        mkdir -p ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}
        
        # Copy the appropriate simulated dataset based on mutation type
        if [[ ${mutation_type} == "SBS" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi
        
        # Run NNLS analysis
        python ${workflow.projectDir}/bin/run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
            -p ${signature_prefix} --optimisation_strategy ${params.optimisation_strategy} \\
            -W ${weak_threshold} -S ${strong_threshold} \\
            -i ${params.output_path}/temp/output_tables -s ${params.output_path}/temp/signature_tables \\
            -o "./" -x --add_suffix ${gpu_flag}
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
    seed                    // run-wide RNG seed (resolved in the main script)
    
    main:
    // Bootstrap NNLS process. Each task runs one chunk of bootstrap iterations
    // in-process via run_NNLS.py --n_bootstrap, writing one indexed output set per
    // iteration (indices start_index .. start_index + n_per_task - 1) directly into
    // the bootstrap_output/ folder. Chunks come from bootstrap_chunks() (see below).
    process run_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${weak_threshold}/${strong_threshold}/${chunk[0]}"
        label 'gpu_nnls'
        publishDir "${params.output_path}/outputs_optimisation"

        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        val weak_threshold
        val strong_threshold
        each chunk              // [start_index, n_per_task]
        val seed

        output:
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_mutations_table.csv", emit: mutations_tables
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_weights_table.csv", emit: weights_tables
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_stat_info.csv", emit: stat_infos

        when:
        !params.run_only_simulations

        script:
        def (start_index, n_per_task) = chunk
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        # Create output directory
        mkdir -p ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}

        # Copy the appropriate simulated dataset (truth weights) based on mutation type
        if [[ ${mutation_type} == "SBS" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "DBS" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "ID" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
            cp ${params.output_path}/temp/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv \\
               ${params.output_path}/outputs_optimisation/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
        fi

        # Run this task's slice of bootstrap iterations in-process. run_NNLS.py writes
        # one indexed output set per iteration directly into the bootstrap_output/ folder.
        python ${workflow.projectDir}/bin/run_NNLS.py -B --n_bootstrap ${n_per_task} --bootstrap_start_index ${start_index} \\
            --bootstrap_output_suffix ${weak_threshold}_${strong_threshold} \\
            -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            --bootstrap_method ${params.bootstrap_method} \\
            -W ${weak_threshold} -S ${strong_threshold} --add_suffix \\
            -p ${signature_prefix} -i ${params.output_path}/temp/output_tables -s ${params.output_path}/temp/signature_tables -o "./" ${gpu_flag} --seed ${seed}
        """
    }

    run_bootstrap_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files,
        weak_threshold,
        strong_threshold,
        bootstrap_chunks(num_bootstrap_samples, params.bootstrap_chunk_size),
        seed
    )

    emit:
    mutations_tables = run_bootstrap_NNLS.out.mutations_tables
    weights_tables = run_bootstrap_NNLS.out.weights_tables
    stat_infos = run_bootstrap_NNLS.out.stat_infos
    bootstrap_done = run_bootstrap_NNLS.out.stat_infos
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
        label 'gpu_nnls'
        publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true
        
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
        path "./${dataset}/gpu_provenance.json", optional: true, emit: gpu_provenance
        tuple val(dataset), val(mutation_type), emit: dataset_mutation_pairs
        
        when:
        !params.run_only_optimisation
        
        script:
        def optimised_flag = params.optimised ? "-x" : ""
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        python ${workflow.projectDir}/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
            --optimisation_strategy ${params.optimisation_strategy} \\
            -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
            -p ${signature_prefix} -i ${params.output_path}/temp/input_tables -s ${params.output_path}/temp/signature_tables -o "./" ${gpu_flag}

        # Copy residuals and fitted values back to input tables
        cp ${dataset}/output_${dataset}_${mutation_type}_residuals.csv ${params.output_path}/temp/input_tables/${dataset}/
        cp ${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv ${params.output_path}/temp/input_tables/${dataset}/
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
    seed                         // run-wide RNG seed (resolved in the main script)
    
    main:
    // Bootstrap process using penalty files directly. Each task runs one chunk of
    // bootstrap iterations in-process, writing indices start_index .. start_index +
    // n_per_task - 1. Chunks come from bootstrap_chunks() (see below).
    process run_final_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${chunk[0]}"
        label 'gpu_nnls'
        publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true

        input:
        tuple val(dataset), val(mutation_type)
        path weak_penalty
        path strong_penalty
        path input_files
        path signature_files
        each chunk                   // [start_index, n_per_task]
        val seed

        output:
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_mutations_table.csv", emit: mutations_tables
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_stat_info.csv", emit: stat_infos
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_weights_table.csv", emit: weights_tables

        when:
        !params.run_only_optimisation

        script:
        def (start_index, n_per_task) = chunk
        def optimised_flag = params.optimised ? "-x" : ""
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        # Run this task's slice of bootstrap iterations in-process; run_NNLS.py writes
        # one indexed output set per iteration into ${dataset}/bootstrap_output/.
        python ${workflow.projectDir}/bin/run_NNLS.py -B --n_bootstrap ${n_per_task} --bootstrap_start_index ${start_index} \\
            -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
            --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \\
            -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
            -p ${signature_prefix} -i ${params.output_path}/temp/input_tables -s ${params.output_path}/temp/signature_tables -o "./" ${gpu_flag} --seed ${seed}
        """
    }

    run_final_bootstrap_NNLS(
        penalties_for_attribution,
        optimal_weak_penalty_files,
        optimal_strong_penalty_files,
        input_files,
        signature_files,
        bootstrap_chunks(num_bootstrap_samples, params.bootstrap_chunk_size),
        seed
    )

    emit:
    mutations_tables = run_final_bootstrap_NNLS.out.mutations_tables
    stat_infos = run_final_bootstrap_NNLS.out.stat_infos
    weights_tables = run_final_bootstrap_NNLS.out.weights_tables
    bootstrap_done = run_final_bootstrap_NNLS.out.stat_infos
}