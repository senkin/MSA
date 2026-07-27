// modules/nnls.nf

// In GPU mode, run_NNLS.py runs the greedy signature optimisation transposed
// across samples (one large masked NNLS batch per greedy step), offloaded to the
// cuML batched NNLS solver via batched_solve_and_score().
def gpu_flag = params.use_GPU ? "--use_gpu --gpu_precision ${params.gpu_precision} --gpu_batch_size ${params.gpu_batch_size}" : ''
// Optional fixed RNG seed for bootstrap resampling (reproducible across CPU/GPU runs).
def seed_flag = params.seed != null ? "--seed ${params.seed}" : ''

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
    
    main:
    // Bootstrap NNLS process. Each task runs `n_per_task` bootstrap iterations
    // in-process via run_NNLS.py --n_bootstrap, writing one indexed output set per
    // iteration (indices start_index .. start_index + n_per_task - 1) directly into
    // the bootstrap_output/ folder. The number of tasks vs iterations-per-task is
    // chosen by the caller depending on params.use_GPU (see below).
    process run_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${weak_threshold}/${strong_threshold}/${start_index}"
        label 'gpu_nnls'
        publishDir "${params.output_path}/outputs_optimisation"

        input:
        val dataset
        val mutation_type
        path input_files
        path signature_files
        val weak_threshold
        val strong_threshold
        val n_per_task
        val start_index

        output:
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_mutations_table.csv", emit: mutations_tables
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_weights_table.csv", emit: weights_tables
        path "SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_*_stat_info.csv", emit: stat_infos

        when:
        !params.run_only_simulations

        script:
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
            -p ${signature_prefix} -i ${params.output_path}/temp/output_tables -s ${params.output_path}/temp/signature_tables -o "./" ${gpu_flag} ${seed_flag}
        """
    }

    // GPU: one task runs all iterations (a single reused context).
    // CPU (default): one task per iteration for maximum parallel width on the scheduler.
    def n_per_task = params.use_GPU ? num_bootstrap_samples : 1
    def start_indices = params.use_GPU ? Channel.value(1) : Channel.of(1..num_bootstrap_samples)

    run_bootstrap_NNLS(
        dataset,
        mutation_type,
        input_files,
        signature_files,
        weak_threshold,
        strong_threshold,
        n_per_task,
        start_indices
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
    
    main:
    // Bootstrap process using penalty files directly. Each task runs `n_per_task`
    // bootstrap iterations in-process, writing indices start_index .. start_index +
    // n_per_task - 1. Tasks-vs-iterations split is chosen by params.use_GPU below.
    process run_final_bootstrap_NNLS {
        tag "${mutation_type}/${dataset}/${start_index}"
        label 'gpu_nnls'
        publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true

        input:
        tuple val(dataset), val(mutation_type)
        path weak_penalty
        path strong_penalty
        path input_files
        path signature_files
        val n_per_task
        each start_index

        output:
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_mutations_table.csv", emit: mutations_tables
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_stat_info.csv", emit: stat_infos
        path "${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_*_weights_table.csv", emit: weights_tables

        when:
        !params.run_only_optimisation

        script:
        def optimised_flag = params.optimised ? "-x" : ""
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        # Run this task's slice of bootstrap iterations in-process; run_NNLS.py writes
        # one indexed output set per iteration into ${dataset}/bootstrap_output/.
        python ${workflow.projectDir}/bin/run_NNLS.py -B --n_bootstrap ${n_per_task} --bootstrap_start_index ${start_index} \\
            -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \\
            --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \\
            -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \\
            -p ${signature_prefix} -i ${params.output_path}/temp/input_tables -s ${params.output_path}/temp/signature_tables -o "./" ${gpu_flag} ${seed_flag}
        """
    }

    // GPU: one task runs all iterations (a single reused context).
    // CPU (default): one task per iteration (each) for maximum parallel width.
    def n_per_task = params.use_GPU ? num_bootstrap_samples : 1
    def start_indices = params.use_GPU ? [1] : (1..num_bootstrap_samples).toList()

    run_final_bootstrap_NNLS(
        penalties_for_attribution,
        optimal_weak_penalty_files,
        optimal_strong_penalty_files,
        input_files,
        signature_files,
        n_per_task,
        start_indices
    )

    emit:
    mutations_tables = run_final_bootstrap_NNLS.out.mutations_tables
    stat_infos = run_final_bootstrap_NNLS.out.stat_infos
    weights_tables = run_final_bootstrap_NNLS.out.weights_tables
    bootstrap_done = run_final_bootstrap_NNLS.out.stat_infos
}