// def signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow OPTIMAL_PENALTIES_workflow {
    take:
    bootstrap_tables        // all bootstrap tables from previous step
    dataset_mutation_pairs  // tuple of (dataset, mutation_type)
    weak_thresholds         // list of weak thresholds
    strong_thresholds       // list of strong thresholds
    
    main:
    // Process to calculate optimal penalties
    process calculate_optimal_penalties {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.optimisation_NNLS_output_path}"
        
        input:
        path '*.csv'  // This mimics the DSL1 file ('*.csv') pattern
        tuple val(dataset), val(mutation_type)
        val weak_thresholds_list
        val strong_thresholds_list
        
        output:
        tuple val(dataset), val(mutation_type), emit: penalties_for_optimisation_plotting
        tuple val(dataset), val(mutation_type), emit: penalties_for_central_NNLS_attribution  
        tuple val(dataset), val(mutation_type), emit: penalties_for_bootstrap_NNLS_attribution
        path "*/*/*.csv", emit: penalty_analysis_csv
        path "*/*/*.json", emit: penalty_analysis_json
        path "SIM_${dataset}/${mutation_type}/optimal_weak_penalty", emit: optimal_weak_penalty
        path "SIM_${dataset}/${mutation_type}/optimal_strong_penalty", emit: optimal_strong_penalty
        
        script:
        // Handle optional flags
        def no_CI_for_penalties_flag = params.no_CI_for_penalties ? "--no_CI" : ""
        def calculate_penalty_on_average_flag = params.calculate_penalty_on_average ? "--average" : ""
        def prioritised_signatures_flag = params.signatures_to_prioritise ? "--signatures_to_prioritise " + params.signatures_to_prioritise.join(' ') : ''
        
        """
        python $baseDir/bin/calculate_optimal_penalties.py -d SIM_${dataset} -t ${mutation_type} \\
            -I ${params.temp_path}/output_tables/SIM_${dataset} \\
            -i ${params.optimisation_NNLS_output_path} -o "./" \\
            -c ${params.SBS_context} \\
            ${no_CI_for_penalties_flag} \\
            ${calculate_penalty_on_average_flag} \\
            ${prioritised_signatures_flag} \\
            -M ${params.metric_to_prioritise} \\
            -T ${params.metric_threshold} \\
            -W ${weak_thresholds_list.join(' ')} \\
            -S ${strong_thresholds_list.join(' ')} \\
            --signature_path ${params.temp_path}/signature_tables \\
            -p ${params.signature_prefix}
        """
    }
    
    // Collect all bootstrap tables
    bootstrap_tables
        .collect()
        .set { all_bootstrap_tables }
    
    // Since we're now within the mutation_type loop, we know exactly what 
    // dataset and mutation_type we're dealing with - no need for unique()
    calculate_optimal_penalties(
        all_bootstrap_tables,
        dataset_mutation_pairs,  // This is already the specific pair for this mutation type
        Channel.value(weak_thresholds),
        Channel.value(strong_thresholds)
    )
    
    emit:
    penalties_for_optimisation_plotting = calculate_optimal_penalties.out.penalties_for_optimisation_plotting
    penalties_for_central_NNLS_attribution = calculate_optimal_penalties.out.penalties_for_central_NNLS_attribution
    penalties_for_bootstrap_NNLS_attribution = calculate_optimal_penalties.out.penalties_for_bootstrap_NNLS_attribution
    penalty_analysis_csv = calculate_optimal_penalties.out.penalty_analysis_csv
    penalty_analysis_json = calculate_optimal_penalties.out.penalty_analysis_json
    optimal_weak_penalty = calculate_optimal_penalties.out.optimal_weak_penalty
    optimal_strong_penalty = calculate_optimal_penalties.out.optimal_strong_penalty
}