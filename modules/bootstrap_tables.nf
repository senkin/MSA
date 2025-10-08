// modules/bootstrap_tables.nf

// def signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow BOOTSTRAP_TABLES_workflow {
    take:
    dataset_mutation_pairs  // tuple of (dataset, mutation_type)
    bootstrap_outputs       // all bootstrap outputs from previous step
    weak_thresholds         // list of weak thresholds
    strong_thresholds       // list of strong thresholds
    num_bootstrap_samples   // number of bootstrap samples
    
    main:
    // Process to make bootstrap tables for each threshold combination
    process make_optimisation_bootstrap_tables {
        tag "${mutation_type}/${dataset}/${weak_threshold}/${strong_threshold}"
        publishDir "${params.optimisation_NNLS_output_path}"
        
        input:
        tuple val(dataset), val(mutation_type)
        val bootstrap_outputs_ready  // signal that all bootstrap outputs are complete
        each weak_threshold
        each strong_threshold
        val num_bootstrap_samples
        
        output:
        path "*/*.csv", emit: bootstrap_tables
        path "*/truth_studies/*.csv", optional: true, emit: truth_studies_csv
        path "*/truth_studies/*.json", optional: true, emit: truth_studies_json
        
        script:
        """
        python $baseDir/bin/make_bootstrap_tables.py -d SIM_${dataset} -t ${mutation_type} -p ${params.signature_prefix} \\
            --suffix ${weak_threshold}_${strong_threshold} -l ${params.confidence_level} \\
            -c ${params.SBS_context} -S ${params.temp_path}/signature_tables \\
            -T ${params.signature_attribution_thresholds.join(' ')} \\
            -i ${params.optimisation_NNLS_output_path} -o "./" -n ${num_bootstrap_samples}
        """
    }
    
    // Wait for all bootstrap outputs to complete
    bootstrap_outputs
        .collect()
        .set { bootstrap_ready_signal }
    
    // Run make_optimisation_bootstrap_tables for each combination
    // The 'each' directive will create all combinations automatically
    make_optimisation_bootstrap_tables(
        dataset_mutation_pairs,
        bootstrap_ready_signal,
        Channel.from(weak_thresholds),    // each weak_threshold
        Channel.from(strong_thresholds),   // each strong_threshold
        num_bootstrap_samples
    )
    
    emit:
    bootstrap_tables = make_optimisation_bootstrap_tables.out.bootstrap_tables
    truth_studies_csv = make_optimisation_bootstrap_tables.out.truth_studies_csv
    truth_studies_json = make_optimisation_bootstrap_tables.out.truth_studies_json
}


// Final bootstrap tables workflow (no penalties needed)
workflow FINAL_BOOTSTRAP_TABLES_workflow {
    take:
    attribution_for_tables  // tuple of (dataset, mutation_type) from final NNLS
    bootstrap_outputs       // bootstrap outputs from final bootstrap NNLS
    num_bootstrap_samples   // number of bootstrap samples
    suffix                  // suffix for output files (e.g., "_abs_mutations" or "_weights")
    
    main:
    // Process to make final bootstrap tables
    process make_bootstrap_tables {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.tables_output_path}", mode: 'copy', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        val bootstrap_outputs_ready
        val num_bootstrap_samples
        val suffix
        
        output:
        path "./${dataset}/CIs_${dataset}_${mutation_type}_bootstrap_output_${suffix}.csv", emit: confidence_intervals
        path "./${dataset}/signatures_prevalences_${dataset}_${mutation_type}.csv", emit: signature_prevalences
        path "./${dataset}/attributions_per_sample_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json", emit: attributions_per_sample
        path "./${dataset}/attributions_per_signature_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json", emit: attributions_per_signature
        path "./${dataset}/stat_metrics_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json", emit: stat_metrics
        path "./${dataset}/pruned_attribution_${dataset}_${mutation_type}_abs_mutations.csv", emit: pruned_attribution
        path "*/truth_studies/*.csv", optional: true, emit: truth_studies_csv
        path "*/truth_studies/*.json", optional: true, emit: truth_studies_json
        
        script:
        def abs_flag = (params.use_absolute_attributions) ? "-a" : ''
        """
        mkdir -p ${params.tables_output_path}/${dataset}
        
        # Copy simulation files if this is a simulated dataset
        if [[ ${dataset} == *"SIM"* ]]; then
            if [[ ${mutation_type} == "SBS" ]]; then
                cp ${workflow.projectDir}/input_mutation_tables/${dataset}/WGS_${dataset}.${params.SBS_context}.weights.csv ${params.tables_output_path}/${dataset}/
            elif [[ ${mutation_type} == "DBS" ]]; then
                cp ${workflow.projectDir}/input_mutation_tables/${dataset}/WGS_${dataset}.dinucs.weights.csv ${params.tables_output_path}/${dataset}/
            elif [[ ${mutation_type} == "ID" ]]; then
                cp ${workflow.projectDir}/input_mutation_tables/${dataset}/WGS_${dataset}.indels.weights.csv ${params.tables_output_path}/${dataset}/
            elif [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]]; then
                cp ${workflow.projectDir}/input_mutation_tables/${dataset}/WGS_${dataset}.${mutation_type}.weights.csv ${params.tables_output_path}/${dataset}/
            fi
        fi
        
        python $baseDir/bin/make_bootstrap_tables.py -d ${dataset} -t ${mutation_type} -p ${params.signature_prefix} ${abs_flag} \\
            -c ${params.SBS_context} -S ${params.temp_path}/signature_tables -l ${params.confidence_level} \\
            -T ${params.signature_attribution_thresholds.join(' ')} \\
            -i ${params.tables_output_path} -o "./" -n ${num_bootstrap_samples}
        """
    }
    
    // Wait for all bootstrap outputs to complete
    bootstrap_outputs
        .collect()
        .set { bootstrap_ready_signal }
    
    // Run final bootstrap tables generation
    make_bootstrap_tables(
        attribution_for_tables,
        bootstrap_ready_signal,
        num_bootstrap_samples,
        suffix
    )
    
    emit:
    confidence_intervals = make_bootstrap_tables.out.confidence_intervals
    signature_prevalences = make_bootstrap_tables.out.signature_prevalences
    attributions_per_sample = make_bootstrap_tables.out.attributions_per_sample
    attributions_per_signature = make_bootstrap_tables.out.attributions_per_signature
    stat_metrics = make_bootstrap_tables.out.stat_metrics
    pruned_attribution = make_bootstrap_tables.out.pruned_attribution
    truth_studies_csv = make_bootstrap_tables.out.truth_studies_csv
    truth_studies_json = make_bootstrap_tables.out.truth_studies_json
}