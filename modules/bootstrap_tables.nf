// modules/bootstrap_tables.nf

workflow BOOTSTRAP_TABLES_workflow {
    take:
    dataset_mutation_pairs  // tuple of (dataset, mutation_type)
    bootstrap_outputs       // all bootstrap outputs from previous step
    weak_thresholds         // channel: list of weak threshold labels
    strong_thresholds       // channel: list of strong threshold labels
    num_bootstrap_samples   // number of bootstrap samples
    
    main:
    // Process to make bootstrap tables for each threshold combination
    process make_optimisation_bootstrap_tables {
        tag "${mutation_type}/${dataset}/${weak_threshold}/${strong_threshold}"
        publishDir "${params.output_path}/outputs_optimisation"
        
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
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        python ${workflow.projectDir}/bin/make_bootstrap_tables.py -d SIM_${dataset} -t ${mutation_type} -p ${signature_prefix} \\
            --suffix ${weak_threshold}_${strong_threshold} -l ${params.confidence_level} \\
            -c ${params.SBS_context} -S ${params.output_path}/temp/signature_tables \\
            -T ${params.signature_attribution_thresholds.join(' ')} \\
            -i ${params.output_path}/outputs_optimisation -o "./" -n ${num_bootstrap_samples}
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
        weak_thresholds,    // each weak_threshold
        strong_thresholds,  // each strong_threshold
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
    
    main:
    // Process to make final bootstrap tables
    process make_bootstrap_tables {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        val bootstrap_outputs_ready
        val num_bootstrap_samples
        
        output:
        path "${dataset}/CIs_${dataset}_${mutation_type}_bootstrap_output_*.csv", emit: confidence_intervals
        path "${dataset}/signatures_prevalences_${dataset}_${mutation_type}.csv", emit: signature_prevalences
        path "${dataset}/attributions_per_sample_${dataset}_${mutation_type}_bootstrap_output_*.json", emit: attributions_per_sample
        path "${dataset}/attributions_per_signature_${dataset}_${mutation_type}_bootstrap_output_*.json", emit: attributions_per_signature
        path "${dataset}/stat_metrics_${dataset}_${mutation_type}_bootstrap_output_*.json", emit: stat_metrics
        path "${dataset}/pruned_attribution_${dataset}_${mutation_type}_abs_mutations.csv", emit: pruned_attribution
        path "*/truth_studies/*.csv", optional: true, emit: truth_studies_csv
        path "*/truth_studies/*.json", optional: true, emit: truth_studies_json
        
        script:
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        def table = [SBS: "${params.SBS_context}", DBS: 'dinucs', ID: 'indels'].get(mutation_type, mutation_type)
        """
        mkdir -p ${params.output_path}/output_tables/${dataset}

        # Simulated datasets come with their truth table (signature weights): when the staged
        # inputs have one, make_bootstrap_tables.py scores the attribution against it
        truth=${params.output_path}/temp/input_tables/${dataset}/WGS_${dataset}.${table}.weights.csv
        if [[ -f "\$truth" ]]; then
            cp "\$truth" ${params.output_path}/output_tables/${dataset}/
        fi

        python ${workflow.projectDir}/bin/make_bootstrap_tables.py -d ${dataset} -t ${mutation_type} -p ${signature_prefix} \\
            -c ${params.SBS_context} -S ${params.output_path}/temp/signature_tables -l ${params.confidence_level} \\
            -T ${params.signature_attribution_thresholds.join(' ')} \\
            -i ${params.output_path}/output_tables -o "./" -n ${num_bootstrap_samples}
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