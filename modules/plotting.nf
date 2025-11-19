// modules/plotting.nf

def abs_flag = (params.use_absolute_attributions) ? "-a" : ''
def strands_flag = (params.show_strands) ? "-b" : ''
def nontranscribed_flag = (params.show_nontranscribed_region) ? "-n" : ''
def error_flag = (params.show_poisson_errors) ? "-e" : ''
def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow plot_spectra_workflow {
    take:
    dataset
    mutation_type
    inputs
    plot_type // 'mutation_spectra' or 'signatures'

    main:
    process plot_spectra {
        tag "${mutation_type}/${dataset}/${plot_type}"
        publishDir "${params.plots_output_path}", mode: 'move', overwrite: true

        input:
        val dataset
        val mutation_type
        path inputs
        val plot_type

        output:
        path '*/*/*.pdf', optional: true
        path '*/*/*/*.pdf', optional: true
        path '*/*/*/*/*.pdf', optional: true
        path '*/*/*/*/*/*.pdf', optional: true

        script:
        if (plot_type == 'mutation_spectra') {
            """
            python ${workflow.projectDir}/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                            -i ${inputs} ${strands_flag} ${nontranscribed_flag} -o "./"
            python ${workflow.projectDir}/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                            -r -i ${inputs} ${strands_flag} ${nontranscribed_flag} -o "./"
            """
        } else if (plot_type == 'signatures') {
            """
            python ${workflow.projectDir}/bin/plot_mutation_spectra.py -S -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \
                                                        -p ${signature_prefix} -s ${inputs} \
                                                        -r ${strands_flag} ${nontranscribed_flag} -o "./"
            """
        } else {
            """
            echo "Skipping ${plot_type} for ${dataset} and ${mutation_type}"
            """
        }
    }
    // Invoke the process
    plot_spectra(dataset, mutation_type, inputs, plot_type)
}

// Bootstrap attributions plotting workflow
workflow BOOTSTRAP_ATTRIBUTIONS_PLOTS_workflow {
    take:
    attribution_for_plots      // tuple of (dataset, mutation_type)
    bootstrap_attributions     // attributions_per_sample files
    
    main:
    process plot_bootstrap_attributions {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.plots_output_path}", mode: 'move', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        path bootstrap_attributions_files
        
        output:
        path "*/*/bootstrap_plots/*.pdf", optional: true, emit: bootstrap_plots_level1
        path "*/*/bootstrap_plots/*/*.pdf", optional: true, emit: bootstrap_plots_level2
        path "*/*/bootstrap_plots/*/*/*.pdf", optional: true, emit: bootstrap_plots_level3
        
        when:
        params.plot_bootstrap_attributions
        
        script:
        """
        python ${workflow.projectDir}/bin/plot_bootstrap_attributions.py -d ${dataset} -t ${mutation_type} -p ${signature_prefix} ${abs_flag} \\
            -c ${params.SBS_context} -S ${params.temp_path}/signature_tables -I ${params.temp_path}/input_tables \\
            -i ${params.tables_output_path} -o "./" -n ${params.number_of_bootstrapped_samples}
        """
    }
    
    // Collect all bootstrap attribution files
    bootstrap_attributions
        .collect()
        .set { all_bootstrap_attributions }
    
    plot_bootstrap_attributions(
        attribution_for_plots,
        all_bootstrap_attributions
    )
    
    emit:
    bootstrap_plots_level1 = plot_bootstrap_attributions.out.bootstrap_plots_level1
    bootstrap_plots_level2 = plot_bootstrap_attributions.out.bootstrap_plots_level2
    bootstrap_plots_level3 = plot_bootstrap_attributions.out.bootstrap_plots_level3
}

// Metrics plotting workflow
workflow METRICS_PLOTS_workflow {
    take:
    attribution_for_plots      // tuple of (dataset, mutation_type)
    signature_prevalences      // signature_prevalences files
    
    main:
    process plot_metrics {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.plots_output_path}", mode: 'move', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        path prevalences_files
        
        output:
        path "*/*/metrics_plots/*.pdf", optional: true, emit: metrics_plots_level1
        path "*/*/metrics_plots/*/*.pdf", optional: true, emit: metrics_plots_level2
        path "*/*/metrics_plots/*/*/*.pdf", optional: true, emit: metrics_plots_level3
        
        when:
        params.plot_metrics
        
        script:
        """
        python ${workflow.projectDir}/bin/plot_metrics.py -d ${dataset} -t ${mutation_type} \\
            -l ${params.metric_threshold} -i ${params.tables_output_path} -o "./"
        """
    }
    
    // Collect all prevalence files
    signature_prevalences
        .collect()
        .set { all_prevalences }
    
    plot_metrics(
        attribution_for_plots,
        all_prevalences
    )
    
    emit:
    metrics_plots_level1 = plot_metrics.out.metrics_plots_level1
    metrics_plots_level2 = plot_metrics.out.metrics_plots_level2
    metrics_plots_level3 = plot_metrics.out.metrics_plots_level3
}

// Fitted spectra plotting workflow
workflow FITTED_SPECTRA_PLOTS_workflow {
    take:
    attribution_for_plots      // tuple of (dataset, mutation_type)
    
    main:
    process plot_fitted_spectra {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.plots_output_path}", mode: 'move', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        
        output:
        path "*/*/*.pdf", optional: true, emit: fitted_plots_level1
        path "*/*/*/*.pdf", optional: true, emit: fitted_plots_level2
        path "*/*/*/*/*.pdf", optional: true, emit: fitted_plots_level3
        path "*/*/*/*/*/*.pdf", optional: true, emit: fitted_plots_level4
        
        when:
        params.plot_fitted_spectra
        
        script:
        """
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -r -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -r -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        """
    }
    
    plot_fitted_spectra(attribution_for_plots)
    
    emit:
    fitted_plots_level1 = plot_fitted_spectra.out.fitted_plots_level1
    fitted_plots_level2 = plot_fitted_spectra.out.fitted_plots_level2
    fitted_plots_level3 = plot_fitted_spectra.out.fitted_plots_level3
    fitted_plots_level4 = plot_fitted_spectra.out.fitted_plots_level4
}

// Residuals plotting workflow
workflow RESIDUALS_PLOTS_workflow {
    take:
    attribution_for_plots      // tuple of (dataset, mutation_type)
    
    main:
    process plot_residuals {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.plots_output_path}", mode: 'move', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        
        output:
        path "*/*/*.pdf", optional: true, emit: residuals_plots_level1
        path "*/*/*/*.pdf", optional: true, emit: residuals_plots_level2
        path "*/*/*/*/*.pdf", optional: true, emit: residuals_plots_level3
        path "*/*/*/*/*/*.pdf", optional: true, emit: residuals_plots_level4
        
        when:
        params.plot_residuals
        
        script:
        """
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -H -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        python ${workflow.projectDir}/bin/plot_mutation_spectra.py -C -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
            -i ${params.tables_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
        """
    }
    
    plot_residuals(attribution_for_plots)
    
    emit:
    residuals_plots_level1 = plot_residuals.out.residuals_plots_level1
    residuals_plots_level2 = plot_residuals.out.residuals_plots_level2
    residuals_plots_level3 = plot_residuals.out.residuals_plots_level3
    residuals_plots_level4 = plot_residuals.out.residuals_plots_level4
}

// Combined plotting workflow for convenience
workflow ALL_FINAL_PLOTS_workflow {
    take:
    attribution_for_plots      // tuple of (dataset, mutation_type)
    bootstrap_attributions     // attributions_per_sample files
    signature_prevalences      // signature_prevalences files
    
    main:
    // Plot bootstrap attributions
    BOOTSTRAP_ATTRIBUTIONS_PLOTS_workflow(
        attribution_for_plots,
        bootstrap_attributions
    )
    
    // Plot metrics
    METRICS_PLOTS_workflow(
        attribution_for_plots,
        signature_prevalences
    )
    
    // Plot fitted spectra
    FITTED_SPECTRA_PLOTS_workflow(
        attribution_for_plots
    )
    
    // Plot residuals
    RESIDUALS_PLOTS_workflow(
        attribution_for_plots
    )
    
    emit:
    bootstrap_plots = BOOTSTRAP_ATTRIBUTIONS_PLOTS_workflow.out.bootstrap_plots_level1
    metrics_plots = METRICS_PLOTS_workflow.out.metrics_plots_level1
    fitted_plots = FITTED_SPECTRA_PLOTS_workflow.out.fitted_plots_level1
    residuals_plots = RESIDUALS_PLOTS_workflow.out.residuals_plots_level1
}