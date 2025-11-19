// modules/optimization_plots.nf

workflow OPTIMISATION_PLOTS_workflow {
    take:
    penalties_for_plotting  // tuple of (dataset, mutation_type) from optimal penalties
    weak_thresholds         // list of weak thresholds
    strong_thresholds       // list of strong thresholds
    
    main:
    // Process to plot optimisation plots
    process plot_optimisation_plots {
        tag "${mutation_type}/${dataset}"
        publishDir "${params.optimisation_plots_output_path}", mode: 'move', overwrite: true
        
        input:
        tuple val(dataset), val(mutation_type)
        val weak_thresholds_list
        val strong_thresholds_list
        
        output:
        path "*/*.pdf", optional: true, emit: plots_level1
        path "*/*/*.pdf", optional: true, emit: plots_level2
        path "*/*/*/*.pdf", optional: true, emit: plots_level3
        
        when:
        params.plot_optimisation_plots
        
        script:
        def signature_prefix = (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + "_conv" : params.signature_prefix
        """
        python ${workflow.projectDir}/bin/plot_metric_heatmaps.py -d SIM_${dataset} -t ${mutation_type} \\
            -i ${params.optimisation_NNLS_output_path} -o "./" \\
            -c ${params.SBS_context} \\
            -l ${params.metric_threshold} \\
            -W ${weak_thresholds_list.join(' ')} \\
            -S ${strong_thresholds_list.join(' ')} \\
            -T ${params.signature_attribution_thresholds.join(' ')} \\
            --signature_path ${params.temp_path}/signature_tables \\
            -p ${signature_prefix}
        """
    }
    
    // Run plotting for each dataset/mutation_type combination
    plot_optimisation_plots(
        penalties_for_plotting,
        Channel.value(weak_thresholds),
        Channel.value(strong_thresholds)
    )
    
    emit:
    plots_level1 = plot_optimisation_plots.out.plots_level1
    plots_level2 = plot_optimisation_plots.out.plots_level2
    plots_level3 = plot_optimisation_plots.out.plots_level3
}