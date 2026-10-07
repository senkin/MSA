// modules/optimisation_plots.nf

include { signature_prefix } from './common'

// Metric heatmaps and confidence probability curves of the threshold scan of a mutation type
process plot_optimisation_plots {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/optimisation_plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path(weak_thresholds_file), path(strong_thresholds_file),
        path("optimisation/SIM_${dataset}/${mutation_type}/*"),   // penalty analysis tables
        path(pairs, stageAs: 'optimisation/*')                    // bundle/<pair> of every threshold pair
    path signature_files, stageAs: 'signature_tables/*'

    output:
    path "*/*.pdf", optional: true, emit: plots_level1
    path "*/*/*.pdf", optional: true, emit: plots_level2
    path "*/*/*/*.pdf", optional: true, emit: plots_level3

    script:
    """
    plot_metric_heatmaps.py -d SIM_${dataset} -t ${mutation_type} \\
        -i optimisation -o "./" \\
        -c ${params.SBS_context} \\
        -l ${params.metric_threshold} \\
        -W \$(cat ${weak_thresholds_file}) \\
        -S \$(cat ${strong_thresholds_file}) \\
        -T ${params.signature_attribution_thresholds.join(' ')} \\
        --signature_path signature_tables \\
        -p ${signature_prefix()}
    """
}
