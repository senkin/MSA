// modules/bootstrap_tables.nf

include { signature_prefix; table_name } from './common'

// Bootstrap summary tables of one threshold pair of the optimisation (simulated data):
// confidence intervals, prevalences and the truth studies the penalty choice is based on.
// Also writes bundle/<pair>: what the penalty choice and its plots read for this pair (truth
// studies and the central attribution tables), hard-linked where possible, not published.
process make_optimisation_bootstrap_tables {
    tag "${mutation_type}/${dataset}/${suffix}"
    publishDir "${params.output_path}/outputs_optimisation", enabled: !params.cleanup_temp,
        saveAs: { f -> f.startsWith('bundle/') ? null : f }

    input:
    tuple val(dataset), val(mutation_type), val(suffix),
        path("attributions/SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}/*"),
        path("attributions/SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}/bootstrap_output/*")
    path signature_files, stageAs: 'signature_tables/*'
    val num_bootstrap_samples

    output:
    path "SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}/*.csv", emit: tables
    path "SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}/truth_studies/*", emit: truth_studies
    tuple val(dataset), val(mutation_type), path("bundle/SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}", type: 'dir'), emit: bundle

    script:
    def pair = "SIM_${dataset}_${params.SBS_context}_NNLS_${suffix}"
    """
    python ${workflow.projectDir}/bin/make_bootstrap_tables.py -d SIM_${dataset} -t ${mutation_type} -p ${signature_prefix()} \\
        --suffix ${suffix} -l ${params.confidence_level} \\
        -c ${params.SBS_context} -S signature_tables \\
        -T ${params.signature_attribution_thresholds.join(' ')} \\
        -i attributions -o "./" -n ${num_bootstrap_samples}

    mkdir -p bundle/${pair}
    cp -r ${pair}/truth_studies bundle/${pair}/
    for table in weights_table stat_info; do
        central=attributions/${pair}/output_SIM_${dataset}_${mutation_type}_\${table}.csv
        ln -L \$central bundle/${pair}/ 2>/dev/null || cp -L \$central bundle/${pair}/
    done
    """
}

// Bootstrap summary tables of the final attribution: confidence intervals, prevalences,
// pruned attribution and, for simulated data, the truth studies
process make_bootstrap_tables {
    tag "${mutation_type}/${dataset}"
    // the bootstrap distributions per sample/signature (JSON) are intermediate files
    publishDir "${params.output_path}/output_tables", mode: 'copy', overwrite: true,
        saveAs: { f -> (params.cleanup_temp && f ==~ /[^\/]+\/[^\/]+\.json/) ? null : f }

    input:
    tuple val(dataset), val(mutation_type), path("attributions/${dataset}/*"), path("attributions/${dataset}/bootstrap_output/*")
    path "input_tables/${dataset}/*"
    path signature_files, stageAs: 'signature_tables/*'
    val num_bootstrap_samples

    output:
    tuple val(dataset), val(mutation_type),
        path("${dataset}/CIs_${dataset}_${mutation_type}_bootstrap_output_*.csv"),
        path("${dataset}/signatures_prevalences_${dataset}_${mutation_type}.csv"),
        path("${dataset}/attributions_per_sample_${dataset}_${mutation_type}_bootstrap_output_*.json"),
        path("${dataset}/attributions_per_signature_${dataset}_${mutation_type}_bootstrap_output_*.json"),
        path("${dataset}/stat_metrics_${dataset}_${mutation_type}_bootstrap_output_*.json"),
        path("${dataset}/pruned_attribution_${dataset}_${mutation_type}_abs_mutations.csv"),
        path("${dataset}/WGS_${dataset}.*.weights.csv", arity: '0..*'),
        path("${dataset}/truth_studies/*", arity: '0..*'), emit: tables

    script:
    def truth = "WGS_${dataset}.${table_name(mutation_type)}.weights.csv"
    """
    # Simulated datasets come with their truth table (signature weights) among the inputs:
    # make_bootstrap_tables.py then scores the attribution against it
    if [[ -f input_tables/${dataset}/${truth} ]]; then
        ln -s ../../input_tables/${dataset}/${truth} attributions/${dataset}/
    fi

    python ${workflow.projectDir}/bin/make_bootstrap_tables.py -d ${dataset} -t ${mutation_type} -p ${signature_prefix()} \\
        -c ${params.SBS_context} -S signature_tables -l ${params.confidence_level} \\
        -T ${params.signature_attribution_thresholds.join(' ')} \\
        -i attributions -o "./" -n ${num_bootstrap_samples}

    # the truth table is kept with the results (and the plots read it from there)
    if [[ -f input_tables/${dataset}/${truth} ]]; then
        cp -L input_tables/${dataset}/${truth} ${dataset}/
    fi
    """
}
