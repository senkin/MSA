// modules/plotting.nf
//
// Plots of the input spectra and signatures, and of the final attribution. Each process stages
// what it reads: input_tables/<dataset>/, signature_tables/ and attributions/<dataset>/ (the
// final attribution tables, bootstrap summaries and, for simulated data, the truth studies).

include { signature_prefix } from './common'

def strands_flag() { return params.show_strands.toString().toBoolean() ? "-b" : '' }
def nontranscribed_flag() { return params.show_nontranscribed_region.toString().toBoolean() ? "-n" : '' }
def error_flag() { return params.show_poisson_errors.toString().toBoolean() ? "-e" : '' }

process plot_signatures {
    tag "${mutation_type}/${dataset}/signatures"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type)
    path signature_files, stageAs: 'signature_tables/*'

    output:
    path '*/*/*.pdf', optional: true
    path '*/*/*/*.pdf', optional: true
    path '*/*/*/*/*.pdf', optional: true
    path '*/*/*/*/*/*.pdf', optional: true

    script:
    """
    plot_mutation_spectra.py -S -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \\
                                                -p ${signature_prefix()} -s signature_tables \\
                                                -r ${strands_flag()} ${nontranscribed_flag()} -o "./"
    """
}

process plot_input_spectra {
    tag "${mutation_type}/${dataset}/mutation_spectra"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type)
    path "input_tables/${dataset}/*"

    output:
    path '*/*/*.pdf', optional: true
    path '*/*/*/*.pdf', optional: true
    path '*/*/*/*/*.pdf', optional: true
    path '*/*/*/*/*/*.pdf', optional: true

    script:
    """
    plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
                                                    -i input_tables ${strands_flag()} ${nontranscribed_flag()} -o "./"
    plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
                                                    -r -i input_tables ${strands_flag()} ${nontranscribed_flag()} -o "./"
    """
}

process plot_bootstrap_attributions {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("attributions/${dataset}/*")
    path "input_tables/${dataset}/*"
    path signature_files, stageAs: 'signature_tables/*'

    output:
    path "*/*/bootstrap_plots/*.pdf", optional: true, emit: bootstrap_plots_level1
    path "*/*/bootstrap_plots/*/*.pdf", optional: true, emit: bootstrap_plots_level2
    path "*/*/bootstrap_plots/*/*/*.pdf", optional: true, emit: bootstrap_plots_level3

    script:
    def common_args = "-d ${dataset} -t ${mutation_type} -p ${signature_prefix()} " +
                      "-c ${params.SBS_context} -S signature_tables " +
                      "-I input_tables " +
                      "-i attributions -o ./ -n ${params.number_of_bootstrapped_samples}"
    """
    # Produce both relative (weights) and absolute (mutation count) attribution plots
    plot_bootstrap_attributions.py ${common_args}
    plot_bootstrap_attributions.py ${common_args} -a
    """
}

process plot_metrics {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("attributions/${dataset}/*"), path("attributions/${dataset}/truth_studies/*", arity: '0..*')

    output:
    path "*/*/metrics_plots/*.pdf", optional: true, emit: metrics_plots_level1
    path "*/*/metrics_plots/*/*.pdf", optional: true, emit: metrics_plots_level2
    path "*/*/metrics_plots/*/*/*.pdf", optional: true, emit: metrics_plots_level3

    script:
    """
    plot_metrics.py -d ${dataset} -t ${mutation_type} \\
        -l ${params.metric_threshold} -i attributions -o "./"
    """
}

process plot_fitted_spectra {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("attributions/${dataset}/*")

    output:
    path "*/*/*.pdf", optional: true, emit: fitted_plots_level1
    path "*/*/*/*.pdf", optional: true, emit: fitted_plots_level2
    path "*/*/*/*/*.pdf", optional: true, emit: fitted_plots_level3
    path "*/*/*/*/*/*.pdf", optional: true, emit: fitted_plots_level4

    script:
    """
    plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -r -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -r -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    """
}

process plot_residuals {
    tag "${mutation_type}/${dataset}"
    publishDir "${params.output_path}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(dataset), val(mutation_type), path("attributions/${dataset}/*")

    output:
    path "*/*/*.pdf", optional: true, emit: residuals_plots_level1
    path "*/*/*/*.pdf", optional: true, emit: residuals_plots_level2
    path "*/*/*/*/*.pdf", optional: true, emit: residuals_plots_level3
    path "*/*/*/*/*/*.pdf", optional: true, emit: residuals_plots_level4

    script:
    """
    plot_mutation_spectra.py -H -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    plot_mutation_spectra.py -C -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \\
        -i attributions ${error_flag()} ${strands_flag()} ${nontranscribed_flag()} -o "./"
    """
}
