// modules/plot_spectra.nf

// parameters for helper flags
params.strands_flag = ''
params.nontranscribed_flag = ''
params.signature_prefix = ''

workflow plot_spectra_workflow {
    take:
    dataset
    mutation_type
    inputs
    plot_type // 'mutation_spectra' or 'signatures'

    main:
    process plot_spectra {
        tag "${mutation_type}/${dataset}/${plot_type}"
        publishDir "${params.plots_output_path}", mode: 'move'

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
            python $baseDir/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                            -i ${params.input_tables} ${params.strands_flag} ${params.nontranscribed_flag} -o "./"
            python $baseDir/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                            -r -i ${params.input_tables} ${params.strands_flag} ${params.nontranscribed_flag} -o "./"
            """
        } else if (plot_type == 'signatures') {
            """
            python $baseDir/bin/plot_mutation_spectra.py -S -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \
                                                        -p ${params.signature_prefix} -s ${params.signature_tables} \
                                                        -r ${params.strands_flag} ${params.nontranscribed_flag} -o "./"
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