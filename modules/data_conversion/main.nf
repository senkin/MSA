// ./modules/data_conversion/main.nf

// parameters for helper flags
params.strands_flag = ''
params.nontranscribed_flag = ''
params.error_flag = ''
params.COSMIC_flag = ''
params.signature_prefix = 'sigProfiler'

workflow convert_data_workflow {
    take:
    dataset
    input_path
    convert_type // 'matrix_generator_matrices', 'extractor_matrices', or 'signature_tables'

    main:
    // Define and invoke the combined process
    process convert_data {
        tag "${convert_type}/${dataset}"
        publishDir "${params.input_tables}", mode: 'move', overwrite: true

        input:
        val dataset
        path input_path
        val convert_type

        output:
        path '*.csv', emit: converted_signatures, optional: true
        path '*/*.csv', emit: converted_output, optional: true

        script:
        if (convert_type == 'matrix_generator_matrices') {
            """
            python $baseDir/bin/convert_SP_to_MSA.py -d ${dataset} -t ${params.mutation_types.join(' ')} \
                                                     -i ${input_path} -s ${params.signature_tables} -o "./"
            """
        } else if (convert_type == 'extractor_matrices') {
            """
            python $baseDir/bin/convert_SP_to_MSA.py -E -d ${dataset} -t ${params.mutation_types.join(' ')} \
                                                     -i ${input_path} -s ${params.signature_tables} -o "./"
            """
        } else if (convert_type == 'signature_tables') {
            """
            python $baseDir/bin/convert_SP_to_MSA.py -S -t ${params.mutation_types.join(' ')} \
                                                     -n ${params.signature_prefix} ${params.COSMIC_flag} \
                                                     -i ${input_path} -s ${params.signature_tables} -o "./"
            """
        } else {
            """
            echo "Invalid convert_type: ${convert_type}"
            exit 1
            """
        }
    }

    // Invoke the process
    convert_data(dataset, input_path, convert_type)

    // emit:
    // // Emit the converted output
    // converted_output
}