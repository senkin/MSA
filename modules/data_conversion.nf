// modules/data_conversion.nf

// parameters for helper flags
// params.strands_flag = ''
// params.nontranscribed_flag = ''
// params.error_flag = ''
// params.COSMIC_flag = ''
def signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow convert_data_workflow {
    take:
    dataset
    input_path
    convert_type // 'matrix_generator_matrices', 'extractor_matrices', or 'signature_tables'

    main:
    // Define and invoke the combined process
    process convert_data {
        tag "${convert_type}/${dataset}"
        // publishDir "${params.input_tables}", mode: 'move', overwrite: true
        // Dynamic publishDir based on convert_type
        publishDir (
            convert_type == 'signature_tables' ? "${params.signature_tables}" : "${params.input_tables}",
            mode: 'move', 
            overwrite: true
        )

        input:
        val dataset
        path input_path
        val convert_type

        output:
        path '*.csv', emit: signatures_for_spectra, optional: true
        path '*.csv', emit: signatures_for_unoptimised_NNLS, optional: true
        path '*/*.csv', emit: converted_SP_to_MSA_for_spectra, optional: true
        path '*/*.csv', emit: converted_SP_to_MSA_for_unoptimised_NNLS, optional: true

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
                                                     -n ${signature_prefix} ${params.COSMIC_flag} \
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

    // Emit the converted output
    emit:
    signatures_for_spectra = convert_data.out.signatures_for_spectra
    signatures_for_unoptimised_NNLS = convert_data.out.signatures_for_unoptimised_NNLS
    converted_SP_to_MSA_for_spectra = convert_data.out.converted_SP_to_MSA_for_spectra
    converted_SP_to_MSA_for_unoptimised_NNLS = convert_data.out.converted_SP_to_MSA_for_unoptimised_NNLS
}