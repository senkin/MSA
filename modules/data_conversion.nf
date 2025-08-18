// modules/data_conversion.nf

// def signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix

workflow convert_data_workflow {
    take:
    dataset
    input_path
    convert_type // 'matrix_generator_matrices', 'extractor_matrices', 'signature_tables', 'specific_mutation_files', 'specific_signature_files'

    main:
    // Define and invoke the combined process
    process convert_data {
        tag "${convert_type}/${dataset}"
        
        // Dynamic publishDir based on convert_type
        publishDir (
            convert_type.contains('signature') ? "${params.temp_path}/signature_tables" : "${params.temp_path}/input_tables",
            mode: 'copy',
            overwrite: true
        )

        input:
        val dataset
        val input_path  // Can be a directory path or list of specific files
        val convert_type

        output:
        path '*.csv', emit: signatures_for_spectra, optional: true
        path '*.csv', emit: signatures_for_unoptimised_NNLS, optional: true
        path '*/*.csv', emit: converted_SP_to_MSA_for_spectra, optional: true
        path '*/*.csv', emit: converted_SP_to_MSA_for_unoptimised_NNLS, optional: true

        script:
        // Base command components
        def base_cmd = "python ${workflow.projectDir}/bin/convert_SP_to_MSA.py"
        def mutation_types_arg = "-t ${params.mutation_types.join(' ')}"
        def signature_tables_arg = "-s ${params.signature_tables}"
        def output_arg = "-o ./"
        
        if (convert_type == 'matrix_generator_matrices') {
            """
            ${base_cmd} -d ${dataset} ${mutation_types_arg} \\
                -i ${input_path} ${signature_tables_arg} ${output_arg}
            """
        } else if (convert_type == 'extractor_matrices') {
            """
            ${base_cmd} -E -d ${dataset} ${mutation_types_arg} \\
                -i ${input_path} ${signature_tables_arg} ${output_arg}
            """
        } else if (convert_type == 'signature_tables') {
            """
            ${base_cmd} -S ${mutation_types_arg} \\
                -n ${params.signature_prefix} ${params.COSMIC_flag} \\
                -i ${input_path} ${signature_tables_arg} ${output_arg}
            """
        } else if (convert_type == 'specific_mutation_files') {
            // Handle specific mutation files - input_path should be a list of files
            """
            ${base_cmd} -d ${dataset} ${mutation_types_arg} \\
                -I ${input_path} ${signature_tables_arg} ${output_arg}
            """
        } else if (convert_type == 'specific_signature_files') {
            // Handle specific signature files - input_path should be a list of files
            """
            ${base_cmd} -S ${mutation_types_arg} \\
                -n ${params.signature_prefix} ${params.COSMIC_flag} \\
                -I ${input_path} ${signature_tables_arg} ${output_arg}
            """
        } else {
            """
            echo "Invalid convert_type: ${convert_type}"
            echo "Supported types: matrix_generator_matrices, extractor_matrices, signature_tables, specific_mutation_files, specific_signature_files"
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

// Additional workflow for handling specific files
workflow convert_specific_files_workflow {
    take:
    dataset
    input_files_list  // List of specific files to convert
    file_type        // 'mutation' or 'signature'

    main:
    if (file_type == 'mutation') {
        convert_data_workflow(dataset, input_files_list, 'specific_mutation_files')
    } else if (file_type == 'signature') {
        convert_data_workflow(dataset, input_files_list, 'specific_signature_files')
    } else {
        error "Invalid file_type: ${file_type}. Must be 'mutation' or 'signature'"
    }

    emit:
    signatures_for_spectra = convert_data_workflow.out.signatures_for_spectra
    signatures_for_unoptimised_NNLS = convert_data_workflow.out.signatures_for_unoptimised_NNLS
    converted_SP_to_MSA_for_spectra = convert_data_workflow.out.converted_SP_to_MSA_for_spectra
    converted_SP_to_MSA_for_unoptimised_NNLS = convert_data_workflow.out.converted_SP_to_MSA_for_unoptimised_NNLS
}