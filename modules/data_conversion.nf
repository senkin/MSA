// modules/data_conversion.nf

workflow convert_data_workflow {
    take:
    dataset
    input_path
    convert_type // 'matrix_generator_matrices', 'extractor_matrices', 'signature_tables', 'specific_mutation_files', 'specific_signature_files'

    main:
    // Define and invoke the combined process
    process convert_data {
        tag "${convert_type}/${dataset}"
        
        // Use temp_path for all converted outputs
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
        path '*.csv', emit: signature_files, optional: true
        path '**/*.csv', emit: input_files, optional: true

        script:
        // Base command components
        def base_cmd = "python ${workflow.projectDir}/bin/convert_SP_to_MSA.py"
        def mutation_types_arg = "-t ${params.mutation_types.join(' ')}"
        def signature_tables_arg = "-s ${params.signature_tables}"
        def temp_output_arg = "-o ./"  // Output to process work directory first, then publishDir handles the move
        
        if (convert_type == 'matrix_generator_matrices') {
            """
            ${base_cmd} -d ${dataset} ${mutation_types_arg} \\
                -i ${input_path} ${signature_tables_arg} ${temp_output_arg}
            """
        } else if (convert_type == 'extractor_matrices') {
            """
            ${base_cmd} -E -d ${dataset} ${mutation_types_arg} \\
                -i ${input_path} ${signature_tables_arg} ${temp_output_arg}
            """
        } else if (convert_type == 'signature_tables') {
            """
            ${base_cmd} -S ${mutation_types_arg} \\
                -n ${params.signature_prefix} ${params.COSMIC_flag} \\
                -i ${input_path} ${signature_tables_arg} ${temp_output_arg}
            """
        } else if (convert_type == 'specific_mutation_files') {
            // Handle specific mutation files - input_path should be a list of files
            """
            ${base_cmd} -d ${dataset} ${mutation_types_arg} \\
                -I ${input_path} ${signature_tables_arg} ${temp_output_arg}
            """
        } else if (convert_type == 'specific_signature_files') {
            // Handle specific signature files - input_path should be a list of files
            """
            ${base_cmd} -S ${mutation_types_arg} \\
                -n ${params.signature_prefix} ${params.COSMIC_flag} \\
                -I ${input_path} ${signature_tables_arg} ${temp_output_arg}
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
    signature_files = convert_data.out.signature_files
    input_files = convert_data.out.input_files
}