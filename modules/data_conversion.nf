// modules/data_conversion.nf

// Convert SigProfiler outputs, or a single user-provided file, to MSA tables: input mutation
// tables (<dataset>/WGS_<dataset>.<table>.csv) or signature tables (<prefix>_conv_*.csv).
// convert_type: 'matrix_generator_matrices', 'extractor_matrices', 'signature_tables',
// 'specific_mutation_files' or 'specific_signature_files'.
process convert_data {
    tag "${convert_type}/${dataset}"
    // converted tables are intermediate files: kept in <output_path>/temp with cleanup_temp = false
    // (a closure, since the path uses an input: as a plain string it fails with the strict parser)
    publishDir { "${params.output_path}/temp/${convert_type.contains('signature') ? 'signature_tables' : 'input_tables'}" },
        mode: 'copy', overwrite: true, enabled: !params.cleanup_temp.toBoolean()

    input:
    val dataset
    path source, stageAs: 'source/*'                   // SigProfiler output folder or a single file
    val convert_type
    val mutation_types_list
    path reference_signatures, stageAs: 'reference_signatures/*'   // index templates

    output:
    path '*.csv', emit: signature_files, optional: true
    path "${dataset}/*.csv", emit: input_files, optional: true

    script:
    def base_cmd = "convert_SP_to_MSA.py"
    def mutation_types_arg = "-t ${mutation_types_list.join(' ')}"
    def signature_tables_arg = "-s reference_signatures"
    def output_arg = "-o ./"
    def COSMIC_arg = params.COSMIC_signatures.toBoolean() ? "--COSMIC" : ""
    def signature_prefix_arg = "${params.signature_prefix}_conv"

    if (convert_type == 'matrix_generator_matrices') {
        """
        ${base_cmd} -d ${dataset} ${mutation_types_arg} \\
            -i ${source} ${signature_tables_arg} ${output_arg}
        """
    } else if (convert_type == 'extractor_matrices') {
        """
        ${base_cmd} -E -d ${dataset} ${mutation_types_arg} \\
            -i ${source} ${signature_tables_arg} ${output_arg}
        """
    } else if (convert_type == 'signature_tables') {
        """
        ${base_cmd} -S ${mutation_types_arg} \\
            -n ${signature_prefix_arg} ${COSMIC_arg} \\
            -i ${source} ${signature_tables_arg} ${output_arg}
        """
    } else if (convert_type == 'specific_mutation_files') {
        """
        ${base_cmd} -d ${dataset} ${mutation_types_arg} -c ${params.SBS_context} \\
            -I ${source} ${signature_tables_arg} ${output_arg}
        """
    } else if (convert_type == 'specific_signature_files') {
        """
        ${base_cmd} -S ${mutation_types_arg} -c ${params.SBS_context} \\
            -n ${signature_prefix_arg} ${COSMIC_arg} \\
            -I ${source} ${signature_tables_arg} ${output_arg}
        """
    } else {
        """
        echo "Invalid convert_type: ${convert_type}"
        echo "Supported types: matrix_generator_matrices, extractor_matrices, signature_tables, specific_mutation_files, specific_signature_files"
        exit 1
        """
    }
}
