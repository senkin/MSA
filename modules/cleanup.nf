#!/usr/bin/env nextflow

process cleanup_temp_files {
    when:
    params.cleanup_temp

    input:
    val ready

    output:
    stdout

    script:
    """
    echo "Cleaning up temporary files..."
    rm -rf ${params.output_path}/temp
    rm -rf ${params.output_path}/outputs_optimisation/*_NNLS_*
    rm -rf ${params.output_path}/output_tables/*/*json
    rm -rf ${params.output_path}/output_tables/*/bootstrap_output
    echo "Cleanup complete."
    """
}

workflow cleanup_workflow {
    take:
    completion_signal
    
    main:
    cleanup_temp_files(completion_signal)
    
    emit:
    cleanup_status = cleanup_temp_files.out
}