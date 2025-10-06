// ./modules/stage_default_inputs.nf

process stage_default_signatures_proc {
    tag "stage_default_signatures"
    publishDir "${params.temp_path}", mode: 'copy', overwrite: true
    
    output:
    path 'signature_tables/*.csv', emit: staged_signature_tables
    
    script:
    """
    mkdir -p signature_tables
    if [ -d "${workflow.projectDir}/signature_tables" ]; then
        cp -a "${workflow.projectDir}/signature_tables/"* signature_tables/ || true
        echo "Staged repository signature_tables -> ${params.temp_path}/signature_tables"
    else
        echo "Warning: no repository signature_tables found at ${workflow.projectDir}/signature_tables"
    fi
    """
}

process stage_default_inputs_proc {
    tag "stage_default_inputs/${dataset}"
    publishDir "${params.temp_path}", mode: 'copy', overwrite: true
    
    input:
    val dataset
    
    output:
    path "input_tables/${dataset}/*.csv", emit: staged_input_tables
    
    script:
    """
    mkdir -p input_tables/${dataset}
    if [ -d "${workflow.projectDir}/input_mutation_tables/${dataset}" ]; then
        cp -a "${workflow.projectDir}/input_mutation_tables/${dataset}/"* input_tables/${dataset}/ || true
        echo "Staged repository input_mutation_tables/${dataset} -> ${params.temp_path}/input_tables/${dataset}"
    else
        echo "Warning: no repository test inputs found at ${workflow.projectDir}/input_mutation_tables/${dataset}"
    fi
    """
}

workflow stage_default_signatures_workflow {
    main:
    stage_default_signatures_proc()
    
    emit:
    staged_signature_tables = stage_default_signatures_proc.out.staged_signature_tables
}

workflow stage_default_inputs_workflow {
    take:
    dataset
    
    main:
    stage_default_inputs_proc(dataset)
    
    emit:
    staged_input_tables = stage_default_inputs_proc.out.staged_input_tables
}