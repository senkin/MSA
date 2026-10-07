// modules/common.nf: helpers shared by the processes and workflows

// Signature tables converted from SigProfiler outputs (or --signatures_file) carry a '_conv' prefix
def signature_prefix() {
    return (params.SP_extractor_output_path || params.signatures_file) ? params.signature_prefix + '_conv' : params.signature_prefix
}

// In GPU mode, run_NNLS.py runs the greedy signature optimisation transposed across samples
// (one large masked NNLS batch per greedy step), offloaded to the cuML batched NNLS solver
def gpu_flag() {
    return params.use_GPU ? "--use_gpu --gpu_precision ${params.gpu_precision} --gpu_batch_size ${params.gpu_batch_size}" : ''
}

// Name of a mutation type's tables: WGS_<dataset>.<table>.csv (and .weights.csv for simulated data)
def table_name(mutation_type) {
    return [SBS: "${params.SBS_context}", DBS: 'dinucs', ID: 'indels'].get(mutation_type, mutation_type)
}

// Split bootstrap iterations 1..n into [start_index, n_iterations] chunks, one chunk per task.
// run_NNLS.py seeds every iteration from (seed, index), so results do not depend on the split.
// chunk_size <= 0 keeps the mode default: one iteration per task on CPU (maximum scheduler
// width), all iterations in a single task on GPU (one reused CUDA context). No chunks (no
// tasks) for n < 1, which the main script only allows where the bootstrap is not run.
def bootstrap_chunks(n, chunk_size) {
    def total = n as int
    if (total < 1) {
        return []
    }
    def size = ((chunk_size ?: 0) as int) > 0 ? (chunk_size as int) : (params.use_GPU ? total : 1)
    return (1..total).step(size).collect { start -> [start, Math.min(size, total - start + 1)] }
}

// Threshold list for run_NNLS.py -W/-S: 'auto', a list from the config, a single number,
// or a comma/space separated string from the command line (e.g. --weak_thresholds 0,0.001)
def threshold_args(value) {
    return value instanceof Collection ? value.join(' ') : value.toString().replace(',', ' ')
}
