import argparse
import os
import warnings
import pandas as pd
import numpy as np
import time
from scipy.optimize import nnls
from numba import jit
from common_methods import make_folder_if_not_exists, calculate_similarity, format_threshold

# Vectorized bootstrap methods
def bootstrap_mutation_table(input_dataframe, method="classic", fitted=None, residuals=None):
    """Optimized bootstrap with vectorized operations."""
    if method == "classic":
        return input_dataframe.sample(n=len(input_dataframe), replace=True)

    elif method == "bootstrap_residuals":
        if residuals is None or fitted is None:
            raise ValueError(f"Residuals and fitted values required for {method}")
        return fitted + bootstrap_mutation_table(residuals, method="classic")

    elif method == "poisson":
        # Vectorized Poisson sampling
        return input_dataframe * np.random.poisson(1, size=input_dataframe.shape)

    elif method in ["binomial", "multinomial", "multinomial_weight"]:
        sums = input_dataframe.sum(axis=0).values
        norm_mutations = input_dataframe.div(sums, axis=1).fillna(0).values
        bootstrap_data = np.zeros_like(input_dataframe.values)

        for i, col_sum in enumerate(sums):
            if method == "binomial":
                bootstrap_data[:, i] = np.random.binomial(int(col_sum), norm_mutations[:, i])
            elif method == "multinomial":
                bootstrap_data[:, i] = np.random.multinomial(int(col_sum), norm_mutations[:, i])
            elif method == "multinomial_weight":
                n = len(norm_mutations[:, i])
                counts = np.random.multinomial(n, [1/n] * n)
                bootstrap = np.repeat(input_dataframe.values[:, i], counts)
                np.random.shuffle(bootstrap)
                bootstrap_data[:, i] = bootstrap

        return pd.DataFrame(bootstrap_data, index=input_dataframe.index, columns=input_dataframe.columns)

    else:
        raise ValueError(f"Unknown bootstrap method: {method}")


@jit(nopython=True, cache=True)
def calculate_stats_numba(observed, fitted):
    """JIT-compiled statistics calculation for speed."""
    residuals = observed - fitted
    rss = np.sum(residuals ** 2)
    chi2 = np.sum(residuals ** 2 / (fitted + 1e-10))  # Avoid division by zero

    mean_obs = np.mean(observed)
    tss = np.sum((observed - mean_obs) ** 2)
    r2 = 1 - rss / tss if tss > 0 else 0

    return rss, chi2, r2


# Mapping from a stat_info index (see the layout assembled in
# perform_signature_attribution) to the (metric, normalise) pair understood by
# calculate_similarity. Index 8 (L2_normalised_by_first) is handled directly for
# speed and is the metric used by the optimisation loops throughout the pipeline.
_IDX_TO_METRIC = {
    4: ('Cosine', False),
    5: ('Correlation', False),
    6: ('Chebyshev', True),
    7: ('L1', True),
    9: ('L3', True),
    10: ('jensen-shannon', False),
}


def batched_nnls_shared_b(A, b, column_sets):
    """Solve NNLS(A[:, cols], b) for each ``cols`` in ``column_sets``.

    Every problem in the batch shares the same target ``b`` and selects columns
    from the same matrix ``A``. In the leave-one-out (removal) and add-one
    (addition) cases all problems have the same number of columns, so the batch
    is uniform.

    THIS IS THE GPU SWAP POINT FOR THE OPTIMISATION HOT PATH. The CPU
    implementation is a tight scipy loop; a batched GPU/cuML NNLS solver can
    replace this function wholesale (it receives ``A``, ``b`` and the per-problem
    column sets, i.e. everything a masked batched kernel needs).

    Parameters
    ----------
    A : ndarray (n_channels, n_signatures_total)
    b : ndarray (n_channels,)
    column_sets : sequence of int index sequences

    Returns
    -------
    list of ndarray : weights[i] has length len(column_sets[i]).
    """
    return [nnls(A[:, cols], b)[0] for cols in column_sets]


def batched_nnls_shared_A(A, B):
    """Solve NNLS(A, b) for every column ``b`` of ``B`` (shared design matrix).

    THIS IS THE GPU SWAP POINT FOR NON-OPTIMISED ATTRIBUTION, where every sample
    is fitted against the same signature matrix: one ``A``, many targets. The CPU
    implementation is a tight scipy loop; a batched GPU/cuML NNLS solver can
    replace it wholesale.

    Parameters
    ----------
    A : ndarray (n_channels, n_signatures)
    B : ndarray (n_channels, n_samples)

    Returns
    -------
    ndarray (n_signatures, n_samples) : weights for each sample (column).
    """
    weights = np.zeros((A.shape[1], B.shape[1]))
    for i in range(B.shape[1]):
        weights[:, i], _ = nnls(A, B[:, i])
    return weights


def _similarity_from_fitted(observed, fitted, idx, norm_obs):
    """Return ONLY the single similarity metric (position ``idx`` in the
    stat_info layout) from an already-computed fitted vector.

    Shared by the single-solve and batched paths so they use identical
    arithmetic. The default optimisation metric, L2_normalised_by_first (idx 8),
    is computed directly from the residual without scipy.spatial.distance.
    """
    if idx == 8:  # L2_normalised_by_first: 1 - ||observed - fitted|| / ||observed||
        return 1.0 - np.linalg.norm(observed - fitted) / norm_obs
    if idx in _IDX_TO_METRIC:
        metric, normalise = _IDX_TO_METRIC[idx]
        return calculate_similarity(observed, fitted, metric=metric, normalise=normalise)
    # Fallback for any other (non-similarity) index: reproduce the exact value the
    # original stat_info list would have held at that position.
    rss, chi2, r2 = calculate_stats_numba(observed, fitted)
    stat_info = [observed.sum(), rss, chi2, r2,
                 calculate_similarity(observed, fitted),
                 calculate_similarity(observed, fitted, metric='Correlation'),
                 calculate_similarity(observed, fitted, metric='Chebyshev', normalise=True),
                 calculate_similarity(observed, fitted, metric='L1', normalise=True),
                 calculate_similarity(observed, fitted, metric='L2_normalised_by_first'),
                 calculate_similarity(observed, fitted, metric='L3', normalise=True),
                 calculate_similarity(observed, fitted, metric='jensen-shannon')]
    return stat_info[idx]


def _attribution_similarity(sig_sub, observed, idx, norm_obs):
    """Solve a single NNLS problem and return only the optimisation metric.

    Used for the base/final solves; the per-step candidate batches go through
    batched_nnls_shared_b() instead.
    """
    weights, _ = nnls(sig_sub, observed)
    fitted = sig_sub @ weights
    return _similarity_from_fitted(observed, fitted, idx, norm_obs)


def perform_signature_attribution(selected_mutations, signatures, normalise_mutations=False,
                                  verbose=False, weights=None):
    """Optimized NNLS attribution with cached computations.

    Computes the full stat_info list (burden + all similarities). This is used
    for the final per-sample attribution that is written to disk; the iterative
    optimisation loops use _attribution_similarity() instead, which computes only
    the single metric they need.

    If ``weights`` is supplied (e.g. from a batched NNLS solve), the NNLS step is
    skipped and those weights are used directly. This lets the non-optimised path
    solve all samples through batched_nnls_shared_A() and still reuse this
    function for the per-sample post-processing.
    """
    if signatures.empty:
        if verbose:
            print('Zero signatures provided to NNLS.')
        n_stats = 11
        return np.nan, np.nan, np.nan, np.nan, [sum(selected_mutations)] + [np.nan] * (n_stats - 1)

    # Convert to numpy for speed
    sig_array = signatures.values
    mut_array = np.array(selected_mutations)

    # NNLS solve (skipped when weights are supplied by a batched solver)
    if weights is None:
        weights, _ = nnls(sig_array, mut_array)

    # Normalize weights
    weight_sum = weights.sum()
    normalised_weights = weights / weight_sum if weight_sum > 0 else weights

    mutation_numbers = (normalised_weights * mut_array.sum() if normalise_mutations
                       else weights)

    # Calculate fitted and residuals
    fitted = sig_array @ weights
    residuals = mut_array - fitted

    # Fast statistics with numba
    rss, chi2, r2 = calculate_stats_numba(mut_array, fitted)
    input_mutational_burden = mut_array.sum()

    # Similarity metrics (these are already optimized in common_methods)
    cosine_similarity = calculate_similarity(mut_array, fitted)
    correlation = calculate_similarity(mut_array, fitted, metric='Correlation')
    chebyshev_similarity = calculate_similarity(mut_array, fitted, metric='Chebyshev', normalise=True)
    L1_similarity = calculate_similarity(mut_array, fitted, metric='L1', normalise=True)
    L2_similarity = calculate_similarity(mut_array, fitted, metric='L2_normalised_by_first')
    L3_similarity = calculate_similarity(mut_array, fitted, metric='L3', normalise=True)
    jensenshannon_similarity = calculate_similarity(mut_array, fitted, metric='jensen-shannon')

    stat_info = [input_mutational_burden, rss, chi2, r2, cosine_similarity,
                 correlation, chebyshev_similarity, L1_similarity, L2_similarity,
                 L3_similarity, jensenshannon_similarity]

    if verbose:
        if chi2 > 1e10:
            print('************* High chi2 sample *************')
        print(f'Signatures: {signatures.columns.tolist()}')
        print(f'Stats: {stat_info}')

    return normalised_weights, mutation_numbers, fitted, residuals, stat_info


def remove_weak_signatures(observed, norm_obs, sig_values, active_cols, col_names,
                           weak_threshold=0.01, idx=8, verbose=False):
    """Greedy leave-one-out removal of weak signatures.

    Operates on integer column indices into ``sig_values`` (no pandas inside the
    loop) and evaluates a single similarity metric per candidate. Column order in
    ``active_cols`` is preserved, so the surviving set matches the original
    DataFrame-based implementation.

    Returns (active_cols, final_similarity).
    """
    active_cols = list(active_cols)
    base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)

    if verbose:
        print(f'Starting removal loop. Base similarity: {base_similarity}')
        print(f'Current signatures: {[col_names[c] for c in active_cols]}')

    while len(active_cols) > 1:
        # Batched leave-one-out solve: each candidate drops one signature, all
        # sharing the same target b (the GPU swap point for the hot path).
        column_sets = [active_cols[:p] + active_cols[p + 1:] for p in range(len(active_cols))]
        candidate_weights = batched_nnls_shared_b(sig_values, observed, column_sets)

        # Contribution of each signature = drop in similarity if it is removed.
        contributions = np.empty(len(active_cols))
        for p, cols in enumerate(column_sets):
            fitted = sig_values[:, cols] @ candidate_weights[p]
            contributions[p] = base_similarity - _similarity_from_fitted(observed, fitted, idx, norm_obs)

        p_weakest = int(np.argmin(contributions))
        if contributions[p_weakest] < weak_threshold:
            dropped = active_cols.pop(p_weakest)
            base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
            if verbose:
                print(f'Dropped {col_names[dropped]}, {len(active_cols)} signatures left')
        else:
            if verbose:
                print('All weak signatures removed')
            break

    if len(active_cols) <= 1 and verbose:
        print(f'Only one signature left: {[col_names[c] for c in active_cols]}')

    final_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    return active_cols, final_similarity


def add_strong_signatures(observed, norm_obs, sig_values, active_cols, candidate_cols, col_names,
                          strong_threshold=0.05, idx=8, verbose=False):
    """Greedy addition of strong signatures drawn from ``candidate_cols``.

    Like remove_weak_signatures, this works on integer column indices and a
    single metric. The surviving set is returned sorted by signature name, to
    match the original implementation's reindex(sorted(columns)).

    Returns (active_cols, final_similarity).
    """
    active_cols = list(active_cols)
    active_set = set(active_cols)
    # Candidates not already in the model, keeping the original column ordering.
    remaining = [c for c in candidate_cols if c not in active_set]

    if not remaining:
        sim = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
        return active_cols, sim

    base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    if verbose:
        print(f'Starting addition loop. Base similarity: {base_similarity}')

    while remaining:
        # Batched add-one solve: each candidate adds one signature, all sharing
        # the same target b (the GPU swap point for the hot path).
        column_sets = [active_cols + [c] for c in remaining]
        candidate_weights = batched_nnls_shared_b(sig_values, observed, column_sets)

        contributions = np.empty(len(remaining))
        for p, cols in enumerate(column_sets):
            fitted = sig_values[:, cols] @ candidate_weights[p]
            contributions[p] = _similarity_from_fitted(observed, fitted, idx, norm_obs) - base_similarity

        p_strongest = int(np.argmax(contributions))
        if contributions[p_strongest] > strong_threshold:
            chosen = remaining.pop(p_strongest)
            active_cols.append(chosen)
            base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
            if verbose:
                print(f'Added {col_names[chosen]}')
        else:
            if verbose:
                print('All strong signatures added')
            break

    # Match the original: reindex(sorted(columns)) -> order surviving sigs by name.
    active_cols = sorted(active_cols, key=lambda c: col_names[c])
    final_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
    return active_cols, final_similarity


def optimise_signatures(selected_mutations, initial_signatures, all_available_signatures,
                       strategy='removal', weak_threshold=0.01, strong_threshold=0.05,
                       similarity_index=-3, loops_limit=100, verbose=False):
    """Optimise the signature set for one sample.

    The candidate search runs entirely on numpy arrays / integer column indices
    and evaluates a single similarity metric per NNLS solve (see
    _attribution_similarity). The result is rebuilt as a DataFrame once, with the
    same column ordering as the original implementation, so downstream output is
    unchanged.
    """
    observed = np.asarray(selected_mutations, dtype=float)
    norm_obs = np.linalg.norm(observed)
    idx = similarity_index % 11  # default -3 maps to 8 (L2_normalised_by_first)

    all_names = list(all_available_signatures.columns)
    name_to_col = {name: i for i, name in enumerate(all_names)}
    sig_values = np.asarray(all_available_signatures.values, dtype=float)
    candidate_cols = list(range(len(all_names)))

    active_cols = [name_to_col[name] for name in initial_signatures.columns]

    if not active_cols:
        return all_available_signatures.iloc[:, :0]

    if verbose:
        base_similarity = _attribution_similarity(sig_values[:, active_cols], observed, idx, norm_obs)
        print(f'Initial signatures: {[all_names[c] for c in active_cols]}')
        print(f'Initial similarity: {base_similarity}')

    if strategy == 'removal':
        active_cols, _ = remove_weak_signatures(
            observed, norm_obs, sig_values, active_cols, all_names,
            weak_threshold, idx, verbose)

    elif strategy == 'addition':
        active_cols, _ = add_strong_signatures(
            observed, norm_obs, sig_values, active_cols, candidate_cols, all_names,
            strong_threshold, idx, verbose)

    elif strategy == 'add-remove':
        converging_similarity = _attribution_similarity(
            sig_values[:, active_cols], observed, idx, norm_obs)

        for loop_counter in range(1, loops_limit + 1):
            prev_similarity = converging_similarity

            active_cols, converging_similarity = add_strong_signatures(
                observed, norm_obs, sig_values, active_cols, candidate_cols, all_names,
                strong_threshold, idx, verbose)

            active_cols, converging_similarity = remove_weak_signatures(
                observed, norm_obs, sig_values, active_cols, all_names,
                weak_threshold, idx, verbose)

            convergence_delta = abs(converging_similarity - prev_similarity)

            if verbose:
                print(f'Loop {loop_counter} done. Convergence delta: {convergence_delta}')

            if convergence_delta <= 0:
                break

            if loop_counter >= loops_limit:
                warnings.warn(f"Maximum iterations ({loops_limit}) reached")
                break
    else:
        raise ValueError(f'Unknown strategy: {strategy}')

    final_names = [all_names[c] for c in active_cols]
    if verbose:
        print(f'Final signatures: {final_names}')

    return all_available_signatures.loc[:, final_names]


def process_samples_batch(input_mutations, signatures, sel_sig_nums, args):
    """Process all samples, accumulating results in preallocated numpy arrays.

    Per-sample results are written into numpy arrays and assembled into output
    DataFrames once at the end. The original implementation assigned into
    DataFrames row-by-row via .loc, which is slow at scale.
    """
    samples = list(input_mutations.columns)
    signature_columns_list = signatures.columns[sel_sig_nums].tolist()
    col_pos = {name: i for i, name in enumerate(signature_columns_list)}

    n_samples = len(samples)
    n_sigs = len(signature_columns_list)
    n_channels = len(input_mutations.index)

    # Pre-allocate output arrays (rows kept zero for skipped/empty samples)
    weights_arr = np.zeros((n_samples, n_sigs))
    mutations_arr = np.zeros((n_samples, n_sigs))
    stat_arr = np.zeros((n_samples, 11))
    residuals_arr = np.zeros((n_channels, n_samples))
    fitted_arr = np.zeros((n_channels, n_samples))

    initial_signatures = signatures.iloc[:, sel_sig_nums]

    def _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted):
        weights_arr[s_i, cols] = normalised_weights
        mutations_arr[s_i, cols] = mutation_numbers
        stat_arr[s_i, :] = stat_info
        residuals_arr[:, s_i] = residuals
        fitted_arr[:, s_i] = fitted

    if not args.optimise_signatures and not initial_signatures.empty:
        # Non-optimised: every sample shares the same signature matrix, so all
        # NNLS solves go through a single batched call (the shared-A GPU swap
        # point). Per-sample post-processing reuses perform_signature_attribution
        # via its weights= argument, keeping output identical to the per-sample path.
        sig_array = np.asarray(initial_signatures.values, dtype=float)
        cols = [col_pos[name] for name in initial_signatures.columns]

        valid_s_i, valid_b = [], []
        for s_i, sample in enumerate(samples):
            selected_mutations = input_mutations[sample].values
            if selected_mutations.sum() <= 0:
                warnings.warn(f"Sample {sample}: Zero mutations, skipping")
                continue
            valid_s_i.append(s_i)
            valid_b.append(np.asarray(selected_mutations, dtype=float))

        if valid_b:
            B = np.column_stack(valid_b)
            weights_batch = batched_nnls_shared_A(sig_array, B)
            for j, s_i in enumerate(valid_s_i):
                normalised_weights, mutation_numbers, fitted, residuals, stat_info = \
                    perform_signature_attribution(B[:, j], initial_signatures,
                                                normalise_mutations=args.normalise_mutations,
                                                verbose=args.verbose, weights=weights_batch[:, j])
                _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted)
    else:
        for s_i, sample in enumerate(samples):
            selected_mutations = input_mutations[sample].values

            if selected_mutations.sum() <= 0:
                warnings.warn(f"Sample {sample}: Zero mutations, skipping")
                continue

            if args.optimise_signatures:
                final_signatures = optimise_signatures(
                    selected_mutations, initial_signatures, signatures,
                    strategy=args.optimisation_strategy,
                    weak_threshold=args.weak_threshold,
                    strong_threshold=args.strong_threshold,
                    verbose=args.verbose)
            else:
                final_signatures = initial_signatures

            if not final_signatures.empty:
                normalised_weights, mutation_numbers, fitted, residuals, stat_info = \
                    perform_signature_attribution(selected_mutations, final_signatures,
                                                normalise_mutations=args.normalise_mutations,
                                                verbose=args.verbose)

                cols = [col_pos[name] for name in final_signatures.columns]
                _store(s_i, cols, normalised_weights, mutation_numbers, stat_info, residuals, fitted)

    stat_columns = ['Mutational burden', 'RSS', 'Chi2', 'R2',
                    'Cosine similarity', 'Correlation', 'Chebyshev similarity',
                    'L1 similarity', 'L2 similarity', 'L3 similarity',
                    'Jensen-Shannon similarity']

    output_weights = pd.DataFrame(weights_arr, index=samples, columns=signature_columns_list)
    output_mutations = pd.DataFrame(mutations_arr, index=samples, columns=signature_columns_list)
    output_stat_info = pd.DataFrame(stat_arr, index=samples, columns=stat_columns)
    residuals_dataframe = pd.DataFrame(residuals_arr, index=input_mutations.index, columns=samples)
    fitted_dataframe = pd.DataFrame(fitted_arr, index=input_mutations.index, columns=samples)

    output_weights.index.name = output_mutations.index.name = output_stat_info.index.name = 'Sample'

    return output_weights, output_mutations, output_stat_info, residuals_dataframe, fitted_dataframe


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM')
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='')
    parser.add_argument("-c", "--context", dest="context", default=96, type=int)
    parser.add_argument("-i", "--input_path", dest="input_path", default='input_mutation_tables/')
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/')
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler')
    parser.add_argument("-o", "--output_path", dest="output_path", default='output_tables/')
    parser.add_argument("-x", "--optimise_signatures", dest="optimise_signatures", action="store_true")
    parser.add_argument("-W", "--weak_threshold", dest="weak_threshold", default=0.01, type=float)
    parser.add_argument("-S", "--strong_threshold", dest="strong_threshold", default=0.01, type=float)
    parser.add_argument("-N", "--normalise_mutations", dest="normalise_mutations", action="store_true")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int)
    parser.add_argument("-B", "--bootstrap", dest="bootstrap", action="store_true")
    parser.add_argument("--bootstrap_method", dest="bootstrap_method", default='binomial')
    parser.add_argument("--add_suffix", dest="add_suffix", action="store_true")
    parser.add_argument("--optimisation_strategy", dest="optimisation_strategy", default='removal')

    args = parser.parse_args()

    # Setup paths
    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    output_path = f"{args.output_path}/{dataset_name}"

    if args.add_suffix:
        output_path += f'_{context}_NNLS'
        if args.optimise_signatures:
            weak_str = format_threshold(args.weak_threshold)
            strong_str = format_threshold(args.strong_threshold)
            output_path += f'_{weak_str}_{strong_str}'
        else:
            output_path += '_unoptimised'

    make_folder_if_not_exists(output_path)

    # Validate mutation type
    if not mutation_type:
        parser.error("Please specify mutation type using -t option (SBS, DBS, ID, SV, CNV)")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError(f"Unknown mutation type: {mutation_type}")

    # Load data based on mutation type
    sig_path = args.signature_tables_path
    inp_path = args.input_path
    sig_prefix = args.signatures_prefix

    if mutation_type == 'SBS':
        if context == 96:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                    sep=',', index_col=[0,1])
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=[0,1])
        elif context in [192, 288]:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_{context}_signatures.csv',
                                    sep=',', index_col=[0,1,2])
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=[0,1,2])
        elif context in [1536, 4608]:
            signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_{context}_signatures.csv',
                                    sep=',', index_col=0)
            input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{context}.csv',
                                         sep=',', index_col=0)
        else:
            raise ValueError(f"Context {context} is not supported")
    elif mutation_type == 'DBS':
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.dinucs.csv',
                                     sep=',', index_col=0)
    elif mutation_type == 'ID':
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.indels.csv',
                                     sep=',', index_col=0)
    else:  # SV and CNV
        signatures = pd.read_csv(f'{sig_path}/{sig_prefix}_{mutation_type}_signatures.csv',
                                sep=',', index_col=0)
        input_mutations = pd.read_csv(f'{inp_path}/{dataset_name}/WGS_{dataset_name}.{mutation_type}.csv',
                                     sep=',', index_col=0)

    print(f"Performing NNLS for {dataset_name} dataset, {mutation_type} mutation type.")
    if mutation_type == 'SBS':
        print(f"SBS context: {context}")
    if args.optimise_signatures:
        print(f"Optimised NNLS: weak/strong thresholds = {args.weak_threshold}/{args.strong_threshold}")
    if args.bootstrap:
        print(f"Bootstrap method: {args.bootstrap_method}")

    # Bootstrap if requested
    if args.bootstrap and args.bootstrap_method != "bootstrap_residuals":
        input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method)

    # Limit samples if specified
    if args.number_of_samples != -1:
        input_mutations = input_mutations.iloc[:, :args.number_of_samples]

    # Align indices
    input_mutations = input_mutations.reindex(signatures.index)

    # Handle bootstrap_residuals method
    if args.bootstrap and args.bootstrap_method == "bootstrap_residuals":
        residuals_df = pd.read_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_residuals.csv',
                                   index_col=[0,1,2])
        fitted_df = pd.read_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_fitted_values.csv',
                               index_col=[0,1,2])
        input_mutations = bootstrap_mutation_table(input_mutations, method=args.bootstrap_method,
                                                   fitted=fitted_df, residuals=residuals_df)

    num_ref_sigs = signatures.shape[1]
    sel_sig_nums = list(range(num_ref_sigs))

    print(f"Analyzing signatures: {signatures.columns[sel_sig_nums].tolist()}")

    # Process samples
    start_time = time.process_time()

    output_weights, output_mutations, output_stat_info, residuals_df, fitted_df = \
        process_samples_batch(input_mutations, signatures, sel_sig_nums, args)

    end_time = time.process_time()
    elapsed = end_time - start_time
    n_samples = len(input_mutations.columns)
    print(f"Attribution took {elapsed:.2f}s ({elapsed/n_samples:.2f}s per sample)")

    # Save outputs
    output_weights.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_weights_table.csv')
    output_mutations.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_mutations_table.csv')
    output_stat_info.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_stat_info.csv')
    residuals_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_residuals.csv')
    fitted_df.to_csv(f'{output_path}/output_{dataset_name}_{mutation_type}_fitted_values.csv')

    if not args.bootstrap:
        residuals_df.to_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_residuals.csv')
        fitted_df.to_csv(f'{inp_path}/{dataset_name}/output_{dataset_name}_{mutation_type}_fitted_values.csv')
