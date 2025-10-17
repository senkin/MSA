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


def perform_signature_attribution(selected_mutations, signatures, normalise_mutations=False, verbose=False):
    """Optimized NNLS attribution with cached computations."""
    if signatures.empty:
        if verbose:
            print('Zero signatures provided to NNLS.')
        n_stats = 11
        return np.nan, np.nan, np.nan, np.nan, [sum(selected_mutations)] + [np.nan] * (n_stats - 1)
    
    # Convert to numpy for speed
    sig_array = signatures.values
    mut_array = np.array(selected_mutations)
    
    # NNLS solve
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


def remove_weak_signatures(selected_mutations, initial_signatures, weak_threshold=0.01, 
                           similarity_index=-3, verbose=False):
    """Optimized signature removal with early stopping."""
    significant_signatures = initial_signatures.copy()
    
    # Calculate base similarity once
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[similarity_index]
    
    if verbose:
        print(f'Starting removal loop. Base similarity: {base_similarity}')
        print(f'Current signatures: {significant_signatures.columns.tolist()}')
    
    while True:
        sig_cols = significant_signatures.columns
        
        # Early exit conditions
        if len(sig_cols) <= 1:
            if verbose:
                print(f'Only one signature left: {sig_cols.tolist()}')
            break
        
        # Vectorized contribution calculation
        contributions = {}
        for signature in sig_cols:
            sigs_without = significant_signatures.drop(signature, axis=1)
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, sigs_without)
            contributions[signature] = base_similarity - stat_info[similarity_index]
        
        # Find weakest signature
        weakest_sig = min(contributions, key=contributions.get)
        min_contribution = contributions[weakest_sig]
        
        if min_contribution < weak_threshold:
            significant_signatures = significant_signatures.drop(weakest_sig, axis=1)
            # Update base similarity
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
            base_similarity = stat_info[similarity_index]
            
            if verbose:
                print(f'Dropped {weakest_sig}, {len(significant_signatures.columns)} signatures left')
        else:
            if verbose:
                print('All weak signatures removed')
            break
    
    # Final similarity
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    final_similarity = stat_info[similarity_index]
    
    return significant_signatures, final_similarity


def add_strong_signatures(selected_mutations, initial_signatures, all_available_signatures,
                         strong_threshold=0.05, similarity_index=-3, verbose=False):
    """Optimized signature addition with early stopping."""
    significant_signatures = initial_signatures.copy()
    remaining_signatures = all_available_signatures.drop(initial_signatures.columns, axis=1)
    
    if remaining_signatures.empty:
        _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
        return significant_signatures, stat_info[similarity_index]
    
    # Calculate base similarity
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[similarity_index]
    
    if verbose:
        print(f'Starting addition loop. Base similarity: {base_similarity}')
    
    while not remaining_signatures.empty:
        # Calculate contributions for all remaining signatures
        contributions = {}
        for signature in remaining_signatures.columns:
            sigs_with = significant_signatures.copy()
            sigs_with[signature] = remaining_signatures[signature]
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, sigs_with)
            contributions[signature] = stat_info[similarity_index] - base_similarity
        
        # Find strongest signature
        strongest_sig = max(contributions, key=contributions.get)
        max_contribution = contributions[strongest_sig]
        
        if max_contribution > strong_threshold:
            significant_signatures[strongest_sig] = remaining_signatures[strongest_sig]
            remaining_signatures = remaining_signatures.drop(strongest_sig, axis=1)
            
            # Update base similarity
            _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
            base_similarity = stat_info[similarity_index]
            
            if verbose:
                print(f'Added {strongest_sig}')
        else:
            if verbose:
                print('All strong signatures added')
            break
    
    # Sort and return
    significant_signatures = significant_signatures.reindex(sorted(significant_signatures.columns), axis=1)
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    
    return significant_signatures, stat_info[similarity_index]


def optimise_signatures(selected_mutations, initial_signatures, all_available_signatures,
                       strategy='removal', weak_threshold=0.01, strong_threshold=0.05,
                       similarity_index=-3, loops_limit=100, verbose=False):
    """Optimized signature optimization with convergence tracking."""
    significant_signatures = initial_signatures.copy()
    
    _, _, _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[similarity_index]
    
    if verbose:
        print(f'Initial signatures: {significant_signatures.columns.tolist()}')
        print(f'Initial similarity: {base_similarity}')
    
    if strategy == 'removal':
        significant_signatures, converging_similarity = remove_weak_signatures(
            selected_mutations, significant_signatures, weak_threshold, similarity_index, verbose)
    
    elif strategy == 'addition':
        significant_signatures, converging_similarity = add_strong_signatures(
            selected_mutations, significant_signatures, all_available_signatures,
            strong_threshold, similarity_index, verbose)
    
    elif strategy == 'add-remove':
        converging_similarity = base_similarity
        
        for loop_counter in range(1, loops_limit + 1):
            prev_similarity = converging_similarity
            
            # Add then remove
            significant_signatures, converging_similarity = add_strong_signatures(
                selected_mutations, significant_signatures, all_available_signatures,
                strong_threshold, similarity_index, verbose)
            
            significant_signatures, converging_similarity = remove_weak_signatures(
                selected_mutations, significant_signatures, weak_threshold, 
                similarity_index, verbose)
            
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
    
    if verbose:
        print(f'Final signatures: {significant_signatures.columns.tolist()}')
        print(f'Final similarity: {converging_similarity}')
    
    return significant_signatures


def process_samples_batch(input_mutations, signatures, sel_sig_nums, args):
    """Batch process samples for efficiency."""
    samples = input_mutations.columns
    signature_columns_list = signatures.columns[sel_sig_nums].tolist()
    
    # Pre-allocate output dataframes
    output_weights = pd.DataFrame(0.0, index=samples, columns=signature_columns_list)
    output_mutations = pd.DataFrame(0.0, index=samples, columns=signature_columns_list)
    output_stat_info = pd.DataFrame(0.0, index=samples, 
                                    columns=['Mutational burden', 'RSS', 'Chi2', 'R2', 
                                            'Cosine similarity', 'Correlation', 
                                            'Chebyshev similarity', 'L1 similarity', 
                                            'L2 similarity', 'L3 similarity', 
                                            'Jensen-Shannon similarity'])
    
    residuals_dataframe = pd.DataFrame(0.0, index=input_mutations.index, columns=samples)
    fitted_dataframe = pd.DataFrame(0.0, index=input_mutations.index, columns=samples)
    
    output_weights.index.name = output_mutations.index.name = output_stat_info.index.name = 'Sample'
    
    initial_signatures = signatures.iloc[:, sel_sig_nums]
    
    for sample in samples:
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
            
            signatures_to_fill = final_signatures.columns
            output_weights.loc[sample, signatures_to_fill] = normalised_weights
            output_mutations.loc[sample, signatures_to_fill] = mutation_numbers
            output_stat_info.loc[sample] = stat_info
            residuals_dataframe[sample] = residuals
            fitted_dataframe[sample] = fitted
    
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
