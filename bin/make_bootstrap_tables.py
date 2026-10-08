#!/usr/bin/env python3
""" make_bootstrap_tables.py
Module to make bootstrap tables and calculate metrics/scores for simulated data.
Refactored for improved efficiency and readability.
"""

from argparse import ArgumentParser
import os
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from common_methods import (
    make_folder_if_not_exists, 
    write_data_to_JSON, 
    calculate_confidence_interval,
    calculate_stat_scores, 
    calculate_sensitivity_CI, 
    calculate_specificity_CI,
    truth_table_path
)


def load_signature_table(signature_tables_path, signatures_prefix, mutation_type, context):
    """Load appropriate signature table based on mutation type and context."""
    if mutation_type == 'SBS':
        if context == 96:
            return pd.read_csv(
                f'{signature_tables_path}/{signatures_prefix}_{mutation_type}_signatures.csv',
                index_col=[0, 1]
            )
        elif context in [192, 288]:
            return pd.read_csv(
                f'{signature_tables_path}/{signatures_prefix}_{mutation_type}_{context}_signatures.csv',
                index_col=[0, 1, 2]
            )
        elif context in [1536, 4608]:
            return pd.read_csv(
                f'{signature_tables_path}/{signatures_prefix}_{mutation_type}_{context}_signatures.csv',
                index_col=0
            )
        else:
            raise ValueError(f"Context {context} is not supported.")
    else:
        return pd.read_csv(
            f'{signature_tables_path}/{signatures_prefix}_{mutation_type}_signatures.csv',
            index_col=0
        )


def calculate_sensitivity_thresholds(signatures_CPs_dict, signatures_to_consider, 
                                     signature_attribution_thresholds, confidence_level):
    """Calculate sensitivity thresholds for each signature."""
    sensitivity_thresholds = pd.DataFrame(columns=signatures_to_consider, dtype=float)
    
    lowest_threshold = signature_attribution_thresholds[0]
    highest_threshold = signature_attribution_thresholds[-1]
    
    for signature in signatures_to_consider:
        lowest_CP = signatures_CPs_dict[lowest_threshold].loc[0, signature]
        highest_CP = signatures_CPs_dict[highest_threshold].loc[0, signature]
        
        if lowest_CP >= confidence_level:
            sensitivity_thresholds.loc[0, signature] = lowest_threshold
        elif highest_CP < confidence_level:
            sensitivity_thresholds.loc[0, signature] = np.nan
        else:
            # Interpolate to find threshold
            for threshold_id, threshold in enumerate(signature_attribution_thresholds[:-1]):
                next_threshold = signature_attribution_thresholds[threshold_id + 1]
                lower_CP = signatures_CPs_dict[threshold].loc[0, signature]
                higher_CP = signatures_CPs_dict[next_threshold].loc[0, signature]
                
                if lower_CP < confidence_level <= higher_CP:
                    if higher_CP == lower_CP:
                        sensitivity_thresholds.loc[0, signature] = threshold
                    else:
                        # Linear interpolation
                        interpolated = threshold + (confidence_level - lower_CP) * \
                                     (next_threshold - threshold) / (higher_CP - lower_CP)
                        sensitivity_thresholds.loc[0, signature] = interpolated
                    break
    
    return sensitivity_thresholds


# Attribution scales produced for every run: relative (weights) and absolute
# (mutation counts). Maps the output-name suffix to the bootstrap file key.
ATTRIBUTION_SCALES = {'weights': 'weights', 'abs_mutations': 'abs'}


def read_bootstrap_tables(number_of_b_samples, bootstrap_files, samples,
                          signatures_to_consider, stat_metrics):
    """
    Read every bootstrap iteration once: BOTH the relative (weights) and absolute
    (mutation-count) attribution tables, so that confidence intervals and attribution
    distributions can be produced on both scales, and the stat metrics.

    Returns:
        attributions: {scale: array of shape (samples, signatures, iterations)}
        stat_values:  array of shape (samples, stat metrics, iterations)
    where scale is one of ATTRIBUTION_SCALES ('weights', 'abs_mutations'). The iterations of
    each sample and signature are contiguous, so the statistics over them are computed exactly
    as over a single array of values.
    """
    attributions = {scale: np.empty((len(samples), len(signatures_to_consider), number_of_b_samples))
                    for scale in ATTRIBUTION_SCALES}
    stat_values = np.empty((len(samples), len(stat_metrics), number_of_b_samples))

    for i in range(number_of_b_samples):
        for scale, file_key in ATTRIBUTION_SCALES.items():
            table = pd.read_csv(bootstrap_files[file_key](i + 1), index_col=0)
            attributions[scale][:, :, i] = table.loc[samples, signatures_to_consider].to_numpy(dtype=float)
        stat_table = pd.read_csv(bootstrap_files['stat'](i + 1), index_col=0)
        # a sample missing from a stat table gets NaN metrics
        stat_values[:, :, i] = stat_table[stat_metrics].reindex(samples).to_numpy(dtype=float)

    return attributions, stat_values


def distributions_per_sample(attributions, samples, signatures):
    """{sample: DataFrame(iterations x signatures)} of a (samples, signatures, iterations) array"""
    return {sample: pd.DataFrame(attributions[s].T, columns=signatures)
            for s, sample in enumerate(samples)}


def distributions_per_column(values, samples, columns):
    """{column: DataFrame(iterations x samples)} of a (samples, columns, iterations) array,
    e.g. per signature or per stat metric"""
    return {column: pd.DataFrame(values[:, c, :].T, columns=samples)
            for c, column in enumerate(columns)}


def calculate_confidence_intervals(attributions, central_attribution, samples, signatures,
                                   confidence_level):
    """
    Confidence intervals for all samples and signatures from bootstrap attributions of shape
    (samples, signatures, iterations). Returns the table of [central value, [lower, upper]]
    cells, and the lower and upper bounds as (samples, signatures) arrays.
    """
    lower, upper = calculate_confidence_interval(attributions, confidence=confidence_level * 100, axis=-1)

    # central values keep the type of their column in the table
    central = [central_attribution[signature].loc[samples].to_numpy() for signature in signatures]
    cells = np.empty((len(samples), len(signatures)), dtype=object)
    for s in range(len(samples)):
        for g in range(len(signatures)):
            cells[s, g] = [central[g][s], [lower[s, g], upper[s, g]]]

    return pd.DataFrame(cells, index=samples, columns=signatures), lower, upper


def calculate_signature_CPs(lower, upper, truth, samples, signatures_to_consider, acting_signatures):
    """
    Signature confidence probabilities: the fraction of samples whose confidence interval
    contains the truth, per acting signature, and pooled over the signatures absent from the
    simulation ('Others'). lower, upper and truth are (samples, signatures) arrays.
    """
    covered = (lower <= truth) & (truth <= upper)
    acting = np.array([signature in acting_signatures for signature in signatures_to_consider], dtype=bool)

    signatures_CPs = {signature: np.where(acting[g], covered[:, g], np.nan)
                      for g, signature in enumerate(signatures_to_consider)}
    if not acting.all():
        # pooled over all signatures absent from the simulation
        signatures_CPs['Others'] = covered[:, ~acting].mean(axis=1)
    signatures_CPs = pd.DataFrame(signatures_CPs, index=samples)

    # Calculate overall CPs
    return signatures_CPs.sum(axis=0).to_frame().T / len(samples)


def calculate_CPs_for_thresholds(lower, upper, truth, samples, signatures_to_consider, thresholds):
    """Calculate CPs for different sensitivity thresholds (lower, upper and truth as in
    calculate_signature_CPs)."""
    signatures_CPs_dict = {}

    for threshold in thresholds:
        # a confidence interval from zero reaches up to the threshold at least
        upper_bound = np.where((lower == 0) & (upper < threshold), threshold, upper)

        # covered if the CI contains the truth, or if both are positive
        covered = ((lower <= truth) & (truth <= upper_bound)) | \
                  ((lower > 0) & (upper_bound > 0) & (truth > 0))
        CPs = pd.DataFrame(covered.astype(float), index=samples, columns=signatures_to_consider)

        # Normalize to get probabilities
        signatures_CPs_dict[threshold] = CPs.sum(axis=0).to_frame().T / len(samples)

    return signatures_CPs_dict


def main():
    start_time = datetime.now()
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_attributions_folder", dest="input_attributions_folder", 
                       default='output_tables/')
    parser.add_argument("-S", "--signature_path", dest="signature_tables_path", 
                       default='signature_tables/')
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", 
                       default='sigProfiler')
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM')
    parser.add_argument("-o", "--output_folder", dest="output_folder", 
                       default='output_tables/')
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='')
    parser.add_argument("-c", "--context", dest="context", default=192, type=int)
    parser.add_argument("-l", "--confidence_level", dest="confidence_level", 
                       default=0.95, type=float)
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("-n", "--number_of_b_samples", dest="number_of_b_samples", 
                       default=1000, type=int)
    parser.add_argument("-N", "--number_of_samples", dest="number_of_samples", 
                       default=-1, type=int)
    parser.add_argument("-T", "--signature_attribution_thresholds", nargs='+', 
                       dest="signature_attribution_thresholds")
    parser.add_argument("--suffix", dest="suffix", default='')

    args = parser.parse_args()

    # Validate mutation type
    if not args.mutation_type:
        parser.error("Please specify mutation type using -t option (SBS, DBS, ID, SV, CNV).")
    if args.mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError(f"Unknown mutation type: {args.mutation_type}")

    # Set up paths
    if not args.suffix:
        input_folder = f'{args.input_attributions_folder}/{args.dataset_name}/'
        output_folder = f'{args.output_folder}/{args.dataset_name}'
    else:
        base_name = f'{args.dataset_name}_{args.context}_NNLS_{args.suffix}'
        input_folder = f'{args.input_attributions_folder}/{base_name}/'
        output_folder = f'{args.output_folder}/{base_name}/'
    
    make_folder_if_not_exists(output_folder)
    # simulated data come with their truth table (in the attributions folder)
    truth_file = truth_table_path(input_folder, args.dataset_name, args.mutation_type, args.context)
    is_simulated = os.path.exists(truth_file)
    if is_simulated:
        make_folder_if_not_exists(f'{output_folder}/truth_studies')

    print("*" * 50)
    print(f"Processing {args.mutation_type} mutation type, {args.dataset_name} dataset")
    print(f"Truth table {truth_file}: " + ("found, scoring against it" if is_simulated else "not found (real data)"))

    # Load central attribution tables
    central_attr_abs = pd.read_csv(
        f'{input_folder}/output_{args.dataset_name}_{args.mutation_type}_mutations_table.csv',
        index_col=0
    )
    central_attr_weights = pd.read_csv(
        f'{input_folder}/output_{args.dataset_name}_{args.mutation_type}_weights_table.csv',
        index_col=0
    )
    central_stat_table = pd.read_csv(
        f'{input_folder}/output_{args.dataset_name}_{args.mutation_type}_stat_info.csv',
        index_col=0
    )

    # Determine bootstrap file patterns
    bootstrap_files = {
        'abs': lambda i: f'{input_folder}/bootstrap_output/output_{args.dataset_name}_{args.mutation_type}{"_" + args.suffix if args.suffix else ""}_{i}_mutations_table.csv',
        'weights': lambda i: f'{input_folder}/bootstrap_output/output_{args.dataset_name}_{args.mutation_type}{"_" + args.suffix if args.suffix else ""}_{i}_weights_table.csv',
        'stat': lambda i: f'{input_folder}/bootstrap_output/output_{args.dataset_name}_{args.mutation_type}{"_" + args.suffix if args.suffix else ""}_{i}_stat_info.csv'
    }

    # Load signatures
    signatures = load_signature_table(
        args.signature_tables_path, 
        args.signatures_prefix, 
        args.mutation_type, 
        args.context
    )

    # Load truth table for simulated data
    truth_attribution_table = None
    if is_simulated:
        truth_attribution_table = pd.read_csv(truth_file, index_col=0)

    # Limit samples if specified
    if args.number_of_samples != -1:
        central_attr_weights = central_attr_weights.head(args.number_of_samples)
        central_attr_abs = central_attr_abs.head(args.number_of_samples)
        if is_simulated:
            truth_attribution_table = truth_attribution_table.head(args.number_of_samples)

    # Central attribution tables on both scales: relative (weights) and absolute (counts)
    central_by_scale = {
        'weights': central_attr_weights,
        'abs_mutations': central_attr_abs,
    }

    # Get samples and signatures
    samples = central_attr_weights.index.tolist()
    signatures_to_consider = list(central_attr_abs.columns)
    signatures = signatures[signatures_to_consider]
    stat_metrics = central_stat_table.columns.tolist()

    # Find acting signatures for simulated data
    acting_signatures = []
    if is_simulated:
        acting_signatures = [
            sig for sig in truth_attribution_table.columns
            if truth_attribution_table[sig].max() > 0
        ]

    print('Processing bootstrap samples...')
    # Read all bootstrap samples once, both relative and absolute tables
    attributions_by_scale, stat_values = read_bootstrap_tables(
        args.number_of_b_samples,
        bootstrap_files,
        samples,
        signatures_to_consider,
        stat_metrics
    )
    stat_metrics_dict = distributions_per_column(stat_values, samples, stat_metrics)

    print(f'Bootstrap processing complete. Elapsed: {datetime.now() - start_time}')

    # Confidence intervals and attribution distributions, written for BOTH scales
    print('Calculating confidence intervals...')
    CI_bounds_by_scale = {}
    for scale in ATTRIBUTION_SCALES:
        ci, lower, upper = calculate_confidence_intervals(
            attributions_by_scale[scale],
            central_by_scale[scale],
            samples,
            signatures_to_consider,
            args.confidence_level
        )
        CI_bounds_by_scale[scale] = (lower, upper)

        scale_filename = f'{args.dataset_name}_{args.mutation_type}_bootstrap_output_{scale}'
        ci.to_csv(f'{output_folder}/CIs_{scale_filename}.csv')
        write_data_to_JSON(
            distributions_per_sample(attributions_by_scale[scale], samples, signatures_to_consider),
            f'{output_folder}/attributions_per_sample_{scale_filename}.json'
        )
        write_data_to_JSON(
            distributions_per_column(attributions_by_scale[scale], samples, signatures_to_consider),
            f'{output_folder}/attributions_per_signature_{scale_filename}.json'
        )
        # stat metrics are scale-independent, but written under both names so the
        # plotting scripts can read them for either scale
        write_data_to_JSON(
            stat_metrics_dict,
            f'{output_folder}/stat_metrics_{scale_filename}.json'
        )

    print(f'Confidence intervals calculated. Elapsed: {datetime.now() - start_time}')

    # The relative (weights) scale is canonical for pruning and the truth studies below
    lower_CI, upper_CI = CI_bounds_by_scale['weights']
    attributions = attributions_by_scale['weights']
    central_attribution = central_attr_weights

    # Pruned attribution: prune on the CI lower bound, emit absolute counts.
    # The prune decision is scale-invariant - a signature is zero in the same
    # bootstrap iterations on both scales, and "percentile == 0" depends only on
    # the zero pattern, not the magnitudes - so the relative CI is used by
    # convention and yields the same table as the absolute CI would.
    central_abs = central_attr_abs.loc[samples, signatures_to_consider].to_numpy(dtype=float)
    pruned_attribution_table = pd.DataFrame(
        np.where(lower_CI == 0, 0, np.around(central_abs)),
        index=samples,
        columns=signatures_to_consider
    )

    # Calculate signature prevalences (presence fraction; identical on both scales)
    print('Calculating signature prevalences...')
    signatures_prevalences = pd.DataFrame(
        (attributions != 0).sum(axis=-1).astype(float),
        index=samples,
        columns=signatures_to_consider
    )
    signatures_prevalences = signatures_prevalences / args.number_of_b_samples

    print(f'Prevalences calculated. Elapsed: {datetime.now() - start_time}')

    # Simulated data specific calculations
    if is_simulated:
        print('Calculating metrics for simulated data...')

        # Lower CI bounds as an attribution table
        lower_CI_attributions = pd.DataFrame(lower_CI, index=samples, columns=signatures_to_consider)

        # Signatures never simulated have zero truth
        for signature in signatures_to_consider:
            if signature not in truth_attribution_table.columns:
                truth_attribution_table[signature] = 0
        truth = truth_attribution_table.loc[samples, signatures_to_consider].to_numpy(dtype=float)

        # Calculate signature CPs
        signatures_CPs = calculate_signature_CPs(
            lower_CI,
            upper_CI,
            truth,
            samples,
            signatures_to_consider,
            acting_signatures
        )

        # Calculate truth-measured differences, summed over the signatures in turn
        measured = {
            'Truth - mean': attributions.mean(axis=-1),
            'Truth - median': np.median(attributions, axis=-1),
            'Truth - NNLS': central_attribution.loc[samples, signatures_to_consider].to_numpy(dtype=float),
        }
        differences = {name: np.zeros(len(samples)) for name in measured}
        for g in range(len(signatures_to_consider)):
            for name, values in measured.items():
                differences[name] += truth[:, g] - values[:, g]
        truth_and_measured_difference = pd.DataFrame(differences, index=samples)

        # Calculate CPs for different thresholds
        if args.signature_attribution_thresholds:
            thresholds = [float(t) / 100 for t in args.signature_attribution_thresholds]
            signatures_CPs_dict = calculate_CPs_for_thresholds(
                lower_CI,
                upper_CI,
                truth,
                samples,
                signatures_to_consider,
                thresholds
            )

            sensitivity_thresholds = calculate_sensitivity_thresholds(
                signatures_CPs_dict,
                signatures_to_consider,
                thresholds,
                args.confidence_level
            )

        # Calculate statistical scores
        scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']
        
        stat_scores_tables = pd.DataFrame(columns=scores, dtype=float)
        stat_scores_from_CI_tables = pd.DataFrame(columns=scores, dtype=float)
        sensitivity_CI_tables = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
        specificity_CI_tables = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
        sensitivity_CI_from_CI_tables = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
        specificity_CI_from_CI_tables = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
        
        stat_scores_tables.loc[0] = calculate_stat_scores(
            signatures, central_attribution, truth_attribution_table
        )
        stat_scores_from_CI_tables.loc[0] = calculate_stat_scores(
            signatures, lower_CI_attributions, truth_attribution_table
        )
        sensitivity_CI_tables.loc[0] = calculate_sensitivity_CI(
            signatures, central_attribution, truth_attribution_table
        )
        sensitivity_CI_from_CI_tables.loc[0] = calculate_sensitivity_CI(
            signatures, lower_CI_attributions, truth_attribution_table
        )
        specificity_CI_tables.loc[0] = calculate_specificity_CI(
            signatures, central_attribution, truth_attribution_table
        )
        specificity_CI_from_CI_tables.loc[0] = calculate_specificity_CI(
            signatures, lower_CI_attributions, truth_attribution_table
        )
        
        # Per-signature scores
        stat_scores_per_sig = {}
        stat_scores_from_CI_per_sig = {}
        sensitivity_CI_per_sig = {}
        sensitivity_CI_from_CI_per_sig = {}
        specificity_CI_per_sig = {}
        specificity_CI_from_CI_per_sig = {}
        signatures_scores = {score: pd.DataFrame(columns=acting_signatures, dtype=float) 
                            for score in scores}
        signatures_scores_from_CI = {score: pd.DataFrame(columns=acting_signatures, dtype=float) 
                                     for score in scores}
        
        for signature in signatures_to_consider:
            stat_scores_per_sig[signature] = pd.DataFrame(columns=scores, dtype=float)
            stat_scores_from_CI_per_sig[signature] = pd.DataFrame(columns=scores, dtype=float)
            sensitivity_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
            sensitivity_CI_from_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
            specificity_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
            specificity_CI_from_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL', 'upper_CL'], dtype=float)
            
            stat_scores_per_sig[signature].loc[0] = calculate_stat_scores(
                [signature], central_attribution, truth_attribution_table
            )
            stat_scores_from_CI_per_sig[signature].loc[0] = calculate_stat_scores(
                [signature], lower_CI_attributions, truth_attribution_table
            )
            sensitivity_CI_per_sig[signature].loc[0] = calculate_sensitivity_CI(
                [signature], central_attribution, truth_attribution_table
            )
            sensitivity_CI_from_CI_per_sig[signature].loc[0] = calculate_sensitivity_CI(
                [signature], lower_CI_attributions, truth_attribution_table
            )
            specificity_CI_per_sig[signature].loc[0] = calculate_specificity_CI(
                [signature], central_attribution, truth_attribution_table
            )
            specificity_CI_from_CI_per_sig[signature].loc[0] = calculate_specificity_CI(
                [signature], lower_CI_attributions, truth_attribution_table
            )
            
            for score in scores:
                signatures_scores[score].loc[0, signature] = \
                    stat_scores_per_sig[signature].loc[0, score]
                signatures_scores_from_CI[score].loc[0, signature] = \
                    stat_scores_from_CI_per_sig[signature].loc[0, score]
        
        print(f'Metrics calculated. Elapsed: {datetime.now() - start_time}')

    # Write outputs
    print('Writing outputs...')

    # Shared, scale-independent outputs (the per-scale CIs / attribution
    # distributions / stat_metrics for both 'weights' and 'abs_mutations' were
    # already written above).
    signatures_prevalences.to_csv(
        f'{output_folder}/signatures_prevalences_{args.dataset_name}_{args.mutation_type}.csv'
    )
    pruned_attribution_table.to_csv(
        f'{output_folder}/pruned_attribution_{args.dataset_name}_{args.mutation_type}_abs_mutations.csv'
    )

    # Simulated data outputs
    if is_simulated:
        truth_studies_path = f'{output_folder}/truth_studies'
        truth_and_measured_difference.to_csv(
            f'{truth_studies_path}/truth_and_measured_difference_{args.mutation_type}.csv'
        )
        signatures_CPs.to_csv(
            f'{truth_studies_path}/signatures_CPs_{args.mutation_type}.csv'
        )
        
        if args.signature_attribution_thresholds:
            sensitivity_thresholds.to_csv(
                f'{truth_studies_path}/sensitivity_thresholds_{args.mutation_type}.csv'
            )
            write_data_to_JSON(
                signatures_CPs_dict,
                f'{truth_studies_path}/signatures_CPs_dict_{args.mutation_type}.json'
            )
        
        stat_scores_tables.to_csv(
            f'{truth_studies_path}/stat_scores_tables_{args.mutation_type}.csv'
        )
        stat_scores_from_CI_tables.to_csv(
            f'{truth_studies_path}/stat_scores_from_CI_tables_{args.mutation_type}.csv'
        )
        sensitivity_CI_tables.to_csv(
            f'{truth_studies_path}/sensitivity_CI_tables_{args.mutation_type}.csv'
        )
        sensitivity_CI_from_CI_tables.to_csv(
            f'{truth_studies_path}/sensitivity_CI_from_CI_tables_{args.mutation_type}.csv'
        )
        specificity_CI_tables.to_csv(
            f'{truth_studies_path}/specificity_CI_tables_{args.mutation_type}.csv'
        )
        specificity_CI_from_CI_tables.to_csv(
            f'{truth_studies_path}/specificity_CI_from_CI_tables_{args.mutation_type}.csv'
        )
        
        write_data_to_JSON(signatures_scores, 
                          f'{truth_studies_path}/signatures_scores_{args.mutation_type}.json')
        write_data_to_JSON(signatures_scores_from_CI,
                          f'{truth_studies_path}/signatures_scores_from_CI_{args.mutation_type}.json')
        write_data_to_JSON(stat_scores_per_sig,
                          f'{truth_studies_path}/stat_scores_per_sig_{args.mutation_type}.json')
        write_data_to_JSON(stat_scores_from_CI_per_sig,
                          f'{truth_studies_path}/stat_scores_from_CI_per_sig_{args.mutation_type}.json')
        write_data_to_JSON(sensitivity_CI_per_sig,
                          f'{truth_studies_path}/sensitivity_CI_per_sig_{args.mutation_type}.json')
        write_data_to_JSON(sensitivity_CI_from_CI_per_sig,
                          f'{truth_studies_path}/sensitivity_CI_from_CI_per_sig_{args.mutation_type}.json')
        write_data_to_JSON(specificity_CI_per_sig,
                          f'{truth_studies_path}/specificity_CI_per_sig_{args.mutation_type}.json')
        write_data_to_JSON(specificity_CI_from_CI_per_sig,
                          f'{truth_studies_path}/specificity_CI_from_CI_per_sig_{args.mutation_type}.json')

    print(f'All outputs written. Total elapsed: {datetime.now() - start_time}')


if __name__ == '__main__':
    main()
