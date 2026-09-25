"""Compare the signature attributions of two runs (e.g. CPU vs GPU, fp64 vs mixed).

Finds every output_*_weights_table.csv under the reference directory (central
attributions and bootstrap_output/ iterations alike), pairs it with the file at the
same relative path under the candidate directory and reports, per file and overall,
how many samples have any signature weight differing by more than --tolerance.

Exit status: 0 if the fraction of differing samples is within --max-differing-fraction
and every reference file has a candidate counterpart, 1 otherwise. --exact requires
identical weights (tolerance 0, no differing samples), e.g. for determinism checks.

Examples:
    python bin/compare_attributions.py cpu/SIM_test_GPU gpu_fp64/SIM_test_GPU
    python bin/compare_attributions.py gpu_fp64 gpu_fp64_repeat --exact
    python bin/compare_attributions.py run_cpu/output_tables run_gpu/output_tables --json report.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def compare_tables(reference, candidate, tolerance, support_threshold):
    """Per-sample comparison of two weights tables (samples x signatures)."""
    ref = pd.read_csv(reference, index_col=0)
    cand = pd.read_csv(candidate, index_col=0)
    if set(ref.index) != set(cand.index):
        raise ValueError(f'{candidate}: samples differ from {reference}')
    signatures = ref.columns.union(cand.columns)
    ref = ref.reindex(columns=signatures, fill_value=0.0)
    cand = cand.reindex(index=ref.index, columns=signatures, fill_value=0.0)
    diff = np.abs(ref.values - cand.values)
    support_changed = (ref.values > support_threshold) != (cand.values > support_threshold)
    max_per_sample = diff.max(axis=1)
    differing = max_per_sample > tolerance
    return {
        'samples': int(len(ref)),
        'differing_samples': int(differing.sum()),
        'support_changed_samples': int(support_changed.any(axis=1).sum()),
        'max_abs_difference': float(diff.max()) if diff.size else 0.0,
        'worst_sample': str(ref.index[int(np.argmax(max_per_sample))]) if diff.size else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--tolerance', type=float, default=1e-3,
                        help='max |weight difference| for a sample to count as matching (default 1e-3)')
    parser.add_argument('--max-differing-fraction', type=float, default=0.02,
                        help='max fraction of differing samples, pooled over all files (default 0.02)')
    parser.add_argument('--support-threshold', type=float, default=1e-6,
                        help='weight above which a signature counts as attributed (default 1e-6)')
    parser.add_argument('--exact', action='store_true',
                        help='require identical weights (tolerance 0, no differing samples)')
    parser.add_argument('--json', type=Path, help='also write the report to this JSON file')
    args = parser.parse_args()
    if args.exact:
        args.tolerance, args.max_differing_fraction = 0.0, 0.0

    reference_files = sorted(args.reference.rglob('output_*_weights_table.csv'))
    if not reference_files:
        sys.exit(f'No output_*_weights_table.csv files found under {args.reference}')

    per_file, missing = {}, []
    for ref_file in reference_files:
        rel = ref_file.relative_to(args.reference)
        cand_file = args.candidate / rel
        if not cand_file.exists():
            missing.append(str(rel))
            continue
        per_file[str(rel)] = compare_tables(ref_file, cand_file, args.tolerance, args.support_threshold)

    total = sum(r['samples'] for r in per_file.values())
    differing = sum(r['differing_samples'] for r in per_file.values())
    fraction = differing / total if total else 0.0
    passed = not missing and fraction <= args.max_differing_fraction
    report = {
        'reference': str(args.reference), 'candidate': str(args.candidate),
        'tolerance': args.tolerance, 'max_differing_fraction': args.max_differing_fraction,
        'files_compared': len(per_file), 'files_missing': missing,
        'samples_compared': total, 'differing_samples': differing, 'differing_fraction': fraction,
        'support_changed_samples': sum(r['support_changed_samples'] for r in per_file.values()),
        'max_abs_difference': max((r['max_abs_difference'] for r in per_file.values()), default=0.0),
        'passed': passed, 'per_file': per_file,
    }

    print(f'{args.candidate} vs {args.reference}: {len(per_file)} weights tables, {total} sample attributions')
    print(f'  differing (> {args.tolerance:g}): {differing} ({100 * fraction:.2f}%, '
          f'allowed {100 * args.max_differing_fraction:.2f}%)')
    print(f"  attributed-signature set changed: {report['support_changed_samples']}; "
          f"max |weight difference|: {report['max_abs_difference']:.3g}")
    worst = sorted(per_file.items(), key=lambda kv: kv[1]['max_abs_difference'], reverse=True)[:3]
    for rel, r in worst:
        if r['differing_samples']:
            print(f"  {rel}: {r['differing_samples']}/{r['samples']} differing, "
                  f"worst sample {r['worst_sample']} ({r['max_abs_difference']:.3g})")
    for rel in missing:
        print(f'  MISSING in candidate: {rel}')
    print('  PASS' if passed else '  FAIL')

    if args.json:
        args.json.write_text(json.dumps(report, indent=2))
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
