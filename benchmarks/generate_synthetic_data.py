#!/usr/bin/env python3
import argparse
import os
import random
import sys
from pathlib import Path
from itertools import combinations
from typing import Callable, Any, ContextManager
import gzip
from collections import defaultdict


def get_pair_from_index(k: int, n_samples: int) -> tuple[int, int]:
    """
    Map a 1D combination index k to a unique pair (i, j) with 0 <= i < j < n_samples.
    Uses binary search over the combinations sequence.
    """
    low = 0
    high = n_samples - 2
    i = 0
    while low <= high:
        mid = (low + high) // 2
        s_mid = mid * n_samples - mid * (mid + 1) // 2
        if s_mid <= k:
            i = mid
            low = mid + 1
        else:
            high = mid - 1

    s_i = i * n_samples - i * (i + 1) // 2
    j = i + 1 + (k - s_i)
    return i, j


def get_indices_to_write(n_pairs: int, total_possible_pairs: int):
    """
    Generator yielding index values to write.
    Uses O(1) auxiliary memory for all pairs, and O(min(n_pairs, total_possible_pairs - n_pairs)) for subsets.
    """
    if n_pairs == total_possible_pairs:
        for idx in range(total_possible_pairs):
            yield idx
    elif n_pairs < total_possible_pairs // 2:
        chosen = set()
        while len(chosen) < n_pairs:
            chosen.add(random.randint(0, total_possible_pairs - 1))
        for idx in sorted(chosen):
            yield idx
    else:
        n_exclude = total_possible_pairs - n_pairs
        excluded = set()
        while len(excluded) < n_exclude:
            excluded.add(random.randint(0, total_possible_pairs - 1))
        for idx in range(total_possible_pairs):
            if idx not in excluded:
                yield idx


def get_stats_from_freq(freq: dict[int, int]) -> tuple[int, float, int]:
    """
    Calculate the min, median, and max of segment counts from a frequency map in O(1) memory.
    """
    total = sum(freq.values())
    if total == 0:
        return 0, 0.0, 0
    min_val = min(k for k, v in freq.items() if v > 0)
    max_val = max(k for k, v in freq.items() if v > 0)

    sorted_freq = sorted(freq.items())
    if total % 2 == 1:
        target_idx = total // 2
        current_idx = 0
        median_val = 0.0
        for val, count in sorted_freq:
            if current_idx <= target_idx < current_idx + count:
                median_val = float(val)
                break
            current_idx += count
        return min_val, median_val, max_val
    else:
        target_idx1 = total // 2 - 1
        target_idx2 = total // 2
        val1, val2 = None, None
        current_idx = 0
        for val, count in sorted_freq:
            if val1 is None and current_idx <= target_idx1 < current_idx + count:
                val1 = val
            if val2 is None and current_idx <= target_idx2 < current_idx + count:
                val2 = val
            if val1 is not None and val2 is not None:
                break
            current_idx += count
        median_val = (val1 + val2) / 2.0
        return min_val, median_val, max_val


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic IBD data (.genome and segments.txt) for benchmarking."
    )
    parser.add_argument(
        "--samples", type=int, default=100, help="Number of unique individuals (N)."
    )
    parser.add_argument(
        "--pairs",
        type=int,
        default=0,
        help="Total pair records to write. Defaults to all N(N-1)/2 pairs.",
    )
    parser.add_argument(
        "--related-fraction",
        type=float,
        default=0.05,
        help="Fraction of pairs that are related (have PI_HAT >= threshold).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="synthetic_data",
        help="Target output directory.",
    )
    parser.add_argument(
        "--gzip",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="gzip the output files to save space",
    )

    args = parser.parse_args()

    args.output_dir.mkdir(exist_ok=True)

    n_samples = args.samples
    # Generate IIDs
    iids = [f"id{i}" for i in range(1, n_samples + 1)]
    fid = "sim_family"

    # We are just going to shadow the open function if the user wishes to write gzipped files
    open_func: Callable[..., ContextManager[Any]] = gzip.open if args.gzip else open

    mode = "wt" if args.gzip else "w"

    suffix = ".gz" if args.gzip else ""
    total_possible_pairs = n_samples * (n_samples - 1) // 2
    n_pairs = args.pairs if args.pairs > 0 else total_possible_pairs
    assert (
        n_pairs <= total_possible_pairs
    ), f"Error: maximum number of possible pairs: {total_possible_pairs}, but user requested {n_pairs}. Make sure the pairs argument is lower than the number of total possible pairs"
    n_pairs = min(n_pairs, total_possible_pairs)

    n_related = int(n_pairs * args.related_fraction)

    genome_file = (
        args.output_dir / f"synthetic_input_{n_samples}_samples.genome{suffix}"
    )
    segments_file = (
        args.output_dir / f"synthetic_segments_{n_samples}_samples.txt{suffix}"
    )

    print(f"Generating synthetic dataset with:")
    print(f"  - Individuals: {n_samples}")
    print(f"  - Total Pairs: {n_pairs}")
    print(f"  - Related Pairs (approx): {n_related}")
    print(f"  - Output genome file: {genome_file}")
    print(f"  - Output segments file: {segments_file}")

    # Related configurations
    rel_types = [
        # (Z0, Z1, Z2, PI_HAT, min_segments, max_segments)
        (0.0, 1.0, 0.0, 0.50, 1, 2),  # Parent-Offspring
        (0.25, 0.50, 0.25, 0.50, 2, 4),  # Full Siblings
        (0.50, 0.50, 0.0, 0.25, 1, 3),  # Second Degree
        (0.75, 0.25, 0.0, 0.125, 1, 2),  # Third Degree
    ]

    # Header format matches PLINK output
    genome_header = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n"
    segments_header = "iid1\tiid2\tstart\tend\tcmlen\tchrom\tibd\n"
    segment_counts_freq = defaultdict(int)

    # Set up Knuth's Algorithm S to select exactly n_related related pairs
    related_selected = 0
    pairs_written = 0

    with open_func(genome_file, mode=mode) as genome_fh, open_func(
        segments_file, mode=mode
    ) as segment_fh:
        genome_fh.write(genome_header)
        segment_fh.write(segments_header)

        for k in get_indices_to_write(n_pairs, total_possible_pairs):
            rem_pairs = n_pairs - pairs_written
            rem_related = n_related - related_selected

            if rem_related > 0 and random.random() < rem_related / rem_pairs:
                is_related = True
                related_selected += 1
            else:
                is_related = False

            pairs_written += 1

            iid1, iid2 = get_pair_from_index(k, n_samples)

            if is_related:
                # Select a random relationship profile
                z0, z1, z2, pi_hat, min_seg, max_seg = random.choice(rel_types)
                phe = -1
                dst = round(random.uniform(0.70, 0.85), 6)
                ppc = 1.0000
                ratio = round(random.uniform(2.0, 15.0), 4)

                # Write matching segment records
                n_segs = random.randint(min_seg, max_seg)
                segment_counts_freq[n_segs] += 1
                for _ in range(n_segs):
                    chrom = random.randint(1, 22)
                    start = random.randint(100000, 100000000)
                    end = start + random.randint(500000, 20000000)
                    cmlen = round((end - start) / 1000000 * random.uniform(1.2, 1.8), 2)
                    ibd_state = 1 if z2 == 0.0 else random.choice([1, 2])

                    segment_fh.write(
                        f"{iid1}\t{iid2}\t{start}\t{end}\t{cmlen}\t{chrom}\t{ibd_state}\n"
                    )
            else:
                # Unrelated pair profile
                z0, z1, z2, pi_hat = 1.0, 0.0, 0.0, 0.0
                phe = -1
                dst = round(random.uniform(0.67, 0.69), 6)
                ppc = round(random.uniform(0.1, 0.5), 4)
                ratio = "NA"

            # Format genome file spacing
            ratio_str = f"{ratio:>8}" if ratio != "NA" else "      NA"
            line = f"  {fid} {iid1}  {fid} {iid2} OT     0  {z0:.4f}  {z1:.4f}  {z2:.4f}  {pi_hat:.4f}  {phe:2d}  {dst:.6f}  {ppc:.4f} {ratio_str}\n"
            genome_fh.write(line)

    min_segs, median_segs, max_segs = get_stats_from_freq(segment_counts_freq)

    print("Done generating synthetic benchmarking files.")
    if sum(segment_counts_freq.values()) > 0:
        print(f"Segment Statistics (per related pair):")
        print(f"  - Min segments:    {min_segs}")
        print(f"  - Median segments: {median_segs}")
        print(f"  - Max segments:    {max_segs}")


if __name__ == "__main__":
    main()
