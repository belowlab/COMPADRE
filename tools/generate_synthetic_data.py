#!/usr/bin/env python3
import argparse
import os
import random
import sys
from pathlib import Path
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic IBD data (.genome and segments.txt) for benchmarking.")
    parser.add_argument("--samples", type=int, default=100, help="Number of unique individuals (N).")
    parser.add_argument("--pairs", type=int, default=0, help="Total pair records to write. Defaults to all N(N-1)/2 pairs.")
    parser.add_argument("--related-fraction", type=float, default=0.05, help="Fraction of pairs that are related (have PI_HAT >= threshold).")
    parser.add_argument("--output-dir", type=Path, default="synthetic_data", help="Target output directory.")
    
    args = parser.parse_args()

    args.output_dir.mkdir(exist_ok=True)

    n_samples = args.samples
    # Generate IIDs
    iids = [f"id{i}" for i in range(1, n_samples + 1)]
    fid = "sim_family"

    # Generate list of all possible unique pairs
    all_pairs = list(combinations(range(0,args.samples), 2))

    total_possible_pairs = len(all_pairs)
    n_pairs = args.pairs if args.pairs > 0 else total_possible_pairs
    assert n_pairs <= total_possible_pairs, f"Error: maximum number of possible pairs: {total_possible_pairs}, but user requested {n_pairs}. Make sure the pairs argument is lower than the number of total possible pairs"
    n_pairs = min(n_pairs, total_possible_pairs)
    
    # Sample pairs if limit requested
    if n_pairs < total_possible_pairs:
        pairs_to_write = random.sample(all_pairs, n_pairs)
    else:
        pairs_to_write = all_pairs
        
    n_related = int(n_pairs * args.related_fraction)
    
    # Shuffle so related pairs are distributed
    random.shuffle(pairs_to_write)
    
    genome_file = args.output_dir / "synthetic_input.genome"
    segments_file = args.output_dir / "synthetic_segments.txt"
    
    print(f"Generating synthetic dataset with:")
    print(f"  - Individuals: {n_samples}")
    print(f"  - Total Pairs: {n_pairs}")
    print(f"  - Related Pairs (approx): {n_related}")
    print(f"  - Output genome file: {genome_file}")
    print(f"  - Output segments file: {segments_file}")
    
    # Related configurations
    rel_types = [
        # (Z0, Z1, Z2, PI_HAT, min_segments, max_segments)
        (0.0, 1.0, 0.0, 0.50, 1, 2),    # Parent-Offspring
        (0.25, 0.50, 0.25, 0.50, 2, 4), # Full Siblings
        (0.50, 0.50, 0.0, 0.25, 1, 3),   # Second Degree
        (0.75, 0.25, 0.0, 0.125, 1, 2)  # Third Degree
    ]
    
    # Header format matches PLINK output
    genome_header = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n"
    segments_header = "iid1\tiid2\tstart\tend\tcmlen\tchrom\tibd\n"
    
    with open(genome_file, "w") as genome_fh, open(segments_file, "w") as segment_fh:
        genome_fh.write(genome_header)
        segment_fh.write(segments_header)
        
        for idx, (iid1, iid2) in enumerate(pairs_to_write):
            is_related = idx < n_related
            
            if is_related:
                # Select a random relationship profile
                z0, z1, z2, pi_hat, min_seg, max_seg = random.choice(rel_types)
                phe = -1
                dst = round(random.uniform(0.70, 0.85), 6)
                ppc = 1.0000
                ratio = round(random.uniform(2.0, 15.0), 4)
                
                # Write matching segment records
                n_segs = random.randint(min_seg, max_seg)
                for _ in range(n_segs):
                    chrom = random.randint(1, 22)
                    start = random.randint(100000, 100000000)
                    end = start + random.randint(500000, 20000000)
                    cmlen = round((end - start) / 1000000 * random.uniform(1.2, 1.8), 2)
                    ibd_state = 1 if z2 == 0.0 else random.choice([1, 2])
                    
                    segment_fh.write(f"{iid1}\t{iid2}\t{start}\t{end}\t{cmlen}\t{chrom}\t{ibd_state}\n")
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
            
    print("Done generating synthetic benchmarking files.")

if __name__ == "__main__":
    main()
