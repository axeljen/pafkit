#!/usr/bin/env python3
import modules.paf as pafparser

import sys
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Plot synteny map from PAF file.")
    parser.add_argument("-i", "--input", required=True, help="Input PAF file")
    parser.add_argument("-o", "--output", help="Output paf file (defaults to stdout if not specified)", default=sys.stdout)
    parser.add_argument("-l", "--min-aln-length", type=int, default=1000, help="Minimum alignment length to consider (default: 1000)")
    parser.add_argument("--reference", type=str, default=None, help="Reference genome file (one sequence name per line) (default: None)")
    parser.add_argument("--query", type=str, default=None, help="Query genome file (one sequence name per line) (default: None)")
    parser.add_argument("-id", "--min-identity", type=float, default=0, help="Minimum fraction identity to consider (default: 0 (bound between 0-1))")
    parser.add_argument("-qual", "--min-quality", type=int, default=0, help="Minimum mapping quality to consider (default: 0)")
    parser.add_argument("-reflen", "--min-reference-length", type=int, default=0, help="Minimum reference sequence length to consider (default: 0)")
    parser.add_argument("-qlen", "--min-query-length", type=int, default=0, help="Minimum query sequence length to consider (default: 0)")
    parser.add_argument("-r", "--reference-region", type=str, default=None, help="Reference region to consider in the format 'chr or chr:start-end' (default: None). Coordinates are 1-based, inclusive.")
    parser.add_argument("-R", "--reference-regions-file", type=str, default=None, help="File with reference regions to consider, one per line in the format 'chr, chr:start-end or chr\tstart\tend' (default: None). Coordinates are 1-based, inclusive.")
    parser.add_argument("--only-contained", action='store_true', help="Only keep alignments fully contained within the specified reference region(s) (default: False)")
    parser.add_argument("--aligned-fraction", type=float, default=0.0, help="Minimum fraction of the reference region that must be covered by the alignment (default: 0.0 (bound between 0-1))")
    return parser.parse_args()

def parse_region(region_str):
    """ parse a region string in the format 'chr', 'chr:start-end' or 'chr\tstart\tend' """
    try:
        if ':' in region_str and '-' in region_str:
            chrom, positions = region_str.split(':')
            start, end = map(int, positions.split('-'))
            start = start - 1  # convert to 0-based
            end = end
        elif '\t' in region_str:
            chrom, start, end = region_str.split('\t')
            start, end = int(start) - 1, int(end)
        else:
            chrom = region_str
            start, end = None, None
    except ValueError:
        raise ValueError(f"Invalid region format: {region_str}")
    return (chrom, start, end)

def parse_fai_to_regions(fai_file):
    with open(fai_file, 'r') as f:
        regions = []
        for line in f:
            parts = line.strip().split('\t')
            regions.append((parts[0], None, None))
    return regions

# get arguments
args = parse_args()
input_paf = args.input
output_paf = args.output
min_aln_length = args.min_aln_length
min_identity = args.min_identity
min_quality = args.min_quality
min_reference_length = args.min_reference_length
min_query_length = args.min_query_length
reference_region = args.reference_region
reference_regions_file = args.reference_regions_file
# don't allow for both a single region and a file
if reference_region and reference_regions_file:
    raise ValueError("Cannot specify both a single reference region and a reference regions file.")
only_contained = args.only_contained
aligned_fraction = args.aligned_fraction
# read reference regions if specified
reference_regions = []
query=args.query
reference=args.reference

if reference_region:
    reference_regions.append(parse_region(reference_region))
if reference_regions_file:
    if reference_regions_file.endswith('.fai'):
        reference_regions = parse_fai_to_regions(reference_regions_file)
    else:
        with open(reference_regions_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    reference_regions.append(parse_region(line))
                
# read paf file and filter alignments
paf = pafparser.read_paf(input_paf, target_genome=reference, query_genome=query)

# start with pulling the regions of interest
subset_paf = pafparser.PAF()
# transfer header info
subset_paf.target_genome = paf.target_genome
subset_paf.query_genome = paf.query_genome
subset_paf.alignments = []
for region in reference_regions:
    chrom, start, end = region
    subset_paf.alignments.extend(paf.fetch(chrom, start, end, only_contained=only_contained))
    
# now apply the other filters
subset_paf.filter_alignments(minlen=min_aln_length, minpid=min_identity, minmapq=min_quality, min_target_sequence_length=min_reference_length, min_query_sequence_length=min_query_length, aligned_fraction=aligned_fraction)

# and write output
subset_paf.write_paf(output_paf)

