#!/usr/bin/env python3

'''
Axel Jensen, 2025

This script provides a simple fission/fusion calling based on pairwise genome alignments. Note that this assumes that 
both genomes are assembled to chromosome level, and it is a relatively rough calling that doesn't take into account
inversions around breakpoints etc., which could make things a bit more complicated. The script offers a visualization of 
called events, and manual inspection is recommended to confirm events.

'''
import sys
from modules import paf as pafparser
import matplotlib.pyplot as plt
import numpy as np
import argparse
from modules.homology_blocks import identify_blocks,call_fissions, call_fusions, write_events
from modules.pafplotter import synteny_plot

# arguments
parser = argparse.ArgumentParser(description="Call fissions and fusions based on PAF alignments")
parser.add_argument("-i", "--input", required=True, help="Input PAF file with alignments between two genomes")
parser.add_argument("-o", "--output", default=sys.stdout, help="Output file for listing fission/fusion events (default: stdout)")
parser.add_argument("-p", "--plot", help="Optional output plot file showing fissions and fusions (e.g., output.png)")
parser.add_argument("-r", "--reference", required=True, help="Reference genome file with at least two tab separated columns, sequence name and length. Additional columns are ignored") 
parser.add_argument("-q", "--query", required=True, help="Query genome file with at least two tab separated columns, sequence name and length. Additional columns are ignored")
parser.add_argument("-gap", "---max-gap", type=int, default=1000000, help="Maximum gap between consecutive alignments to merge them (default: 1000000)")
parser.add_argument("-l","--min-aln-length", type=int, default=10000, help="Minimum alignment length to consider (default: 10000)")
parser.add_argument("-id", "--min-identity", type=float, default=0, help="Minimum fraction identity to consider (default: 0 (bound between 0-1))")
parser.add_argument("-qual", "--min-quality", type=int, default=0, help="Minimum mapping quality to consider (default: 0)")
args = parser.parse_args()

alignment_file = args.input
reference_file = args.reference
query_file = args.query
output_file = args.output
plot_file = args.plot
min_length = args.min_aln_length
min_identity = args.min_identity
if min_identity > 1 or min_identity < 0:
    raise ValueError("Minimum identity must be between 0 and 1.")
min_quality = args.min_quality
max_gap = args.max_gap

# Parse alignment file
paf = pafparser.read_paf(alignment_file, target_genome = reference_file, query_genome = query_file)
# filter
paf.filter_alignments(minlen=min_length, minmapq=min_quality, minpid=min_identity, only_ref_chroms=True, only_query_chroms=True)
# sort and merge
paf.sort_on_target()
paf.merge_consecutive(max_gap=max_gap)
# filter on length again, now using the max_gap argument since these are essentially raw homology blocks
paf.filter_alignments(minlen=max_gap)
# identify blocks
blocks = identify_blocks(paf, min_length=max_gap)
# call fissions and fusions
fissions = call_fissions(blocks)
fusions = call_fusions(blocks)

events = fissions + fusions
# output results
write_events(output_file, events)

# if plot
if plot_file:
    # convert fissions to segment
    segments = []
    for fission in fissions:
        segments.append({
            'sequence': fission['target_sequence_1'],
            'start': fission['target_min_1'],
            'end': fission['target_max_1'],
            'type': 'fission'
        })
    # reorient query sequences to match reference sequences
    for seq in paf.query_genome.sequences.keys():
        dominating_strand = paf.get_dominating_strand(seq)
        if dominating_strand == '-':
            paf.flip_query_alignments(seq)
    # add cumulative positions
    paf.add_cumulative_query_start(10000000)
    synteny_plot(paf, output_file=plot_file, reference_segments=segments, query_segments=None)
