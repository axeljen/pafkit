#!/usr/bin/env python3

import modules.paf as pafparser
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Wedge, FancyBboxPatch, PathPatch
from matplotlib.path import Path
import numpy as np

"""
Axel Jensen, 2025

Quick syntenymap plotter for PAF files.

Minimal usage:
python3 plot_paf.py -i input.paf -r reference_genome_file.txt -q query_genome_file.txt -o output.png

"""

def synteny_plot(paf, output_file, ref_color="lightgray", query_color="lightgray", aln_color="orange", aln_color_inverted="red",
                 ref_color_alpha=1.0, query_color_alpha=1.0, aln_color_alpha=0.4, aln_color_inverted_alpha=0.4, ref_labels=True, query_labels=True,
                 refgenome_label=None, querygenome_label=None):
    """ plot synteny as a ribbon plot using matplotlib """
    # set up figure
    fig, ax = plt.subplots(figsize=(10, 6), dpi = 300)
    ax.set_xlim(-2, max(paf.target_genome.get_maxpos(), paf.query_genome.get_maxpos()) + 1)
    # y-positions of lower edge of reference and query rectangles
    y_ref = 2
    y_query = 1
    # heights of chromosomes
    rect_height = 0.15
    # this also gives us the ymin and ymax of connecting polygons
    pol_ymin = y_query + rect_height
    pol_ymax = y_ref
    ax.set_ylim(y_query - .5, y_ref + .5)
    # plot reference sequences as rectangles (top)
    for sequence, length in paf.target_genome.sequences.items():
        start = paf.target_genome.cumulative_startpos[sequence]
        length_scaled = length
        rect = Rectangle((start, y_ref), length_scaled, rect_height, facecolor=ref_color, edgecolor="k")
        ax.add_patch(rect)
        if label_refchroms:
            ax.text(start + length_scaled / 2, y_ref + rect_height + 0.1, sequence, ha="center", fontsize=6, rotation=45, va="bottom")    
    for sequence, length in paf.query_genome.sequences.items():
        start = paf.query_genome.cumulative_startpos[sequence] 
        length_scaled = length
        rect = Rectangle((start, y_query), length_scaled, rect_height, facecolor=query_color, edgecolor="k")
        ax.add_patch(rect)
        if label_querychroms:
            ax.text(start + length_scaled / 2, y_query - rect_height - 0.1, sequence, ha="center", fontsize=6, rotation=45, va="bottom")
    if refgenome_label:
        ax.text(-2, y_ref + rect_height / 2, refgenome_label, ha="right", fontsize=8, va="bottom", fontweight='bold')
    if querygenome_label:
        ax.text(-2, y_query + rect_height/ 2 , querygenome_label, ha="right", fontsize=8, va="top", fontweight='bold')
    # plot alignments as ribbons
    for i, aln in enumerate(paf.alignments):  # adapt if paf stores alignments differently
        # reference coords (top)
        x1 = paf.target_genome.cumulative_startpos[aln.target_sequence] + aln.target_sequence_start
        x2 = paf.target_genome.cumulative_startpos[aln.target_sequence] + aln.target_sequence_end
        y1 = pol_ymax
        # query coords (bottom)
        x3 = paf.query_genome.cumulative_startpos[aln.query_sequence] + aln.query_sequence_start
        x4 = paf.query_genome.cumulative_startpos[aln.query_sequence] + aln.query_sequence_end
        y2 = pol_ymin
        # polygon connecting ref block to query block
        if aln.target_sequence_strand == '+':
            poly = Polygon([[x1, y1], [x2, y1], [x4, y2], [x3, y2]],
                           closed=True, facecolor=aln_color, alpha=aln_color_alpha, edgecolor=None)
            #draw_sigmoid_ribbon(ax, x1, x2, y1, x3, x4, y2, color=aln_color, alpha=aln_color_alpha)
        else:
            poly = Polygon([[x1, y1], [x2, y1], [x4, y2], [x3, y2]],
                           closed=True, facecolor=aln_color_inverted, alpha=aln_color_inverted_alpha, edgecolor=None)
            #draw_sigmoid_ribbon(ax, x1, x2, y1, x3, x4, y2, color=aln_color_inverted, alpha=aln_color_inverted_alpha)
        ax.add_patch(poly)    
    # clean up
    ax.axis("off")
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    return True


def parse_args():
    parser = argparse.ArgumentParser(description="Plot synteny map from PAF file.")
    parser.add_argument("-i", "--input", required=True, help="Input PAF file. Use '-' for stdin.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome file (one sequence name per line)")
    parser.add_argument("-q", "--query", required=True, help="Query genome file (one sequence name per line)")
    parser.add_argument("-o", "--output", required=True, help="Output image file (e.g., output.png)")
    parser.add_argument("-l","--min-aln-length", type=int, default=1000, help="Minimum alignment length to consider (default: 1000)")
    parser.add_argument("-id", "--min-identity", type=float, default=0, help="Minimum fraction identity to consider (default: 0 (bound between 0-1))")
    parser.add_argument("-qual", "--min-quality", type=int, default=0, help="Minimum mapping quality to consider (default: 0)")
    parser.add_argument("-reflen", "--min-reference-length", type=int, default=0, help="Minimum reference sequence length to consider (default: 0)")
    parser.add_argument("-qlen","--min-query-length", type=int, default=0, help="Minimum query sequence length to consider (default: 0)")
    parser.add_argument("-gap", "--gap-size", type=int, default=1000000, help="Gap size in bp between sequences in the plot (default: 1000000)")
    parser.add_argument("--chromlabs", choices=['both', 'ref', 'query', 'none'], default='both', help="Show chromosome labels: both, ref, query, none (default: both)")
    parser.add_argument("--skip-reorientation", action='store_true', help="Skip reorienting query sequences to match reference (default: False)")
    parser.add_argument("--refgenome-label", type=str, default=None, help="Label for reference genome (default: None)")
    parser.add_argument("--querygenome-label", type=str, default=None, help="Label for query genome (default: None)")
    return parser.parse_args()

# Get arguments
args = parse_args()
input_file = args.input
reference_file = args.reference
query_file = args.query
output_file = args.output
min_aln_length = args.min_aln_length
min_identity = args.min_identity
if min_identity > 1 or min_identity < 0:
    raise ValueError("Minimum identity must be between 0 and 1.")
min_quality = args.min_quality
min_reference_length = args.min_reference_length
min_query_length = args.min_query_length
gap_size = args.gap_size
chromlabs = args.chromlabs
skip_reorientation = args.skip_reorientation
refgenome_label = args.refgenome_label
querygenome_label = args.querygenome_label

if chromlabs == 'none':
    label_refchroms = False
    label_querychroms = False
elif chromlabs == 'both':
    label_refchroms = True
    label_querychroms = True
elif chromlabs == 'ref':
    label_refchroms = True
    label_querychroms = False
elif chromlabs == 'query':
    label_refchroms = False
    label_querychroms = True

# Parse alignment file
paf = pafparser.read_paf(input_file, target_genome=reference_file, query_genome=query_file)

# filter
paf.filter_alignments(minlen=min_aln_length, minmapq=min_quality, minpid=min_identity, only_ref_chroms=True, only_query_chroms=True, min_target_sequence_length=min_reference_length, min_query_sequence_length=min_query_length)

# make sure that it's sorted
paf.sort_on_target()

# reorient query sequences to match reference sequences
if not skip_reorientation:
    for seq in paf.query_genome.sequences.keys():
        dominating_strand = paf.get_dominating_strand(seq)
        if dominating_strand == '-':
            paf.flip_query_alignments(seq)

# add cumulative positions
paf.add_cumulative_query_start(gap_size)


# plot
synteny_plot(paf, output_file, ref_color="lightgray", query_color="lightgray", aln_color="lightgreen", aln_color_inverted="red",
             refgenome_label=refgenome_label, querygenome_label=querygenome_label)

