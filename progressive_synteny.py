#!/usr/bin/env python3

'''
Axel Jensen, 2025

This script will take one or more PAF alignment files and accompanying genome/fasta index files,
and generate a synteny plot with genomes progressively plotted along the Y axis. The backbone reference
genome is plotted on top, and all query genomes are expected to be aligned to that same reference.

'''

import sys
from modules import paf as pafparser
import matplotlib.pyplot as plt
import argparse
from importlib import reload
import modules.pafplotter as pp

def parse_config(config_file):
    parsed_data = {'queries': []}
    with open(config_file, 'r') as cf:
        active_query = 0
        for line in cf:
            if line.startswith('#') or line.strip() == '':
                continue
            text = line.split("#")[0].strip()  # remove comments
            if text.split('=')[0].strip() in ['reference', 'Reference']:
                parsed_data['reference'] = text.split('=')[1].strip()
            if text.strip() in ['query', 'Query']:
                active_query += 1
                parsed_data['queries'].append({'index': active_query})
            if text.split('=')[0].strip() in ['alignment', 'Alignment', 'paf', 'PAF', 'Paf']:
                parsed_data['queries'][-1]['paf'] = text.split('=')[1].strip()
            if text.split('=')[0].strip() in ['fai', 'Fai', 'genome', 'Genome']:
                parsed_data['queries'][-1]['fai'] = text.split('=')[1].strip()
            if text.split('=')[0].strip() in ['name', 'Name', 'label', 'Label']:
                parsed_data['queries'][-1]['label'] = text.split('=')[1].strip()
    return parsed_data

def parse_genomefile(genome_file, min_sequence_length=0):
    seq_lengths = {}
    with open(genome_file, 'r') as gf:
        for line in gf:
            parts = line.strip().split('\t')
            seq_name = parts[0]
            seq_length = int(parts[1])
            if seq_length >= min_sequence_length:
                seq_lengths[seq_name] = seq_length
    return seq_lengths

# a custum way of estimating median target position for query sequences, since we'll be using some non-standard class attributes here
def get_median_target_position_custom(paf, query_sequence, target_genome):
    median = []
    for aln in paf.fetch_query(query_sequence):
        startpos = aln.query_sequence_start + target_genome[aln.target_sequence_liftover]['']
            


# arguments
parser = argparse.ArgumentParser(description="Generate progressive synteny plot from PAF alignments")
parser.add_argument("-i", "--input-file", required=True, help="Input configuration file, see synmap_config.txt for an example.")
parser.add_argument("-o", "--output", required=True, help="Output plot file (e.g., output.png).")
# parameters for filtering input alignments
parser.add_argument("--min-length", type=int, default=1000, help="Minimum alignment length to consider.")
parser.add_argument("--min-identity", type=float, default=0.8, help="Minimum percent identity for alignments.")
parser.add_argument("--min-sequence-length", type=int, default=0, help="Minimum sequence length to consider from fasta index.")
parser.add_argument("--min-qual", type=int, default=0, help="Minimum mapping quality for alignments.")
parser.add_argument("--plot-gap", type=int, default=10000000, help="Gap size between genomes in the plot.")
parser.add_argument("--equal-widths", action='store_true', help="Use equal widths for all genomes in the plot. If given, the chromosomes of the first reference will be spaced using gapsize, whereas the remaining queries will be spaced such that the last chromosome ends in line with the last in the reference, regardless of gap size. Note that this argument can create some weird patterns if gap size is too small, leading to overlapping chromosomes in the query genomes.")

args = parser.parse_args()

input_file = args.input_file
output_file = args.output
min_sequence_length = args.min_sequence_length
min_alignment_length = args.min_length
min_identity = args.min_identity
min_mapping_quality = args.min_qual
plot_gap = args.plot_gap
equal_widths = args.equal_widths

# teststuff
input_file = "/cfs/klemming/home/a/axeljen/pafkit/config_hominins.txt"
output_file = "test_progressive_synteny.png" 
min_sequence_length = 10000000
min_alignment_length = 10000
min_identity = 0.8
min_mapping_quality = 20
plot_gap = 10000000
equal_widths = True
# read input configuration file
input_files = parse_config(input_file)

# parse reference genome
reference_genome = parse_genomefile(input_files['reference'], min_sequence_length=min_sequence_length)
# prep a dict for storing the parsed pafs
parsed_pafs = []

# parse each query paf
for query in input_files['queries']:
    print("Parsing PAF file for query genome: {}".format(query.get('label', f"query_{query['index']}")))
    query_label = query.get('label', f"query_{query['index']}")
    paf_file = query['paf']
    genome_file = query['fai']
    parsed_paf = pafparser.read_paf(paf_file,target_genome=input_files['reference'],query_genome=genome_file)
    # filter paf
    parsed_paf.filter_alignments(minlen=min_alignment_length, minpid=min_identity, minmapq=min_mapping_quality, min_query_sequence_length=min_sequence_length, min_target_sequence_length=min_sequence_length)
    parsed_pafs.append({'label': query_label, 'paf': parsed_paf})

# next step is to "translate" the pafs so that the target genome coordinates are lifted over to the previous query genome
# we'll store these in a new list
translated_pafs = [parsed_pafs[0]]  # first paf is unchanged

for i in range(1, len(parsed_pafs)):
    translated_paf = pafparser.PAF(alignments = [], target_genome=parsed_pafs[i-1]['paf'].query_genome, query_genome=parsed_pafs[i]['paf'].query_genome)
    print("Liftover alignments from {} to previous query genome {}.".format(parsed_pafs[i].get('label'), parsed_pafs[i-1].get('label')))
    for aln in parsed_pafs[i].get('paf').alignments:
        # liftover ref coordinates to previous query genome
        target_seq = aln.target_sequence
        target_start = aln.target_sequence_start
        target_end = aln.target_sequence_end
        liftover_start_raw = parsed_pafs[i-1].get('paf').get_query_position_from_target(target_seq, target_start)
        liftover_start_seq = liftover_start_raw[0]
        if liftover_start_raw[1] == liftover_start_raw[2]:
            liftover_start = liftover_start_raw[1]
        else:
            liftover_start = None  # ambiguous liftover
        liftover_end_raw = parsed_pafs[i-1].get('paf').get_query_position_from_target(target_seq, target_end)
        liftover_end_seq = liftover_end_raw[0]
        if liftover_end_raw[1] == liftover_end_raw[2]:
            liftover_end = liftover_end_raw[1]
        else:
            liftover_end = None  # ambiguous liftover
        if liftover_start_seq != liftover_end_seq:
            liftover_seq = None
            liftover_start = None
            liftover_end = None  # different sequences, cannot liftover
        else:
            liftover_seq = liftover_start_seq
        # add alignment coordinates lifted over to previous query genome
        aln.target_sequence_start_liftover = liftover_start
        aln.target_sequence_end_liftover = liftover_end
        aln.target_sequence_liftover = liftover_seq
        # copy alignment
        new_aln = pafparser.Alignment()
        new_aln.query_sequence = aln.query_sequence
        new_aln.query_sequence_start = aln.query_sequence_start
        new_aln.query_sequence_end = aln.query_sequence_end
        new_aln.target_sequence = aln.target_sequence_liftover
        new_aln.target_sequence_start = aln.target_sequence_start_liftover
        new_aln.target_sequence_end = aln.target_sequence_end_liftover
        new_aln.aln_mapping_quality = aln.aln_mapping_quality
        new_aln.aln_length_target = aln.aln_length_target
        new_aln.aln_length_query = aln.aln_length_query
        new_aln.aln_base_matches = aln.aln_base_matches
        new_aln.aln_base_total = aln.aln_base_total
        new_aln.target_sequence_strand = aln.target_sequence_strand
        new_aln.query_sequence_length = aln.query_sequence_length
        new_aln.target_sequence_length = aln.target_sequence_length
        # store the original target sequence as well, so we can use this to color the alignment plot later
        new_aln.original_target_sequence = aln.target_sequence
        # if liftover was successful, add alignment
        if new_aln.target_sequence is not None and new_aln.target_sequence_start is not None and new_aln.target_sequence_end is not None:
            translated_paf.alignments.append(new_aln)
    translated_pafs.append({'label': parsed_pafs[i]['label'], 'paf': translated_paf})

reload(pafparser)

# add cumulative start positions, first to the first reference/and query
translated_pafs[0]['paf'].add_cumulative_query_start(sequence_gap=plot_gap, equal_widths=equal_widths)

# if plotting with equal widths, fetch the maximum chromosome position of this one that'll be used in downstream gap calculations
if equal_widths:
    num_chroms_ref = len(translated_pafs[0]['paf'].query_genome.sequences)
    max_y_position = max(translated_pafs[0]['paf'].target_genome.cumulative_startpos[seq] + length for seq, length in translated_pafs[0]['paf'].target_genome.sequences.items())

# next, iteritavely go through the remaining alignments
# we'll update the cumulative starting order for the previous genome as we go, but start with the first query
target_order = list(sorted(translated_pafs[0]['paf'].query_genome.sequences, key=lambda x: translated_pafs[0]['paf'].query_genome.cumulative_startpos[x]))

for paf in translated_pafs[1:]:
    print("Adding cumulative start positions for query genome: {}".format(paf.get('label')))
    # if equal widths is set, calculate the gap for the query chromosomes based on the number of chromosomes
    if equal_widths:
        num_chroms = len(paf['paf'].target_genome.sequences)
        target_gap = max_y_position / num_chroms if num_chroms > 0 else 0
        paf['paf'].add_cumulative_query_start(sequence_gap=target_gap, equal_widths=equal_widths)
    else:
        paf['paf'].add_cumulative_query_start(sequence_gap=plot_gap, equal_widths=equal_widths)
    # update target order for next round
    target_order = list(sorted(paf['paf'].query_genome.sequences, key=lambda x: paf['paf'].query_genome.cumulative_startpos[x]))

# start a plotter object
reload(pp)
paf_plotter = pp.PafPlotter(translated_pafs)
paf_plotter.plot_genomes()
paf_plotter.plot_alignments()
paf_plotter.set_xlim()
paf_plotter.set_ylim()
output_file = "test.png"
paf_plotter.save(output_file)
# first the reference genome
genome_dicts['reference'] = {'label': 'reference', 'sequences': []}
for seq in reference_genome:
    seq_length = reference_genome[seq]
    genome_dicts['reference']['sequences'].append({'name': seq, 'length': seq_length, 'median_start': None, 'cumulative_start': None})

# add cumulative start to the reference, based on argument gap size
cumulative_pos = 0
for seq_dict in genome_dicts[0]['sequences']:
    seq_dict['cumulative_start'] = cumulative_pos
    cumulative_pos += seq_dict['length'] + plot_gap

# do the first query genome
genome_dicts.append({'label': parsed_pafs[0]['label'], 'sequences': []})
for seq in parsed_pafs[0]['paf'].query_genome.sequences:
    median_start = parsed_pafs[0]['paf'].get_median_target_position(seq)
    if seq in parsed_pafs[0]['paf'].query_genome.sequences:
        seq_length = parsed_pafs[0]['paf'].query_genome.sequences[seq]
    else:
        seq_length = 0
    genome_dicts[1]['sequences'].append({'name': seq, 'length': seq_length, 'median_start': median_start, 'cumulative_start': None})