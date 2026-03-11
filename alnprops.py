#!/usr/bin/env python3
import modules.paf as pafparser
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Estimate the proportions of aligned length to the respective target chromosomes for query sequences in the input PAF file.")
    parser.add_argument("-i", "--input", required=True, help="Input PAF file. Use '-' for stdin.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome file (one sequence name per line)")
    parser.add_argument("-q", "--query", required=True, help="Query genome file (one sequence name per line)")
    parser.add_argument("-o", "--output", help="Output tsv file (default: stdout). Format is one line per query/target combination with columns: query_name, target_name, aligned_length, query_length, target_length, aligned_length_to_query_proportion, aligned_length_to_target_proportion")
    return parser.parse_args()

# Get arguments
args = parse_args()
input_file = args.input if args.input != '-' else sys.stdin
reference_file = args.reference
query_file = args.query
output_file = args.output

# dictionary for storing results
results = {}

# Parse alignment file
paf = pafparser.read_paf(input_file, target_genome=reference_file, query_genome=query_file)

# for each query sequence in the alignments, calculate the total aligned length to each target chromosome and the total length of the query sequence
for aln in paf.alignments:
    query_name = aln.query_sequence
    target_name = aln.target_sequence
    aligned_length_query = aln.aln_length_query
    aligned_length_target = aln.aln_length_target
    query_sequence_length = aln.query_sequence_length
    target_sequence_length = aln.target_sequence_length

    if query_name not in results:
        results[query_name] = {}
    if target_name not in results[query_name]:
        results[query_name][target_name] = {'aligned_length': 0, 'query_length': query_sequence_length, 'target_length': target_sequence_length,
        'query_proportion': 0, 'target_proportion': 0}
    
    results[query_name][target_name]['aligned_length'] += aligned_length_query
    results[query_name][target_name]['query_proportion'] = results[query_name][target_name]['aligned_length'] / query_sequence_length
    results[query_name][target_name]['target_proportion'] = results[query_name][target_name]['aligned_length'] / target_sequence_length


# write results to output file
with open(output_file, 'w') if output_file else sys.stdout as f:
    f.write("query_sequence\ttarget_sequence\taligned_length\tquery_sequence_length\ttarget_sequence_length\taln_prop_query\taln_prop_target\n")
    for query_name in results:
        for target_name in results[query_name]:
            aligned_length = results[query_name][target_name]['aligned_length']
            query_length = results[query_name][target_name]['query_length']
            target_length = results[query_name][target_name]['target_length']
            query_proportion = results[query_name][target_name]['query_proportion']
            target_proportion = results[query_name][target_name]['target_proportion']
            f.write(f"{query_name}\t{target_name}\t{aligned_length}\t{query_length}\t{target_length}\t{query_proportion:.4f}\t{target_proportion:.4f}\n")
