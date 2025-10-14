import modules.paf as pafparser
import modules.homology_blocks as hb
from importlib import reload

reload(hb)

alignment_file = "/Users/axeljensen/pafkit_testfiles/Cercopithecus_cephus_hap1_curated.FINAL.chromosomes_to_GCF_003339765.1_Mmul_10_genomic.paf"
query_file = "/Users/axeljensen/pafkit_testfiles/Cercopithecus_cephus_hap1_curated.FINAL.chromosomes.fa.fai"
reference_file = "/Users/axeljensen/pafkit_testfiles/Macaca_mulatta_mmul10_chromosomes.fa.fai"


paf = pafparser.read_paf(alignment_file, target_genome = reference_file, query_genome = query_file)

paf.filter_alignments(only_ref_chroms=True)
# start by sorting and merging consecutive alignments
paf.filter_alignments(only_ref_chroms=True)
paf.sort_on_target()

# merge consecutive alignments, use a generous gap length
max_gap = 1000000
paf.merge_consecutive(max_gap=max_gap)

# remove merged alignments that are shorter than threshold
min_length = 1000000
paf.filter_alignments(minlen=min_length)

# sort alignments into blocks consisting of things thar are on the same chromosomes
blocks = []

current_block = [paf.alignments[0]]

for aln in paf.alignments[1:]:
    if aln.target_sequence == current_block[-1].target_sequence and aln.query_sequence == current_block[-1].query_sequence:
        # append if both target and query sequences are the same
        current_block.append(aln)
    else:
        # otherwise, start a new block
        blocks.append(current_block)
        current_block = [aln]
# append the last block
blocks.append(current_block)

for block in blocks:
    print("{query}\t{qmin}\t{qmax}\t{target}\t{tmin}\t{tmax}\t{n}".format(
        query=block[0].query_sequence,
        qmin=min([aln.query_sequence_start for aln in block]),
        qmax=max([aln.query_sequence_end for aln in block]),
        target=block[0].target_sequence,
        tmin=min([aln.target_sequence_start for aln in block]),
        tmax=max([aln.target_sequence_end for aln in block]),
        n=len(block)
    ))

# clean up blocks
cleaned_blocks = []
for block in blocks:
    min_target = min([aln.target_sequence_start for aln in block])
    max_target = max([aln.target_sequence_end for aln in block])
    min_query = min([aln.query_sequence_start for aln in block])
    max_query = max([aln.query_sequence_end for aln in block])
    length_target = max_target - min_target
    length_query = max_query - min_query
    target_sequence = block[0].target_sequence
    query_sequence = block[0].query_sequence
    cleaned_blocks.append({
        "target_sequence": target_sequence,
        "query_sequence": query_sequence,
        "min_target": min_target,
        "max_target": max_target,
        "length_target": length_target,
        "min_query": min_query,
        "max_query": max_query,
        "length_query": length_query,
        "num_alignments": len(block)
    })

# remove blocks that are too short
min_length = 1000000
cleaned_blocks = [block for block in cleaned_blocks if block["length_target"] >= min_length and block["length_query"] >= min_length]
    
for block in cleaned_blocks:
    print("{query}\t{min_q}\t{max_q}\t{target}\t{min_t}\t{max_t}".format(
        query=block["query_sequence"],
        min_q=block["min_query"],
        max_q=block["max_query"],
        target=block["target_sequence"],
        min_t=block["min_target"],
        max_t=block["max_target"],
    ))