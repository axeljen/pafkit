import sys

# sort alignments into blocks consisting of things thar are on the same chromosomes
def identify_blocks(paf, min_length=1000000):
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
    # clean them up a bit
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
        # if strand is negative, swap start and end coordinates of query
        if block[0].target_sequence_strand == '-':
            min_query, max_query = paf.query_genome.sequences[query_sequence] - max_query, paf.query_genome.sequences[query_sequence] - min_query
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
    cleaned_blocks = [block for block in cleaned_blocks if block["length_target"] >= min_length and block["length_query"] >= min_length]
    print("Identified {} blocks after filtering for minimum length of {}".format(len(cleaned_blocks), min_length))
    print(cleaned_blocks)
    return cleaned_blocks

def call_fissions(blocks):
    # sort blocks by target sequence and start position
    blocks = sorted(blocks, key=lambda x: (x["target_sequence"], x["min_target"]))
    fissions = []
    last_block = blocks[0]
    for block in blocks[1:]:
        if block["target_sequence"] == last_block["target_sequence"]:
            # same target sequence, check if query sequence is different and non-overlapping
            if block["query_sequence"] != last_block["query_sequence"]:
                fissions.append({
                    "event": "fission",
                    "target_sequence_1": block["target_sequence"],
                    "target_min_1": last_block["max_target"],
                    "target_max_1": block["min_target"],
                    "query_sequence_1": last_block["query_sequence"],
                    "query_min_1": last_block["min_query"],
                    "query_max_1": last_block["max_query"],
                    "query_sequence_2": block["query_sequence"],
                    "query_min_2": block["min_query"],
                    "query_max_2": block["max_query"],
                    "target_sequence_2": "NA",
                    "target_min_2": "NA",
                    "target_max_2": "NA"
                })
                last_block = block
            else:
                # update last block if overlapping or same query
                if block["max_query"] > last_block["max_query"]:
                    last_block = block
        else:
            # different target sequence, reset
            last_block = block
    return fissions
                

def call_fusions(blocks):
    # Fusions are called in the same way, only query sorted instead and query sequences checked for being the same
    blocks = sorted(blocks, key=lambda x: (x["query_sequence"], x["min_query"]))
    fusions = []
    last_block = blocks[0]
    for block in blocks[1:]:
        if block["query_sequence"] == last_block["query_sequence"]:
            # same query sequence, check if target sequence is different and non-overlapping
            if block["target_sequence"] != last_block["target_sequence"]:
                fusions.append({
                    "event": "fusion",
                    "target_sequence_1": last_block["target_sequence"],
                    "target_min_1": last_block["min_target"],
                    "target_max_1": last_block["max_target"],
                    "query_sequence_1": block["query_sequence"],
                    "query_min_1": last_block["max_query"],
                    "query_max_1": block["min_query"],
                    "target_sequence_2": block["target_sequence"],
                    "target_min_2": block["min_target"],
                    "target_max_2": block["max_target"],
                    "query_sequence_2": "NA",
                    "query_min_2": "NA",
                    "query_max_2": "NA"
                })
                last_block = block
            else:
                # update last block if overlapping or same target
                if block["max_target"] > last_block["max_target"]:
                    last_block = block
        else:
            # different query sequence, reset
            last_block = block
    return fusions

def write_events(outfile, events):
    headers = ['event', 'target_sequence_1', 'target_min_1', 'target_max_1', 'query_sequence_1', 'query_min_1', 'query_max_1', 'target_sequence_2', 'target_min_2', 'target_max_2', 'query_sequence_2', 'query_min_2', 'query_max_2']
    if outfile == sys.stdout:
        f = sys.stdout
    else:
        f = open(outfile, 'w')
    f.write("\t".join(headers) + "\n")
    for event in events:
        if event['event'] == 'fission':
            f.write("{event}\t{target_sequence}\t{target_min}\t{target_max}\t{query_sequence_1}\t{query_start_1}\t{query_end_1}\tNA\tNA\tNA\t{query_sequence_2}\t{query_start_2}\t{query_end_2}\n".format(
                event='fission',
                target_sequence=event['target_sequence_1'],
                target_min=event['target_min_1'],
                target_max=event['target_max_1'],
                query_sequence_1=event['query_sequence_1'],
                query_start_1=event['query_min_1'],
                query_end_1=event['query_max_1'],
                query_sequence_2=event['query_sequence_2'],
                query_start_2=event['query_min_2'],
                query_end_2=event['query_max_2']
            ))
        elif event['event'] == 'fusion':
            f.write("{event}\t{target_sequence_1}\t{target_start_1}\t{target_end_1}\t{query_sequence}\t{query_min}\t{query_max}\t{target_sequence_2}\t{target_start_2}\t{target_end_2}\tNA\tNA\tNA\n".format(
                event='fusion',
                target_sequence_1=event['target_sequence_1'],
                target_start_1=event['target_min_1'],
                target_end_1=event['target_max_1'],
                query_sequence=event['query_sequence_1'],
                query_min=event['query_min_1'],
                query_max=event['query_max_1'],
                target_sequence_2=event['target_sequence_2'],
                target_start_2=event['target_min_2'],
                target_end_2=event['target_max_2']
            ))
    if f is not sys.stdout:
        f.close()

