import modules.paf as pafparser

class HomologyBlocks:
    def __init__(self, paf):
        self.block_1 = [paf.alignments[0]]
        self.block_2 = []
        self.blocks = []
    def is_consecutive(self, aln, block):
        if len(block) == 0:
            return False
        last_aln = block[-1]
        if aln.is_consecutive(last_aln):
            return True
        return False
    def block_length(self, block):
        length = 0
        for aln in block:
            length += aln.aln_length_target
        return length
    def add_block(self, block, max_gap=10000):
        # check that only one ref and query sequence is present in the block
        ref_seqs = set()
        query_seqs = set()
        for aln in block:
            ref_seqs.add(aln.target_sequence)
            query_seqs.add(aln.query_sequence)
        if len(ref_seqs) > 1 or len(query_seqs) > 1:
            raise ValueError("Block contains multiple reference or query sequences")
        # if length is too short, do not add the block
        if self.block_length(block) < max_gap:
            return
        refchrom = list(ref_seqs)[0]
        querychrom = list(query_seqs)[0]
        block_start_ref = min([aln.target_sequence_start for aln in block])
        block_end_ref = max([aln.target_sequence_end for aln in block])
        block_start_query = min([aln.query_sequence_start for aln in block])
        block_end_query = max([aln.query_sequence_end for aln in block])
        block_length_ref = block_end_ref - block_start_ref
        block_length_query = block_end_query - block_start_query
        self.blocks.append({
            "refchrom": refchrom,
            "querychrom": querychrom,
            "block_start_ref": block_start_ref,
            "block_end_ref": block_end_ref,
            "block_length_ref": block_length_ref,
            "block_start_query": block_start_query,
            "block_end_query": block_end_query,
            "block_length_query": block_length_query,
            "num_alignments": len(block)
        })
    def check_blocks(self, min_block_length=100000):
        len_block_2 = self.block_length(self.block_2)
        if len_block_2 > min_block_length:
            self.add_block(self.block_1, min_block_length)
            self.block_1 = self.block_2
            self.block_2 = []


def identify_homology_blocks(paf, min_block_length=100000):
    homology_blocks = HomologyBlocks(paf)
    for aln in paf.alignments[1:]:
        # case 1: consecutive to block_1 → extend it
        if homology_blocks.is_consecutive(aln, homology_blocks.block_1):
            homology_blocks.block_1.append(aln)
            continue
        # case 2: consecutive to block_2 → extend it, *then* check if long enough to promote
        if len(homology_blocks.block_2) > 0 and homology_blocks.is_consecutive(aln, homology_blocks.block_2):
            homology_blocks.block_2.append(aln)
            # check after appending
            homology_blocks.check_blocks(min_block_length)
            continue
        # case 3: not consecutive to either block
        # first check if block_2 is long enough to replace block_1
        homology_blocks.check_blocks(min_block_length)
        # start a new block_2 with the new alignment
        homology_blocks.block_2 = [aln]
    # finalize at end
    homology_blocks.add_block(homology_blocks.block_1, min_block_length)
    homology_blocks.add_block(homology_blocks.block_2, min_block_length)

    return homology_blocks.blocks
