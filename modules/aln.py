# class for storing alignments (as represented within a paf file)

# some terminology: 
## the target/reference genome is always referred to as "target", and the query genome as "query"
## scaffolds/contigs/chromosomes are always referred to as "sequence"
## Features related to the query are always prefixed with "query_", and those related to the target with "target_"
## Features that are related to the alignment itself are prefixed with "aln_"
## Coordinates are always 0-based, half-open.
class Alignment():
    def __init__(self):
        self.query_sequence = None
        self.query_sequence_length = None
        self.query_sequence_start = None
        self.query_sequence_end = None
        self.target_sequence_strand = None
        self.target_sequence = None
        self.target_sequence_length = None
        self.target_sequence_start = None
        self.target_sequence_end = None
        self.aln_length_target = None
        self.aln_length_query = None
        self.aln_base_matches = None
        self.aln_base_total = None
        self.aln_mapping_quality = None
        self.aln_samfields = {}  # dictionary to store SAM fields if present
    def __repr__(self):
        # return a string representation of the alignment
        return "{}:{}-{} -> {}:{}-{} ({} bp, strand: {})".format(self.query_sequence, self.query_sequence_start, self.query_sequence_end, self.target_sequence, self.target_sequence_start, self.target_sequence_end, self.aln_length_query, self.target_sequence_strand)
    def parse_paf_alignment(self, pafline):
        pafitems = pafline.strip().split()
        self.query_sequence = pafitems[0]
        self.query_sequence_length = int(pafitems[1])
        self.query_sequence_start = int(pafitems[2])
        self.query_sequence_end = int(pafitems[3])
        self.aln_length_query = self.query_sequence_end - self.query_sequence_start
        self.target_sequence_strand = pafitems[4]
        self.target_sequence = pafitems[5]
        self.target_sequence_length = int(pafitems[6])
        self.target_sequence_start = int(pafitems[7])
        self.target_sequence_end = int(pafitems[8])
        self.aln_length_target = self.target_sequence_end - self.target_sequence_start
        self.aln_base_matches = int(pafitems[9])
        self.aln_base_total = int(pafitems[10])
        self.aln_mapping_quality = int(pafitems[11])
        # if there are more fields, store them in the samfields dictionary
        if len(pafitems) > 12:
            for item in pafitems[12:]:
                key, value = item.split(':')[0], (item.split(':')[1], item.split(':')[2])
                self.aln_samfields[key] = value
    def is_consecutive(self, previous_alignment, max_gap=10000):
        # 1. Check if they are on the same strand
        if previous_alignment.target_sequence_strand != self.target_sequence_strand:
            return False
        # 2. Query coordinates check based on the strand
        if previous_alignment.target_sequence_strand == '+':
            # Query start must be larger in the second alignment
            if self.query_sequence_start <= previous_alignment.query_sequence_start:
                return False
            # Check if query coordinates are within max distance
            if self.query_sequence_start - previous_alignment.query_sequence_end > max_gap:
                return False
        else:  # For minus strand
            # Query start must be smaller in the second alignment
            if self.query_sequence_start >= previous_alignment.query_sequence_start:
                return False
            # Check if query coordinates are within max distance
            if previous_alignment.query_sequence_start - self.query_sequence_end > max_gap:
                return False
        # 3. Target coordinates must be larger in the second alignment and within max_gap
        if self.target_sequence_start <= previous_alignment.target_sequence_start:
            return False
        if self.target_sequence_start - previous_alignment.target_sequence_end > max_gap:
            return False
        # 4. query and target scaffolds must be the same
        if not self.query_sequence == previous_alignment.query_sequence or not self.target_sequence == previous_alignment.target_sequence:
            return False
        return True
    def merge(self, aln2):
        # merge two alignments and return a new alignment
        merged_aln = Alignment()
        merged_aln.query_sequence = self.query_sequence
        merged_aln.query_sequence_length = self.query_sequence_length
        merged_aln.query_sequence_start = min(self.query_sequence_start, aln2.query_sequence_start)
        merged_aln.query_sequence_end = max(self.query_sequence_end, aln2.query_sequence_end)
        merged_aln.query_aln_length = merged_aln.query_sequence_end - merged_aln.query_sequence_start
        merged_aln.target_sequence = self.target_sequence
        merged_aln.target_sequence_length = self.target_sequence_length
        merged_aln.target_sequence_start = min(self.target_sequence_start, aln2.target_sequence_start)
        merged_aln.target_sequence_end = max(self.target_sequence_end, aln2.target_sequence_end)
        merged_aln.target_aln_length = merged_aln.target_sequence_end - merged_aln.target_sequence_start
        merged_aln.aln_base_matches = self.aln_base_matches + aln2.aln_base_matches
        merged_aln.aln_base_total = self.aln_base_total + aln2.aln_base_total
        merged_aln.aln_mapping_quality = (self.aln_mapping_quality + aln2.aln_mapping_quality) / 2
        return merged_aln
    def flip_query_alignment(self):
        self.query_sequence_start, self.query_sequence_end = self.query_sequence_length - self.query_sequence_end, self.query_sequence_length - self.query_sequence_start
        if self.target_sequence_strand == "-":
            self.target_sequence_strand = "+"
        else:
            self.target_sequence_strand = "-"
    def refpos2querypos(self, target_pos):
        # this function returns the query position that corresponds to a target position in alignment based on the cigar string.
        # if there is no cigar string throw an error
        if 'cg' not in self.aln_samfields:
            raise ValueError("No cigar string found in alignment.")
        # define whether the cigar operations consumes query and reference as a boolean tuple
        cigarops = {
            'M': (True, True),  # match or mismatch
            'I': (True, False),  # insertion to query
            'D': (False, True),  # deletion from query
            'N': (False, True),  # skipped region from query
            'S': (True, False),  # soft clipping in query
            'H': (False, False),  # hard clipping in query
            'P': (False, False),  # padding in query
            'X': (True, True),   # mismatch
            '=': (True, True)  # match
        }
        # otherwise, start on the respective start positions for ref and query, and walk the cigar string until we hit the target position
        ref_pos = self.target_sequence_start
        # store query position as relative to allow for strand correction ater
        query_pos = 0
        for n, op in re.findall(r'(\d+)([MIDNSHPX=])', self.aln_samfields['cg']):
            n = int(n)
            if op in cigarops:
                query_consumes, target_consumes = cigarops[op]
                if target_consumes:
                    if ref_pos + n >= target_pos:
                        # then figure out how far into the operation the target position is
                        relative_pos = target_pos - ref_pos
                        if query_consumes:
                            # this means that there's an exact position that we can return
                            if self.target_sequence_strand == "+":
                                return self.query_sequence_start + query_pos + relative_pos
                            else:
                                # if the strand is negative, we need to flip the query position
                                return self.query_sequence_end - (query_pos + relative_pos)
                        else:
                            # otherwise, we'll return the latest available query position
                            if self.target_sequence_strand == "+":
                                return self.query_sequence_start + query_pos
                            else:
                                return self.query_sequence_end - query_pos
                    else:
                        # otherwise, we just move the reference position forward
                        ref_pos += n
                if query_consumes:
                    query_pos += n
                else:
                    continue
    def to_paf(self):
        # convert the alignment back to a paf line
        pafitems = [
            self.query_sequence,
            str(self.query_sequence_length),
            str(self.query_sequence_start),
            str(self.query_sequence_end),
            self.target_sequence_strand,
            self.target_sequence,
            str(self.target_sequence_length),
            str(self.target_sequence_start),
            str(self.target_sequence_end),
            str(self.aln_base_matches),
            str(self.aln_base_total),
            str(self.aln_mapping_quality)
        ]
        # add samfields if present
        for key, value in self.aln_samfields.items():
            pafitems.append(f"{key}:{value[0]}:{value[1]}")
        return "\t".join(pafitems)