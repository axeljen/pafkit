""" 
Axel Jensen, 2025
Class to identify scaffold/chromosome breaks in a query assembly based on alignments to a reference genome. 

If the alignment consists of chromosome level genome assemblies, these breaks should in principle represent 
chromosomal fission events (or fusions in the reference).

"""

class FissionFinder():
    def __init__(self):
        self.left_side_alignments = []
        self.right_side_alignments = []
        self.target_sequence = None
        self.query_sequence = None
        self.target_min = None
        self.target_max = None
        self.fission_points = []
    def get_dominating_query_left(self):
        chrom = None
        fraction = 0
        for seq in set([aln.query_sequence for aln in self.left_side_alignments]):
            seq_size = sum([aln.aln_length_query for aln in self.left_side_alignments if aln.query_sequence == seq])
            total_size = sum([aln.aln_length_query for aln in self.left_side_alignments])
            if seq_size / total_size > fraction:
                fraction = seq_size / total_size
                chrom = seq
        return chrom, fraction
    def get_dominating_query_right(self):
        chrom = None
        fraction = 0
        for seq in set([aln.query_sequence for aln in self.right_side_alignments]):
            seq_size = sum([aln.aln_length_query for aln in self.right_side_alignments if aln.query_sequence == seq])
            total_size = sum([aln.aln_length_query for aln in self.right_side_alignments])
            if seq_size / total_size > fraction:
                fraction = seq_size / total_size
                chrom = seq
        return chrom, fraction
    def get_dominating_target_left(self):
        chrom = None
        fraction = 0
        for seq in set([aln.target_sequence for aln in self.left_side_alignments + self.right_side_alignments]):
            seq_size = sum([aln.aln_length_target for aln in self.left_side_alignments + self.right_side_alignments if aln.target_sequence == seq])
            total_size = sum([aln.aln_length_target for aln in self.left_side_alignments + self.right_side_alignments])
            if seq_size / total_size > fraction:
                fraction = seq_size / total_size
                chrom = seq
        return chrom, fraction
    def get_dominating_target_right(self):
        chrom = None
        fraction = 0
        for seq in set([aln.target_sequence for aln in self.right_side_alignments]):
            seq_size = sum([aln.aln_length_target for aln in self.right_side_alignments if aln.target_sequence == seq])
            total_size = sum([aln.aln_length_target for aln in self.right_side_alignments])
            if seq_size / total_size > fraction:
                fraction = seq_size / total_size
                chrom = seq
        return chrom, fraction
    def check_fission(self, min_margin=1000000, min_fraction = 0.8):
        left_size = sum([aln.aln_length_target for aln in self.left_side_alignments])
        right_size = sum([aln.aln_length_target for aln in self.right_side_alignments])
        if left_size < min_margin or right_size < min_margin:
            return False
        target_chrom_left = None
        left_fraction_target = 0
        for seq in set([aln.target_sequence for aln in self.left_side_alignments]):
            seq_size = sum([aln.aln_length_target for aln in self.left_side_alignments if aln.target_sequence == seq])
            if seq_size / left_size > left_fraction_target:
                left_fraction_target = seq_size / left_size
                target_chrom_left = seq
        target_chrom_right = None
        right_fraction_target = 0
        for seq in set([aln.target_sequence for aln in self.right_side_alignments]):
            seq_size = sum([aln.aln_length_target for aln in self.right_side_alignments if aln.target_sequence == seq])
            if seq_size / right_size > right_fraction_target:
                right_fraction_target = seq_size / right_size
                target_chrom_right = seq
        if target_chrom_left != target_chrom_right:
            return False
        query_chrom_left = None
        left_fraction_query = 0
        for seq in set([aln.query_sequence for aln in self.left_side_alignments]):
            seq_size = sum([aln.aln_length_query for aln in self.left_side_alignments if aln.query_sequence == seq])
            if seq_size / left_size > left_fraction_query:
                left_fraction_query = seq_size / left_size
                query_chrom_left = seq
        query_chrom_right = None
        right_fraction_query = 0
        for seq in set([aln.query_sequence for aln in self.right_side_alignments]):
            seq_size = sum([aln.aln_length_query for aln in self.right_side_alignments if aln.query_sequence == seq])
            if seq_size / right_size > right_fraction_query:
                right_fraction_query = seq_size / right_size
                query_chrom_right = seq
        if query_chrom_left == query_chrom_right:
            return False
        if left_fraction_target < min_fraction or right_fraction_target < min_fraction:
            return False
        if left_fraction_query < min_fraction or right_fraction_query < min_fraction:
            return False
        self.fission_points.append({'reference_chrom': target_chrom_left, 'query_chrom_left': query_chrom_left, 'query_chrom_right': query_chrom_right, 'target_position': (self.target_min, self.target_max)})
        return True
    def reset(self):
        self.left_side_alignments = self.right_side_alignments
        self.right_side_alignments = []
        self.target_min = self.target_max
        self.target_max = None
        self.query_sequence = self.get_dominating_query_left()[0]
    def find_fissions(self, alignments, min_margin=1000000, min_fraction = 0.8):
        for aln in alignments:
            print("Checking aln {}:{}-{} ({}), {}:{}-{} ({})".format(aln.target_sequence, aln.target_sequence_start, aln.target_sequence_end, aln.target_sequence_strand, aln.query_sequence, aln.query_sequence_start, aln.query_sequence_end, aln.query_sequence_strand))
            if len(self.left_side_alignments) == 0:
                self.left_side_alignments.append(aln)
                self.target_sequence = aln.target_sequence
                self.query_sequence = aln.query_sequence
                self.target_min = aln.target_sequence_end
            else:
                if aln.target_sequence == self.target_sequence and aln.query_sequence == self.query_sequence:
                    self.left_side_alignments.append(aln)
                    self.target_min = aln.target_sequence_end
                else:
                    if aln.query_sequence != self.query_sequence:
                        self.right_side_alignments.append(aln)
                        self.target_max = aln.target_sequence_start
                        if self.check_fission(min_margin=min_margin, min_fraction=min_fraction):
                            print("Fission detected at {}:{} between query chromosomes {} and {}".format(self.target_sequence, (self.target_min + self.target_max)//2, self.get_dominating_query_left()[0], self.get_dominating_query_right()[0]))
                        self.reset()
                    else:
                        self.left_side_alignments.append(aln)
                        self.target_min = aln.target_sequence_end
        # check for fission at the end
        if len(self.right_side_alignments) > 0:
            if self.check_fission(min_margin=min_margin, min_fraction=min_fraction):
                print("Fission detected at {}:{} between query chromosomes {} and {}".format(self.target_sequence, (self.target_min + self.target_max)//2, self.get_dominating_query_left()[0], self.get_dominating_query_right()[0]))
        return self.fission_points
                        
            
        