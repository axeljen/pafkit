# simple class for storing pairwise alignment info from PAF files, and a bunch of functions for working with these
from .genome import Genome
from .aln import Alignment
import sys
import numpy as np
import warnings

class PAF():
    def __init__(self):
        self.alignments = []
        self.target_genome = None
        self.query_genome = None
        self.sorted = None
        self.header = None
    def __repr__(self):
        # print the first 10 alignments
        string_repr = "Paf with {} alignments.\nTotal alignment length: {}.".format(len(self.alignments), self.aln_length())
        if self.target_genome:
            string_repr += "\n"
        return string_repr
    # parse paf file and populate alignments list
    def parse_paf(self, paf_file, target_genome = None, query_genome = None, overwrite=False):
        # if paf is not empty, require overwrite
        if len(self.alignments) > 0 and not overwrite:
            raise ValueError("Alignments already present. Set overwrite=True to overwrite.")
        self.alignments = []
        header = ""
        # if we're reading from standard input
        if hasattr(paf_file, 'read'):
            f = paf_file
        else:
            f = open(paf_file)
        for line in f.readlines():
            if line.startswith("@"):
                header += line
                continue
            alignment = Alignment()
            alignment.parse_paf_alignment(line)
            self.alignments.append(alignment)
        if target_genome:
            self.target_genome = Genome()
            self.target_genome.parse_genome(target_genome)
        else:
            warnings.warn("No target genome index file provided, some functions require this to work.")
        if query_genome:
            self.query_genome = Genome()
            self.query_genome.parse_genome(query_genome)
        else:
            warnings.warn("No query genome index file provided, some functions require this to work.")
    def filter_alignments(self, minlen=10000, minmapq=40, minpid=0.9, only_ref_chroms=False, only_query_chroms=False, min_target_sequence_length=0, min_query_sequence_length=0,aligned_fraction=None):
        """
        Filter out alignments based on minimum length, minimum mapping quality, and minimum percent identity.
        Additionally, we can remove alignments to sequences not present in the target/query genome index files, such that one can choose to, e.g., 
        only keep alignmnets mapping to contigs/sequences mapped to chromosomes.
        """
        # filter alignments based on user-defined criteria
        #print("Filtering alignments: minimum length: {}, minimum mapping quality: {}, minimum percent identity: {}.".format(minlen, minmapq, minpid))
        if only_ref_chroms:
            #print("Removing alignments to sequences not present in the reference genome.")
            self.alignments = [a for a in self.alignments if a.target_sequence in self.target_genome.sequences]
        if only_query_chroms:
            #print("Removing alignments from sequences not present in the query genome.")
            self.alignments = [a for a in self.alignments if a.query_sequence in self.query_genome.sequences]
        # remove alignments from/to sequences shorter than the minimum length
        #print("Removing alignments to reference sequences shorter than {} bp.".format(min_target_sequence_length))
        self.alignments = [a for a in self.alignments if a.target_sequence_length >= min_target_sequence_length and a.query_sequence_length >= min_query_sequence_length]
        # then filter based on the remaining criteria
        self.alignments = [a for a in self.alignments if a.aln_length_query >= minlen and a.aln_mapping_quality >= minmapq and a.aln_base_matches/a.aln_base_total >= minpid]
        # last remove sequences without any existing alignments from the genome objects
        #print("Number of target chromosomes before filtering: {}".format(len(self.target_genome.sequences) if self.target_genome else 0))
        #print("Number of query chromosomes before filtering: {}".format(len(self.query_genome.sequences) if self.query_genome else 0))
        if self.target_genome:
            aligned_target_sequences = set([a.target_sequence for a in self.alignments])
            #print("Aligned target sequences: {}".format(aligned_target_sequences))
            self.target_genome.sequences = {seq: length for seq, length in self.target_genome.sequences.items() if seq in aligned_target_sequences}
        if self.query_genome:
            aligned_query_sequences = set([a.query_sequence for a in self.alignments])
            #print("Aligned query sequences: {}".format(aligned_query_sequences))
            self.query_genome.sequences = {seq: length for seq, length in self.query_genome.sequences.items() if seq in aligned_query_sequences}
        if aligned_fraction:
            # this setting specifies the minimun fraction of a query sequence to be covered by alignments to be kept
            if not self.query_genome:
                raise ValueError("Query genome must be provided to filter on aligned fraction.")
            query_coverage = {seq: 0 for seq in self.query_genome.sequences.keys()}
            for aln in self.alignments:
                query_coverage[aln.query_sequence] += aln.aln_length_query
            query_coverage = {seq: cov / self.query_genome.sequences[seq] for seq, cov in query_coverage.items()}
            filtered_alignments = []
            for aln in self.alignments:
                if query_coverage[aln.query_sequence] >= aligned_fraction:
                    filtered_alignments.append(aln)
            self.alignments = filtered_alignments
            # and last remove sequences without any existing alignments from the genome objects
            aligned_query_sequences = set([a.query_sequence for a in self.alignments])
            self.query_genome.sequences = {seq: length for seq, length in self.query_genome.sequences.items() if seq in aligned_query_sequences}
            
        #print("Number of target chromosomes after filtering: {}".format(len(self.target_genome.sequences) if self.target_genome else 0))
        #print("Number of query chromosomes after filtering: {}".format(len(self.query_genome.sequences) if self.query_genome else 0))
    def sort_on_target(self, suppress_warnings=False):
        """
        Sort alignments based on target genome position, really only meaningful if a target genome index file has been provided.
        """
		# if there's a target genome index, we'll use this to sort the alignments
        if self.target_genome:
            # sort based on target genome sequences
            chromorder = {chrom: i for i, chrom in enumerate(list(self.target_genome.sequences.keys()))}
            self.alignments = sorted(self.alignments, key=lambda x: (chromorder[x.target_sequence], x.target_sequence_start))
        else:
            if not suppress_warnings:
                print("Warning: No reference genome provided, sorting alphanumerically only.")
                #   continue
            self.alignments = sorted(self.alignments, key=lambda x: (x.target_sequence, x.target_sequence_start))
        self.sorted = "target"
    def sort_on_query(self):
        """
        Sort alignments based on query genome position, really only meaningful if a query genome index file has been provided.
        """
        if self.query_genome:
            # sort based on query genome sequences
            chromorder = {chrom: i for i, chrom in enumerate(list(self.query_genome.sequences.keys()))}
            self.alignments = sorted(self.alignments, key=lambda x: (chromorder[x.query_sequence], x.query_sequence_start))
        else:
            print("Warning: No query genome provided, sorting alphanumerically only.")
            self.alignments = sorted(self.alignments, key=lambda x: (x.query_sequence, x.query_sequence_start))
        self.sorted = "query"
    def merge_consecutive(self, max_gap=0):
        """
        Merge consecutive alignments, that are less than max_gap apart on both target and query.
        """
        if self.sorted is not "target":
            raise ValueError("Alignments must be sorted before merging.")
        merged_alignments = [self.alignments[0]]
        for aln in self.alignments[1:]:
            if self.sorted is None or self.sorted == "query":
                raise ValueError("Alignments must be sorted by reference position before merging.")
            if aln.is_consecutive(merged_alignments[-1], max_gap=max_gap):
                merged_alignments[-1] = aln.merge(merged_alignments[-1])
            else:
                merged_alignments.append(aln)
        self.alignments = merged_alignments
    def get_dominating_strand(self, query_sequence):
        """ 
        Get the dominating strand for a query sequence.
        Can be used to get rid of things that look like whole-chromosome inversions, which are just caused be the query being in the opposite orientation.
        """
        query_alignments = [a for a in self.alignments if a.query_sequence == query_sequence]
        plus_strand = sum([a.aln_length_query for a in query_alignments if a.target_sequence_strand == "+"])
        minus_strand = sum([a.aln_length_query for a in query_alignments if a.target_sequence_strand == "-"])
        if plus_strand > minus_strand:
            return "+"
        else:
            return "-"
    def flip_query_alignments(self, query_sequence):
        """ 
        Based on the above function, this one can be used to swap strand of all alignments for a given query sequence.
        """
        flipped_alignments = []
        for aln in self.alignments:
            if aln.query_sequence == query_sequence:
                aln.flip_query_alignment()
            flipped_alignments.append(aln)
        self.alignments = flipped_alignments
    def get_median_target_position(self, query_chrom):
        target_positions = []
        target_lengths = []
        for aln in self.fetch_query(query_chrom):
            target_positions.append((aln.target_sequence_start + self.target_genome.cumulative_startpos[aln.target_sequence]))
            target_lengths.append(aln.aln_length_target)
        if not target_positions:
            # return infinity if no positions found
            return float('inf')
        median_pos = weighted_median(target_positions, target_lengths)
        return median_pos
    def get_min_target_position(self, query_chrom):
        previous_position = 0
        previous_sequence = None
        for aln in self.alignments:
            if not previous_sequence == aln.target_sequence:
                previous_position = 0
                previous_sequence = aln.target_sequence
            if aln.query_sequence == query_chrom:
                return (aln.target_sequence, previous_position, aln.target_sequence_start)
            previous_position = aln.target_end
        return (None,None, None)
    def get_max_target_position(self, query_chrom):
        maxpos = 0
        for i,aln in enumerate(self.alignments):
            if aln.query_sequence == query_chrom:
                maxpos = aln.target_sequence
                chrom = aln.target_sequence
                try:
                    next_aln = self.alignments[i + 1]
                    if next_aln.query_sequence == query_chrom:
                        maxpos_next = next_aln.target_sequence_end
                    else:
                        maxpos_next = self.target_genome.sequences[aln.target_sequence] 
                except IndexError:
                    maxpos_next = self.target_genome.sequences[aln.target_sequence]
        return (chrom, maxpos, maxpos_next)
    def aln_length(self):
        return sum([a.aln_length_target for a in self.alignments])
    def get_query_position_from_target(self, target_chrom, target_pos):
        """ Very (too) simplistic liftover from target to query. This one should be updated. """
        for aln in self.alignments:
            if aln.target_sequence == target_chrom and aln.target_start <= target_pos <= aln.target_end:
                # check how far from the start of the alignment our target position is
                relative_pos = (target_pos - aln.target_start) / aln.target_aln_length
                # calculate the query position
                if aln.target_sequence_strand == "+":
                    return (aln.query_sequence, round(aln.query_start + relative_pos * aln.aln_length_query, 0))
                else:
                    return (aln.query_sequence, round(aln.query_end - relative_pos * aln.aln_length_query, 0))
        return None
    # def switch_strand_sign(self,chrom):
    #     for aln in self.alignments:
    #         if aln.query_sequence == chrom:
    #             if aln.strand == "+":
    #                 aln.strand = "-"
    #             else:
    #                 aln.strand = "+"
    # def get_alignments(self,query_chrom = None, target_chrom = None):
    #     alignments = []
    #     if query_chrom:
    #         for a in self.alignments:
    #             if a.query_sequence == query_chrom:
    #                 alignments.append(a)
    #     if target_chrom:
    #         for a in self.alignments:
    #             if a.target_sequence == target_chrom:
    #                 alignments.append(a)
    #     return alignments
    def fetch(self, target_sequence = None, target_sequence_start = None, target_sequence_end = None, only_contained = False):
        if target_sequence == None and not target_sequence_start == None:
            raise ValueError("If target_sequence_start is provided, target_sequence must also be provided.")
        if target_sequence_start == None and not target_sequence_end == None:
            raise ValueError("Target end can only be provided if target_sequence_start is also provided.")
        if target_sequence_start != None and target_sequence_end == None:
            raise ValueError("If target_sequence_start is provided, target_sequence_end must also be provided.")
        # fetch alignments that overlap with the target region
        overlapping_alignments = []
        for aln in self.alignments:
            if target_sequence == None:
                overlapping_alignments.append(aln)
            else:
                if target_sequence_start == None:
                    if aln.target_sequence == target_sequence:
                        overlapping_alignments.append(aln)
                else:
                    if only_contained:
                        if aln.target_sequence == target_sequence and aln.target_sequence_start >= target_sequence_start and aln.target_sequence_end <= target_sequence_end:
                            overlapping_alignments.append(aln)
                    else:
                        if aln.target_sequence == target_sequence and aln.target_sequence_start <= target_sequence_end and aln.target_sequence_end >= target_sequence_start:
                            overlapping_alignments.append(aln)
        return overlapping_alignments
    def fetch_query(self, query_sequence = None, query_sequence_start = None, query_sequence_end = None, only_contained = False):
        if query_sequence == None and not query_sequence_start == None:
            raise ValueError("If query_sequence_start is provided, query_sequence must also be provided.")
        if query_sequence_start == None and not query_sequence_end == None:
            raise ValueError("Query end can only be provided if query_sequence_start is also provided.")
        if query_sequence_start != None and query_sequence_end == None:
            raise ValueError("If query_sequence_start is provided, query_sequence_end must also be provided.")
        # fetch alignments that overlap with the query region
        overlapping_alignments = []
        for aln in self.alignments:
            if query_sequence == None:
                overlapping_alignments.append(aln)
            else:
                if query_sequence_start == None:
                    if aln.query_sequence == query_sequence:
                        overlapping_alignments.append(aln)
                else:
                    if only_contained:
                        if aln.query_sequence == query_sequence and aln.query_start >= query_sequence_start and aln.query_end <= query_sequence_end:
                            overlapping_alignments.append(aln)
                    else:
                        if aln.query_sequence == query_sequence and aln.query_start <= query_sequence_end and aln.query_end >= query_sequence_start:
                            overlapping_alignments.append(aln)
        return overlapping_alignments
    def sort_query_chroms(self, by = 'median_target'):
        if by != 'median_target':
            raise ValueError("Only 'median_target' sorting is implemented.")
        if not self.query_genome:
            raise ValueError("Query genome must be provided to sort query chromosomes.")
        if not self.target_genome:
            raise ValueError("Target genome must be provided to sort query chromosomes.")
        # get median target positions of all query sequences
        query_positions = []
        for sequence in self.query_genome.sequences:
            median_start = self.get_median_target_position(sequence)
            if median_start is not None:
                query_positions.append({'sequence': sequence, 'median_start': median_start})
        # sort on median target position
        query_positions = sorted(query_positions, key=lambda x: x['median_start'])
        # get the new order of query sequences
        new_order = [q['sequence'] for q in query_positions]
        # sort the query genome sequences
        self.query_genome.sort_sequences(new_order)
    def add_cumulative_query_start(self, sequence_gap=0):
        # sort on target
        self.sort_on_target()
        if self.target_genome.cumulative_startpos is None or self.target_genome.cumulative_startpos == {}:
            self.target_genome.add_cumulative_startpos(sequence_gap)
        # get median startpos of query sequences
        query_positions = []
        for sequence in self.query_genome.sequences:
            median_start = self.get_median_target_position(sequence)
            if median_start is not None:
                query_positions.append({'sequence': sequence, 'median_start': median_start})
            else:
                query_positions.append({'sequence': sequence, 'median_start': float('inf')})
        # sort on median start position
        query_positions = sorted(query_positions, key=lambda x: x['median_start'])
        # add cumulative query start positions
        cumulative_pos = 0
        for q in query_positions:
            self.query_genome.cumulative_startpos[q['sequence']] = cumulative_pos
            cumulative_pos += self.query_genome.sequences[q['sequence']] + sequence_gap
            # print("Cumulative start position of query sequence {}: {}".format(q['sequence'], self.query_genome.cumulative_startpos[q['sequence']]))
    def syntenymap(self, sequence_gap=0):
        """ 
            Simplistic way of preparing a syntenymap for plotting.
        """
        # add cumulative startpositions of chromosomes
        self.add_cumulative_query_start(sequence_gap=sequence_gap)
        # make a big string
        synmap = ""
        synmap += "{tchrom}\t{tstart_original}\t{tend_original}\t{qchrom}\t{qstart_original}\t{qend_original}\t{tstart_cumulative}\t{tend_cumulative}\t{qstart_cumulative}\t{qend_cumulative}\t{strand}\n".format(tchrom="target_sequence", tstart_original="target_sequence_start", tend_original="target_sequence_end", qchrom="query_sequence", qstart_original="query_sequence_start", qend_original="query_sequence_end", tstart_cumulative="target_genome_start", tend_cumulative="target_genome_end", qstart_cumulative="query_genome_start", qend_cumulative="query_genome_end", strand="strand")
        # print alignments
        for aln in self.alignments:
            synmap += "{tchrom}\t{tstart_original}\t{tend_original}\t{qchrom}\t{qstart_original}\t{qend_original}\t{tstart_cumulative}\t{tend_cumulative}\t{qstart_cumulative}\t{qend_cumulative}\t{strand}\n".format(tchrom=aln.target_sequence, tstart_original=aln.target_sequence_start, tend_original=aln.target_sequence_end, qchrom=aln.query_sequence, qstart_original=aln.query_sequence_start, qend_original=aln.query_sequence_end, tstart_cumulative=aln.target_sequence_start + self.target_genome.cumulative_startpos[aln.target_sequence], tend_cumulative=aln.target_sequence_end + self.target_genome.cumulative_startpos[aln.target_sequence], qstart_cumulative=aln.query_sequence_start + self.query_genome.cumulative_startpos[aln.query_sequence], qend_cumulative=aln.query_sequence_end + self.query_genome.cumulative_startpos[aln.query_sequence], strand=aln.target_sequence_strand)
        return synmap
    def write_paf(self, output_file):
        # if given an open file handle, use that
        if hasattr(output_file, 'write'):
            f = output_file
            close_file = False
        else:
            f = open(output_file, 'w')
            close_file = True
        for aln in self.alignments:
            f.write(aln.to_paf() + "\n")
        if close_file:
            f.close()
        
    
# read/write functions
def read_paf(paf_file, target_genome = None, query_genome = None):
    paf = PAF()
    if paf_file == "-":
        paf.parse_paf(sys.stdin, target_genome=target_genome, query_genome=query_genome)
    else:
        paf.parse_paf(paf_file, target_genome=target_genome, query_genome=query_genome)
    return paf


def weighted_median(positions, weights):
    """Compute weighted median of positions with given weights (alignment lengths)."""
    # Sort by position
    idx = np.argsort(positions)
    positions = np.array(positions)[idx]
    weights = np.array(weights)[idx]

    # Cumulative sum of weights
    cumulative = np.cumsum(weights)
    cutoff = cumulative[-1] / 2.0

    # Find the first position where cumulative weight exceeds half total
    return positions[np.searchsorted(cumulative, cutoff)]