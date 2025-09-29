import os
# class for storing simple genome info as chrom-length pairs
class Genome():
    def __init__(self):
        self.name = None
        self.sequences = {}
        self.cumulative_startpos = {} # to be used in some plotting functions, where we want the sequences lined up one after the other
    def __repr__(self):
        return "Genome object with {} sequences.\n{}\n...\n{}".format(len(self.sequences),list(self.sequences.keys())[0],list(self.sequences.keys())[-1])
    # add sequence to genome
    def add_sequence(self, sequence, length):
        self.sequences[sequence] = length
    def add_cumulative_startpos(self,scaffold_gap):
        # remove any existing cumulative positions, in case this one is called after some filtration
        self.cumulative_startpos = {}
        cumulative_pos = 0
        for sequence, length in self.sequences.items():
            self.cumulative_startpos[sequence] = cumulative_pos
            cumulative_pos += length + scaffold_gap
    def parse_genome(self, genome_file):
        # expects a file with at least two columns: sequence name and length, tab or space delimited without header
        self.name = os.path.basename(genome_file)
        # The idea is to be able to parse a .fai file.
        with open(genome_file) as f:
            for line in f.readlines():
                scaffold = line.strip().split()[0]
                length = line.strip().split()[1]
                self.add_sequence(scaffold, int(length))
    def __str__(self):
        return "\n".join(["{}: {}".format(sequence, length) for sequence, length in self.sequences.items()])
    def get_maxpos(self):
        maxpos = 0
        for sequence, length in self.sequences.items():
            maxpos = self.cumulative_startpos[sequence] + length if self.cumulative_startpos[sequence] + length > maxpos else maxpos
        return maxpos
    def sort_sequences(self, sequence_order):
        # sequence_order is a list of sequence names in the desired order
        sorted_sequences = {}
        for seq in sequence_order:
            if seq in self.sequences:
                sorted_sequences[seq] = self.sequences[seq]
        self.sequences = sorted_sequences

        