import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from modules import genome, paf
import numpy as np

def synteny_plot(paf, output_file, ref_color="lightgray", query_color="lightgray", aln_color="green", aln_color_inverted="red",
                 ref_color_alpha=1.0, query_color_alpha=1.0, aln_color_alpha=0.4, aln_color_inverted_alpha=0.4, ref_labels=True, query_labels=True,
                 refgenome_label=None, querygenome_label=None, reference_segments=None, query_segments=None):
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
        if ref_labels:
            ax.text(start + length_scaled / 2, y_ref + rect_height + 0.1, sequence, ha="center", fontsize=6, rotation=45, va="bottom")    
    for sequence, length in paf.query_genome.sequences.items():
        start = paf.query_genome.cumulative_startpos[sequence] 
        length_scaled = length
        rect = Rectangle((start, y_query), length_scaled, rect_height, facecolor=query_color, edgecolor="k")
        ax.add_patch(rect)
        if query_labels:
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
    # if there are segments to highlight, add these on top of the respective chromosomes
    if reference_segments and len(reference_segments) > 0:
        for seg in reference_segments:
            if seg['sequence'] in paf.target_genome.sequences:
                start = paf.target_genome.cumulative_startpos[seg['sequence']] + seg['start']
                length = seg['end'] - seg['start']
                rect = Rectangle((start, y_ref - 0.1), length, rect_height + 0.1, facecolor=seg.get('color', 'red'), edgecolor="k")
                ax.add_patch(rect)
    if query_segments and len(query_segments) > 0:
        for seg in query_segments:
            if seg['sequence'] in paf.query_genome.sequences:
                start = paf.query_genome.cumulative_startpos[seg['sequence']] + seg['start']
                length = seg['end'] - seg['start']
                rect = Rectangle((start, y_query - 0.1), length, rect_height + 0.1, facecolor=seg.get('color', 'red'), edgecolor="k")
                ax.add_patch(rect)
    # clean up
    ax.axis("off")
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    return True


def plot_chromosome(
    start, length,
    xbase=None, ybase=None,
    size=0.2, color="gray",
    direction="horizontal", radius=0.1,
    npoints=50
):
    """
    Plot a chromosome with adjustable thickness (size) and cap radius.
    
    Parameters
    ----------
    start : float
        Start coordinate (x if horizontal, y if vertical)
    length : float
        Chromosome length (excluding caps if radius>0)
    xbase, ybase : float
        Position of the baseline (one must be provided)
    size : float
        Thickness of the chromosome (height for horizontal, width for vertical)
    radius : float
        Radius of the rounded caps; 0 gives flat ends
    direction : str
        'horizontal' or 'vertical'
    """
    if xbase is None and ybase is None:
        raise ValueError("Either xbase or ybase must be provided")
    if radius < 0:
        raise ValueError("radius must be >= 0")
    # half the thickness, for centering around the baseline
    half = size / 2
    if direction == "horizontal":
        if ybase is None:
            raise ValueError("ybase must be provided for horizontal direction")
        ycenter = ybase + half
        if radius == 0:
            verts = [
                (start, ycenter - half),
                (start + length, ycenter - half),
                (start + length, ycenter + half),
                (start, ycenter + half),
                (start, ycenter - half)
            ]
            codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
        else:
            # arcs on left and right
            theta = np.linspace(np.pi/2, -np.pi/2, npoints)
            left_arc = np.column_stack([
                start + radius * (1 - np.cos(theta)),
                ycenter + half * np.sin(theta)
            ])
            right_arc = np.column_stack([
                start + length - radius + radius * np.cos(theta),
                ycenter + half * np.sin(theta)
            ])
            verts = np.concatenate([left_arc, right_arc[::-1], [left_arc[0]]])
            codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
    else:  # vertical
        if xbase is None:
            raise ValueError("xbase must be provided for vertical direction")
        xcenter = xbase + half
        if radius == 0:
            verts = [
                (xcenter - half, start),
                (xcenter + half, start),
                (xcenter + half, start + length),
                (xcenter - half, start + length),
                (xcenter - half, start)
            ]
            codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
        else:
            theta = np.linspace(0, np.pi, npoints)
            bottom_arc = np.column_stack([
                xcenter + half * np.cos(theta),
                start + radius * (1 - np.sin(theta))
            ])
            top_arc = np.column_stack([
                xcenter + half * np.cos(-theta),
                start + length - radius + radius * np.sin(theta)
            ])
            verts = np.concatenate([bottom_arc, top_arc[::-1], [bottom_arc[0]]])
            codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = PathPatch(path, facecolor=color, edgecolor='k', lw=1, zorder=1)
    return patch

def plot_genome(ax, genome, ybase=0, size=1, color="lightgray", radius=0.2):
    for seq in genome.sequences.keys():
        start = genome.cumulative_startpos[seq]
        length = genome.sequences[seq]
        print("Plotting chromosome {}: start {}, length {}".format(seq, start, length))
        chrom_patch = plot_chromosome(start, length, ybase=ybase, size=size, color=color, direction="horizontal", radius=radius)
        ax.add_patch(chrom_patch)
    return ax

def plot_alignment(ax, tstart, tend, qstart, qend, ystart, yend, color="blue"):
    """ plot a single alignment as a polygon between two y-positions """
    x1 = tstart
    x2 = tend
    x3 = qstart
    x4 = qend
    poly = Polygon([[x1, yend], [x2, yend], [x4, ystart], [x3, ystart]],
                   closed=True, facecolor=color, alpha=0.5, edgecolor=None)
    ax.add_patch(poly)
    return ax

def plot_alignments(ax, paf, ystart, yend, color="blue"):
    for aln in paf.alignments:
        # get the cumulative starting positions for the respective genomes
        tstart = aln.target_sequence_start + paf.target_genome.cumulative_startpos[aln.target_sequence]
        tend = aln.target_sequence_end + paf.target_genome.cumulative_startpos[aln.target_sequence]
        qstart = aln.query_sequence_start + paf.query_genome.cumulative_startpos[aln.query_sequence]
        qend = aln.query_sequence_end + paf.query_genome.cumulative_startpos[aln.query_sequence]
        ax = plot_alignment(ax, tstart, tend, qstart, qend, ystart, yend, color=color)
    return ax



class PafPlotter:
    """ class that stores genomes and alignments, and methods to plot them """
    def __init__(self, paflist, vertical_spacing=10, chrom_height=1):
        self.paflist = paflist  # list of PAF objects
        self.fig, self.ax = plt.subplots(dpi = 300)
        # Precompute the layout once!
        self.vertical_spacing = vertical_spacing
        self.chrom_height = chrom_height
        self._genome_tracks = self.genome_tracks()
        self._alignment_tracks = self.alignment_tracks()
        print(self._genome_tracks)

    def genome_tracks(self, labels = None, xmargin=10000000):
        ytop = (self.vertical_spacing + self.chrom_height) * (len(self.paflist))
        self.ytop = ytop
        self.ymax = ytop + self.vertical_spacing
        self.ymin = 0 - self.vertical_spacing
        self.xmin = 0 - xmargin
        self.xmax = max([paf['paf'].target_genome.get_maxpos() for paf in self.paflist] + 
                        [paf['paf'].query_genome.get_maxpos() for paf in self.paflist]) + xmargin
        if labels == None:
            labels = []
            for i in range(len(self.paflist) + 1):
                labels.append("Genome {}".format(i+1))
        genome_tracks = []
        for i, paf in enumerate(self.paflist):
            if i == 0:
                # for the first alignment, we add both the target and query genomes
                genome_tracks.append({'genome': paf['paf'].target_genome,
                                      'ymin': ytop - self.chrom_height,
                                      'ymax': ytop,
                                      'size': self.chrom_height,
                                      'label': labels[i]})
            # then add the query genomes for the remaining ones
            i += 1
            ytop -= (self.vertical_spacing + self.chrom_height)
            genome_tracks.append({'genome': paf['paf'].query_genome,
                                  'ymin': ytop - self.chrom_height,
                                  'ymax': ytop,
                                  'size': self.chrom_height,
                                  'label': labels[i]})
        print("Added {} genome tracks.".format(len(genome_tracks)))
        return genome_tracks
    def alignment_tracks(self):
        ytop = self.ytop - self.chrom_height
        alignment_tracks = []
        for paf in self.paflist:
            alignment_tracks.append({'alignments': paf['paf'].alignments,
                                     'yend': ytop,
                                     'ystart': ytop - self.vertical_spacing})
            ytop -= (self.vertical_spacing + self.chrom_height)
        return alignment_tracks
    def plot_genomes(self):
        for track in self._genome_tracks:
            plot_genome(
                self.ax,
                track['genome'],
                ybase=track['ymin'],
                size=track['size']
            )
        return self.ax
    def plot_alignments(self):
        for i, paf in enumerate(self.paflist):
            track = self._alignment_tracks[i]
            plot_alignments(
                self.ax,
                paf['paf'],
                ystart=track['ystart'],
                yend=track['yend']
            )
        return self.ax
    def set_xlim(self):
        self.ax.set_xlim(self.xmin, self.xmax)
    def set_ylim(self):
        self.ax.set_ylim(self.ymin, self.ymax)
    def save(self, output_file):
        self.ax.axis("off")
        self.fig.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close(self.fig)
        return True

# ytop = 50

# for i in range(5):
#     y = ytop - i * 10


# # # make a testplot with chromosomes

# # fig, ax = plt.subplots()
# # ax = plot_genome(ax, paf['paf'].target_genome, ybase=0, size=1, color="lightblue", radius=10)
# # ax = plot_genome(ax, paf['paf'].query_genome, ybase=2, size=1, color="lightgreen", radius=10)
# # # plot alignments
# # ax = plot_alignments(ax, paf['paf'], ystart=1, yend=2, color="orange")
# # ax.set_xlim(-1, max(paf['paf'].target_genome.get_maxpos(), paf['paf'].query_genome.get_maxpos()) + 1)
# # ax.set_ylim(-1, 3)
# # ax.axis("off")
# # #ax.relim()
# # #ax.autoscale_view()
# # plt.savefig("test_genome.png", bbox_inches='tight', dpi=300)
# # plt.close()


# # fig, ax = plt.subplots()
# # chr1 = plot_chromosome(0, 8, ybase=0, size=1, color="lightblue", direction="horizontal", radius=0.2)
# # chr2 = plot_chromosome(6, 9, xbase=0, size=1, color="lightgreen", direction="vertical", radius=0.2)
# # ax.add_patch(chr1)
# # ax.add_patch(chr2)
# # ax.set_xlim(-1, 15)
# # ax.set_ylim(-1, 15)
# # ax.axis("off")
# # plt.savefig("test_chromosomes.png", bbox_inches='tight', dpi=300)
# # plt.close()