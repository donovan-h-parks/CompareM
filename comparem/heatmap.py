from __future__ import division

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

from itertools import permutations

import re

class Heatmap(object):
    def __init__(self, infile, outfile='test.png', tree_file=None):
        self.tree_file = tree_file
        self.outfile = outfile
        self.genomes = set()
        self._parse_data(infile)

    def _parse_data(self, infile):
        data = {}
        with open(infile) as fp:
            # waste the first line as that's the header
            next(fp)
            for line in fp:
                fields = line.rstrip().split('\t')
                fields[0] = re.sub(r'_genes$', "", fields[0])
                fields[2] = re.sub(r'_genes$', "", fields[2])
                self.genomes.add(fields[0])
                self.genomes.add(fields[2])
                try:
                    data[fields[0]][fields[2]] = [int(fields[1]), int(fields[3]), int(fields[4]), float(fields[5])]
                except KeyError:
                    data[fields[0]] = {}
                    data[fields[0]][fields[2]] = [int(fields[1]), int(fields[3]), int(fields[4]), float(fields[5])]
                except IndexError, e:
                    print fields
                    raise e

        self.perc_ids = np.zeros([len(self.genomes), len(self.genomes)])
        self.perc_aln = np.zeros([len(self.genomes), len(self.genomes)])
        genome_to_index = {}
        for n, g in enumerate(sorted(self.genomes)):
            genome_to_index[g] = n

        for g1, g2 in permutations(self.genomes, 2):
            try:
                self.perc_ids[genome_to_index[g1], genome_to_index[g2]] = data[g1][g2][3]
                self.perc_aln[genome_to_index[g1], genome_to_index[g2]] = (data[g1][g2][2] / (data[g1][g2][0] + data[g1][g2][1] - data[g1][g2][2])) * 100
            except:
                self.perc_ids[genome_to_index[g1], genome_to_index[g2]] = data[g2][g1][3]
                self.perc_aln[genome_to_index[g1], genome_to_index[g2]] = (data[g2][g1][2] / (data[g2][g1][0] + data[g2][g1][1] - data[g2][g1][2])) * 100



    def plot(self):
        ''' Make a single heatmap using the percentage alignment and identity

            Both of the matrixes are assumed to be 2-D and will be interleaved
            into a single matrix for display.  The lower triangle of the perc_ids
            will be merged with the upper triangle of the perc_aln. The lower and
            upper triangles of this merged matrix will then be colored with
            separate colormaps to make things look pretty.

            perc_ids:   2-D numpy matrix containing the values from the write_table
                        method
            perc_aln:   2-D numpy matrix containing the values from the write_table
                        method
            names:      list of the names for each row and column in the matrix
            outfile:    path to write output image to
            tree_file:  file containing a newick tree used for ordering the rows
                        and columns of the matrix
        '''


        fig = plt.figure(figsize=(0.5*len(self.genomes), 0.5*len(self.genomes)))

        gs = gridspec.GridSpec(2, 2,
                           height_ratios=[40,1]
                           )
        ax = plt.subplot(gs[0, :])
        cba_ax = plt.subplot(gs[1,0])
        cbb_ax = plt.subplot(gs[1,1])
        #ax = fig.add_subplot(111)

        #if tree_file is not None:
        #    try:
        #        from Bio import Phylo
        #    except ImportError:
        #        print "cannot import Bio.Phylo, will not reorder matrix "\
        #              "based on tree file"
        #    else:
        #        ids_pd = pd.DataFrame(perc_ids, index=names, columns=names)
        #        aln_pd = pd.DataFrame(perc_aln, index=names, columns=names)
        #        tree = Phylo.read(tree_file, 'newick')
        #        leaves = map(str, tree.get_terminals())
        #        leaves.reverse()
        #        ids_pd = ids_pd.loc[leaves, leaves]
        #        aln_pd = aln_pd.loc[leaves, leaves]
        #        perc_ids = ids_pd.as_matrix()
        #        perc_aln = aln_pd.as_matrix()
        #        names = leaves

        merged = np.tril(self.perc_ids) + np.triu(self.perc_aln)


        mask_upper   = np.transpose(np.tri(merged.shape[0]))
        mask_lower   = np.tri(merged.shape[0])
        merged_lower = np.ma.masked_array(merged, mask=mask_lower)
        merged_upper = np.ma.masked_array(merged, mask=mask_upper)

        pa      = ax.pcolormesh(merged_lower, cmap=plt.get_cmap('Blues')) #,
                        #'sequential', 9).mpl_colormap)
        pb      = ax.pcolormesh(merged_upper, cmap=plt.get_cmap('Purples')) #,
                        #'sequential', 9).mpl_colormap)

        xticks = np.arange(0.5, merged.shape[1] + 0.5)
        ax.set_xticks(xticks)
        ax.set_xticklabels(sorted(self.genomes), rotation=45, ha='right')

        yticks = np.arange(0.5, merged.shape[1] + 0.5)
        ax.set_yticks(yticks)
        ax.set_yticklabels(sorted(self.genomes))

        spines = ['top', 'bottom', 'right', 'left', 'polar']
        for spine in spines:
            # The try/except is for polar coordinates, which only have a 'polar'
            # spine and none of the others
            try:
                ax.spines[spine].set_visible(False)
            except KeyError:
                pass

        x_pos = set(['top', 'bottom'])
        y_pos = set(['left', 'right'])
        xy_pos = [x_pos, y_pos]
        xy_ax_names = ['xaxis', 'yaxis']

        for ax_name, pos in zip(xy_ax_names, xy_pos):
            axis = ax.__dict__[ax_name]
            axis.set_ticks_position('none')

        N = merged.shape[1]
        verts = [
                [0,1]
                ]
        verts.append([0, N])
        for i in reversed(range(2, N+1)):
            verts.append([i-1,i])
            verts.append([i-1,i-1])
        verts.append([0,1])
        verts = np.array(verts)

        codes = [
                Path.MOVETO
                ]
        for i in range(len(verts)-2):
            codes.append(Path.LINETO)
        codes.append(Path.CLOSEPOLY)

        verts1 = np.copy(verts)
        verts1 = verts1[:,::-1]
        verts1[:,0] = verts1[:,0]
        verts1 = verts1

        path = Path(verts, codes)

        path2 = Path(verts1, codes)
        patch = patches.PathPatch(path, lw=1, fc='none')
        ax.add_patch(patch)
        patch2 = patches.PathPatch(path2, lw=1, fc='none')
        ax.add_patch(patch2)

        cba = plt.colorbar(pa, orientation='horizontal', cax=cbb_ax)
        cbb = plt.colorbar(pb, orientation='horizontal', cax=cba_ax)

        cba.ax.axes.tick_params(labelsize=8)
        cbb.ax.axes.tick_params(labelsize=8)
        cba.set_label('Ortholog %', fontsize=10)
        cbb.set_label('AAI %', fontsize=10)

        plt.tight_layout()
        fig.savefig(self.outfile)
