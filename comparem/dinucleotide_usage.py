###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import logging
import ntpath
import operator
from collections import defaultdict

from biolib.seq_io import SeqIO
from biolib.parallel import Parallel

from numpy import mean, std, array, dot, transpose, sqrt


class DinucleotideUsage(object):
    """Calculate dinucleotide usage over a set of genomes.

    The dinucleotides are formed from the 3rd and succeeding
    1st position nucleotides. This has been suggested for
    identifying lateral gene transfer as this dinucleotide
    patten is minimally restricted by amino acid preference
    and codon usage:

    Hooper SD, Berg OG. 2002. Detection of genes with atypical
        nucleotide sequence in microbial genomes. J. Mol. Evol.
        54:365-75
    """

    def __init__(self, output_dir, cpus=1, keep_ambiguous=False):
        """Initialization.

        Parameters
        ----------
        output_dir : str
            Directory to store results.
        cpus : int
            Number of cpus to use.
        keep_ambiguous: boolean
            Keep codons with ambiguous bases.
        """
        self.logger = logging.getLogger()

        self.cpus = cpus
        self.output_dir = output_dir
        self.keep_ambiguous = keep_ambiguous

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def dinucleotide_usage(self, seqs, genome_id):
        """ Calculate dinucleotide usage within sequences.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        genome_id : str
            Unique id of genome used to create output file.

        Returns
        -------
        dict : d[dinucleotide] -> count
            Occurrence of each dinucleotide.
        """

        # calculate dinucleotide usage for each genome
        # and the genome as a while
        di_usage = defaultdict(lambda: defaultdict(int))
        genome_di_usage = defaultdict(int)
        di_set = set()
        n1 = defaultdict(int)
        n3 = defaultdict(int)
        for gene_id, seq in seqs.iteritems():
            for i in xrange(2, len(seq) - 1, 3):
                dinucleotide = seq[i:i + 2].upper()
                if self.keep_ambiguous or 'N' not in dinucleotide:
                    di_usage[gene_id][dinucleotide] += 1
                    genome_di_usage[dinucleotide] += 1
                    n3[seq[i:i + 1].upper()] += 1
                    n1[seq[i + 1:i + 2].upper()] += 1
                    di_set.add(dinucleotide)

        di_set_sorted = sorted(di_set)

        """
        Moving to method in paper above.
        -> need to properly calculates the z_GA values.
        """

        # calculate mean and std. dev. for each dinucleotide
        di_mean = defaultdict(float)
        di_std = defaultdict(float)
        N = sum(genome_di_usage.values())
        for a in ['A', 'C', 'G', 'T']:
            for b in ['A', 'C', 'G', 'T']:
                m = float(n3[a] * n1[b]) / N
                di_mean[a + b] = m
                var = m * (1.0 - min(n3[a], n1[b]) / N) * ((N - max(n3[a], n1[b])) / (N - 1.0))
                di_std[a + b] = sqrt(var)

        # calculate genomic dinculotide bias
        for di in di_set_sorted:
            genome_di_usage[di] = float(genome_di_usage.get(di, 0) - di_mean[di]) / di_std[di]

        print di_usage['cck10.paired_contig_7524_8']
        for dinucleotides in di_usage.values():
            for di in di_set_sorted:
                dinucleotides[di] = float(dinucleotides.get(di, 0) - di_mean[di]) / di_std[di]  # ************************************* di_mean so be the mean calculated over the individual gene!!
        print di_usage['cck10.paired_contig_7524_8']

        # calculate covariance matrix
        cov_matrix = []
        for di_j in di_set_sorted[0:-1]:
            row = []
            for di_k in di_set_sorted[0:-1]:
                s = 0
                for dinucleotides in di_usage.values():
                    s += (dinucleotides[di_j] - genome_di_usage[di_j]) * (dinucleotides[di_k] - genome_di_usage[di_k])
                s /= (len(di_usage) - 1)

                row.append(s)
            cov_matrix.append(row)
        cov_matrix = transpose(array(cov_matrix))

        # calculate Manhattan distances for each dinucleotide
        for gene_id, dinucleotides in di_usage.iteritems():
            x = []
            for di in di_set_sorted[0:-1]:
                x.append(dinucleotides[di] - genome_di_usage[di])
            x = array(x)

            T = dot(transpose(x), cov_matrix).dot(x)  # ****** NEED TO THROW IN INVERSE HERE!

            print T


        # calculate Manhattan distance for each gene
        dist = {}
        genome_sum_di = sum(genome_di_usage.values())
        for gene_id, dinucleotides in di_usage.iteritems():
            d = 0
            gene_sum_di = sum(dinucleotides.values())
            for di in di_set_sorted:
                d += abs(dinucleotides.get(di, 0) * 100.0 / gene_sum_di - genome_di_usage.get(di, 0) * 100.0 / genome_sum_di)
            dist[gene_id] = d

        # model all distances as a normal distribution
        m = mean(dist.values())
        s = std(dist.values())

        # calculate standard deviations from the mean
        std_mean_dict = {}
        for gene_id, d in dist.iteritems():
            std_mean_dict[gene_id] = (d - m) / s

        gene_ids_sorted = sorted(std_mean_dict.items(), key=operator.itemgetter(1), reverse=True)

        # report dinucleotide usage of each gene
        output_file = os.path.join(self.output_dir, genome_id + '.di_usage.tsv')
        fout = open(output_file, 'w')

        fout.write('Gene Id\tSeq. length (bp)\t# dinucleotides\tManhattan distance\tDeviations from mean')
        for di in di_set_sorted:
            fout.write('\t' + di)
        fout.write('\n')

        fout.write('%s\t%d\t%d' % ('<complete genome>', sum([len(x) for x in seqs.values()]), genome_sum_di))
        fout.write('\t%.1f\t%.1f' % (0, 0))
        for di in di_set_sorted:
            fout.write('\t%.2f' % (genome_di_usage.get(di, 0) * 100.0 / genome_sum_di))
        fout.write('\n')

        for gene_id, std_mean in gene_ids_sorted:
            dinucleotides = di_usage[gene_id]
            sum_di = sum(dinucleotides.values())
            fout.write('%s\t%d\t%d' % (gene_id, len(seqs[gene_id]), sum_di))
            fout.write('\t%.2f\t%.2f' % (dist[gene_id], std_mean))

            for di in di_set_sorted:
                fout.write('\t%.2f' % (dinucleotides.get(di, 0) * 100.0 / sum_di))
            fout.write('\n')
        fout.close()

    def _producer(self, gene_file):
        """Calculates codon usage of a genome.

        This function is intended to be used as a producer
        within a producer/consumer multiprocessing framework.
        It calculates the codon usage for a single genome
        and returns the results for consumption by the
        consumer function.

        Parameters
        ----------
        gene_file : str
            Fasta file containing amino acid sequences.

        Returns
        -------
        str
           Unique identifier of genome.
        dict : d[codon] -> count
            Occurrence of each codon.
        dict : d[codon] -> length
            Average length of genes for a given stop codon.
        """

        genome_id = ntpath.basename(gene_file)
        genome_id = genome_id.replace('.genes.fna', '')
        genome_id = os.path.splitext(genome_id)[0]

        seq_io = SeqIO()
        seqs = seq_io.read_fasta(gene_file)
        self.dinucleotide_usage(seqs, genome_id)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of genomes processed.
        total_items : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return '    Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self, gene_files):
        """Calculate codon usage over a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes in nucleotide space.

        Returns
        -------
        dict of dict : d[genome_id][codon] -> count
           Codon usage of each genome.
        set
           Set with all identified codons.
        dict of dict : d[genome_id][codon] -> length
            Mean length of genes for each stop codon.
        """

        self.logger.info('  Calculating codon usage for each genome.')

        parallel = Parallel(self.cpus)
        parallel.run(self._producer, None, gene_files[0:1], self._progress)

