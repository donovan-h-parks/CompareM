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

import sys
import logging
import ntpath
from collections import defaultdict

from comparem.seq_io import SeqIO


class CodonUsage(object):
    """Calculate codon usage over a set of genomes."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def run(self, gene_files, keep_ambiguous):
        """Calculate codon usage over a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes.
        keep_ambiguous: boolean
            Keep codons with ambiguous bases.

        Returns
        -------
        dict of dict : dict[genomeId][codon]
           Codon usage of each genome.
        set
           Set with all identified codons.
        """

        self.logger.info('  Calculating codon usage for each genome.')

        seqIO = SeqIO()

        codon_set = set()
        genomes_codon_usage = defaultdict(lambda: defaultdict(int))

        processed_items = 0
        for f in gene_files:
            processed_items += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, len(gene_files), float(processed_items) * 100 / len(gene_files))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            genomeId = ntpath.basename(f)
            genomeId = genomeId[0:genomeId.rfind('.')]
            seqs = seqIO.readFasta(f)
            for _seqId, seq in seqs.iteritems():
                for i in xrange(0, len(seq), 3):
                    codon = seq[i:i + 3].upper()
                    if keep_ambiguous or 'N' not in codon:
                        codon_set.add(codon)
                        genomes_codon_usage[genomeId][codon] += 1

        sys.stdout.write('\n')

        return genomes_codon_usage, codon_set
