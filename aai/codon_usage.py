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

    def run(self, geneFiles, bKeepAmbiguous):
        """Calculate codon usage over a set of genomes.

        Parameters
        ----------
        geneFiles : list
            Fasta files containing called genes.
        bKeepAmbiguous: boolean
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

        codonSet = set()
        codonsInGenomes = defaultdict(lambda: defaultdict(int))

        processedItems = 0
        for f in geneFiles:
            processedItems += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genomes.' % (processedItems, len(geneFiles), float(processedItems) * 100 / len(geneFiles))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            genomeId = ntpath.basename(f)
            genomeId = genomeId[0:genomeId.rfind('.')]
            seqs = seqIO.readFasta(f)
            for _seqId, seq in seqs.iteritems():
                for i in xrange(0, len(seq), 3):
                    codon = seq[i:i + 3].upper()
                    if bKeepAmbiguous or 'N' not in codon:
                        codonSet.add(codon)
                        codonsInGenomes[genomeId][codon] += 1

        sys.stdout.write('\n')

        return codonsInGenomes, codonSet
