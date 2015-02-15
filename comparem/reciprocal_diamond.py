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
import sys
import subprocess
import logging

from biolib.common import remove_extension
from biolib.parallel import Parallel

"""
To Do:
-- this should be renamed to indicate it performs reciprocal blast and
    generalized to allow it to be placed in biolib
"""


class ReciprocalDiamond(object):
    """Wrapper for running reciprocal blast in parallel with diamond."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger()

        self._check_for_diamond()

        self.cpus = cpus

    def _check_for_diamond(self):
        """Check to see if diamond is on the system path."""
        try:
            subprocess.call(['diamond', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  Make sure diamond is on your system path.")
            sys.exit()

    def _producer_blast(self, genome_pair):
        """Apply reciprocal blast to a pair of genomes.

        Parameters
        ----------
        genome_pair : list
            Fasta files with genes in amino acid space for each genome.
        """
        aa_gene_fileA, aa_gene_fileB = genome_pair

        genome_idA = remove_extension(aa_gene_fileA, '.genes.faa')
        genome_idB = remove_extension(aa_gene_fileB, '.genes.faa')

        dbA = os.path.join(self.output_dir, genome_idA + '.db')
        dbB = os.path.join(self.output_dir, genome_idB + '.db')

        output_fileAB = os.path.join(self.output_dir, genome_idA + '-' + genome_idB + '.blastp.tsv')
        cmd = "diamond blastp --compress 0 -p %d -q %s -d %s -o %s -k 1 -e %s" % (self.producer_cpus, 
                                                                                    aa_gene_fileA, 
                                                                                    dbB, 
                                                                                    output_fileAB, 
                                                                                    str(self.evalue))
        os.system(cmd)

        output_fileBA = os.path.join(self.output_dir, genome_idB + '-' + genome_idA + '.blastp.tsv')
        cmd = "diamond blastp --compress 0 -p %d -q %s -d %s -o %s -k 1 -e %s" % (self.producer_cpus, 
                                                                                    aa_gene_fileB, 
                                                                                    dbA, 
                                                                                    output_fileBA, 
                                                                                    str(self.evalue))
        os.system(cmd)

        return True

    def _producer_db(self, aa_gene_file):
        """Create diamond database.

        Parameters
        ----------
        aa_gene_files : str
            Fasta file with genes in amino acid space.
        """

        genome_id = remove_extension(aa_gene_file, self.extension)

        blast_DB = os.path.join(self.output_dir, genome_id + '.db')
        cmd = 'diamond makedb -p %d --in %s -d %s' % (self.producer_cpus, aa_gene_file, blast_DB)
        os.system(cmd)

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

    def run(self, aa_gene_files, evalue, output_dir):
        """Apply reciprocal blast to all pairs of genomes in parallel.

        Parameters
        ----------
        aa_gene_files : list of str
            Amino acid fasta files to process via reciprocal blast.
        evalue : float
            E-value threshold used by blast.
        output_dir : str
            Directory to store blast results.
        """
        
        self.evalue = evalue
        self.output_dir = output_dir

        # set CPUs per producer process
        self.producer_cpus = 1
        if self.cpus > len(aa_gene_files):
            self.producer_cpus = self.cpus / aa_gene_files

        # create the blast databases in serial
        self.logger.info('  Creating diamond databases:')

        parallel = Parallel(self.cpus)
        parallel.run(self._producer_db, None, aa_gene_files, self._progress)

        # perform reciprocal blast between all genome pairs
        self.logger.info('')
        self.logger.info('  Identifying hits between all pairs of genomes:')

        genome_pairs = []
        for i in xrange(0, len(aa_gene_files)):
            for j in xrange(i + 1, len(aa_gene_files)):
                genome_pairs.append([aa_gene_files[i], aa_gene_files[j]])

        parallel.run(self._producer_blast, None, genome_pairs, self._progress)
