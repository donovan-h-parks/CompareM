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
import tempfile

from biolib.common import concatenate_files
from biolib.external.diamond import Diamond

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

        self.cpus = cpus

    def run(self, aa_gene_files, evalue, per_identity, output_dir):
        """Apply reciprocal blast to all pairs of genomes in parallel.

        Parameters
        ----------
        aa_gene_files : list of str
            Amino acid fasta files to process via reciprocal blast.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        output_dir : str
            Directory to store blast results.
        """

        # concatenate all gene files and create a single diamond database
        self.logger.info('  Creating diamond database (be patient!).')
        gene_file = os.path.join(output_dir, 'all_genes.faa')
        concatenate_files(aa_gene_files, gene_file)
        diamond_db = os.path.join(output_dir, 'all_genes')

        diamond = Diamond(self.cpus)
        diamond.make_database(gene_file, diamond_db)

        # blast all genes against the database
        self.logger.info('')
        self.logger.info('  Identifying hits between all pairs of genomes (be patient!).')
        hits_daa_file = os.path.join(output_dir, 'all_hits')
        diamond.blastp(gene_file, diamond_db, evalue, per_identity, len(aa_gene_files) * 10, hits_daa_file)

        # create flat hits table
        self.logger.info('  Creating table with hits.')
        hits_table_file = os.path.join(output_dir, 'all_hits.tsv')
        diamond.view(hits_daa_file + '.daa', hits_table_file)
