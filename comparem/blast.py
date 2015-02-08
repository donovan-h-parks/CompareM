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
import multiprocessing as mp

from comparem.common import remove_extension

"""
*****************************************************************************
To do:
 - blast genome against itself to find duplicate genes.
 - consider moving over to using 'diamond blastp'
 -- need to compare blastp vs. diamond in terms of speed and results
 -- diamond is not recommended for small datasets so it may not be ideal here
 -- also it isn't as sensitive
*****************************************************************************
"""


class Blast(object):
    """Wrapper for running reciprocal blast in parallel."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self._check_for_blast()

    def _check_for_blast(self):
        """Check to see if BLAST is on the system before we try to run it."""
        try:
            subprocess.call(['blastp', '-help'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  Make sure blastp is on your system path.")
            sys.exit()

    def __producer(self, output_dir, evalue, extension, producer_queue, consumer_queue):
        """Apply reciprocal blast to a pair of genomes.

        Parameters
        ----------
        output_dir : str
            Directory to store blast results.
        evalue : float
            E-value threshold used by blast.
        extension : str
            Extension of gene files.
        producer_queue : queue
            Queue containing pairs of genomes to process.
        consumer_queue : queue
            Queue to indicate completion of reciprocal blast.
        """

        while True:
            aa_gene_fileA, aa_gene_fileB = producer_queue.get(block=True, timeout=None)
            if aa_gene_fileA == None:
                break

            genome_idA = remove_extension(aa_gene_fileA, extension)
            genome_idB = remove_extension(aa_gene_fileB, extension)

            dbA = os.path.join(output_dir, genome_idA + '.db')
            dbB = os.path.join(output_dir, genome_idB + '.db')

            output_fileAB = os.path.join(output_dir, genome_idA + '-' + genome_idB + '.blastp.tsv')
            cmd = "blastp -query %s -db %s -out %s -max_target_seqs 1 -evalue %s -outfmt '6 qseqid qlen sseqid slen length pident evalue bitscore'" % (aa_gene_fileA, dbB, output_fileAB, str(evalue))
            os.system(cmd)

            output_fileBA = os.path.join(output_dir, genome_idB + '-' + genome_idA + '.blastp.tsv')
            cmd = "blastp -query %s -db %s -out %s -max_target_seqs 1 -evalue %s -outfmt '6 qseqid qlen sseqid slen length pident evalue bitscore'" % (aa_gene_fileB, dbA, output_fileBA, str(evalue))
            os.system(cmd)

            consumer_queue.put(aa_gene_fileA)

    def __consumer(self, num_genome_pairs, consumer_queue):
        """Track completion of reciprocal blast runs.

        Parameters
        ----------
        num_genome_pairs : int
            Number of genome pairs to processed.
        consumer_queue : queue
            Queue used to indicate completion of blast runs.
        """

        processed_items = 0
        while True:
            statusStr = '    Finished processing %d of %d (%.2f%%) genome pairs.' % (processed_items, num_genome_pairs, float(processed_items) * 100 / num_genome_pairs)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            aa_gene_fileA = consumer_queue.get(block=True, timeout=None)
            if aa_gene_fileA == None:
                break

            processed_items += 1

        sys.stdout.write('\n')

    def run(self, aa_gene_files, evalue, extension, output_dir, cpus):
        """Apply reciprocal blast to all pairs of genomes in parallel.

        Parameters
        ----------
        aa_gene_files : list of str
            Amino acid fasta files to process via reciprocal blast.
        evalue : float
            E-value threshold used by blast.
        extension : str
            Extension of gene files.
        output_dir : str
            Directory to store blast results.
        cpus : int
            Number of cpus to use.
        """

        # create the blast databases in serial
        self.logger.info('  Creating blast databases.')

        for aa_gene_file in aa_gene_files:
            genome_id = remove_extension(aa_gene_file, extension)

            blast_DB = os.path.join(output_dir, genome_id + '.db')
            log_file = os.path.join(output_dir, genome_id + '.log')
            if not os.path.exists(blast_DB):
                cmd = 'makeblastdb -dbtype prot -in %s -out %s -logfile %s' % (aa_gene_file, blast_DB, log_file)
                os.system(cmd)

        # perform blast in parallel
        self.logger.info('')
        self.logger.info('  Identifying blast hits between all pairs of genomes:')

        # populate producer queue with data to process
        producer_queue = mp.Queue()
        num_pairs = 0
        for i in xrange(0, len(aa_gene_files)):
            for j in xrange(i + 1, len(aa_gene_files)):
                producer_queue.put((aa_gene_files[i], aa_gene_files[j]))
                num_pairs += 1

        for _ in range(cpus):
            producer_queue.put((None, None))

        try:
            consumer_queue = mp.Queue()
            producer_proc = [mp.Process(target=self.__producer, args=(output_dir, evalue, extension, producer_queue, consumer_queue)) for _ in range(cpus)]
            consumer_proc = mp.Process(target=self.__consumer, args=(num_pairs, consumer_queue))

            consumer_proc.start()

            for p in producer_proc:
                p.start()

            for p in producer_proc:
                p.join()

            consumer_queue.put(None)
            consumer_proc.join()
        except:
            for p in producer_proc:
                p.terminate()
            consumer_proc.terminate()
