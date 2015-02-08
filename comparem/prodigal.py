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
import ntpath
import logging
import tempfile
import shutil
import multiprocessing as mp

from comparem.seq_io import SeqIO
from comparem.common import check_file_exists

import numpy as np


class Prodigal(object):
    """Wrapper for running Prodigal in parallel."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self._check_for_prodigal()

    def _check_for_prodigal(self):
        """Check to see if Prodigal is on the system before we try to run it."""
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.info("  Make sure prodigal is on your system path.")
            sys.exit()

    def __producer(self, called_genes, output_dir, producer_queue, consumer_queue):
        """Apply reciprocal blast to a pair of genomes.

        Parameters
        ----------
        called_genes : boolean
            Flag indicating genes are already called.
        output_dir : str
            Directory to store called genes.
        producer_queue : queue
            Queue containing fasta files to process.
        consumer_queue : queue
            Queue to indicate completion of prodigal.
        """

        while True:
            genome_file = producer_queue.get(block=True, timeout=None)
            if genome_file == None:
                break

            genome_id = ntpath.basename(genome_file)
            genome_id = genome_id[0:genome_id.rfind('.')]

            aa_gene_file = os.path.join(output_dir, genome_id + '.genes.faa')
            nt_gene_file = os.path.join(output_dir, genome_id + '.genes.fna')
            gff_file = os.path.join(output_dir, genome_id + '.gff')

            if called_genes:
                os.system('ln -s %s %s' % (os.path.abspath(genome_file), aa_gene_file))
            else:
                tmp_dir = tempfile.mkdtemp()

                seqIO = SeqIO()
                seqs = seqIO.read_fasta(genome_file)

                # determine number of bases
                total_bases = 0
                for seq in seqs.values():
                    total_bases += len(seq)

                # call genes under different translation tables
                table_coding_density = {}
                for translation_table in [4, 11]:
                    os.makedirs(os.path.join(tmp_dir, str(translation_table)))
                    aa_gene_file_tmp = os.path.join(tmp_dir, str(translation_table), genome_id + '.genes.faa')
                    nt_gene_file_tmp = os.path.join(tmp_dir, str(translation_table), genome_id + '.genes.fna')
                    gff_file_tmp = os.path.join(tmp_dir, str(translation_table), genome_id + '.gff')

                    # check if there is sufficient bases to calculate prodigal parameters
                    if total_bases < 100000:
                        proc_str = 'meta'  # use best precalculated parameters
                    else:
                        proc_str = 'single'  # estimate parameters from data

                    cmd = 'prodigal -p %s -q -f gff -g %d -a %s -d %s -i %s > %s 2> /dev/null' % (proc_str,
                                                                                                  translation_table,
                                                                                                  aa_gene_file_tmp,
                                                                                                  nt_gene_file_tmp,
                                                                                                  genome_file,
                                                                                                  gff_file_tmp)
                    os.system(cmd)

                    # determine coding density
                    prodigalParser = ProdigalGeneFeatureParser(gff_file_tmp)

                    codingBases = 0
                    for seq_id, seq in seqs.iteritems():
                        codingBases += prodigalParser.coding_bases(seq_id)

                    codingDensity = float(codingBases) / total_bases
                    table_coding_density[translation_table] = codingDensity

                # determine best translation table
                best_translation_table = 11
                if (table_coding_density[4] - table_coding_density[11] > 0.05) and table_coding_density[4] > 0.7:
                    best_translation_table = 4

                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table), genome_id + '.genes.faa'), aa_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table), genome_id + '.genes.fna'), nt_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table), genome_id + '.gff'), gff_file)

                # clean up temporary files
                shutil.rmtree(tmp_dir)

            consumer_queue.put(genome_id)

    def __consumer(self, num_genomes, consumer_queue):
        """Track completion of reciprocal blast runs.

        Parameters
        ----------
        num_genomes : int
            Number of genomes to processed.
        consumer_queue : queue
            Queue used to indicate completion of blast runs.
        """

        processed_items = 0
        while True:
            status = '    Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, num_genomes, float(processed_items) * 100 / num_genomes)
            sys.stdout.write('%s\r' % status)
            sys.stdout.flush()

            genome_id = consumer_queue.get(block=True, timeout=None)
            if genome_id == None:
                break

            processed_items += 1

        sys.stdout.write('\n')

    def run(self, genome_files, called_genes, output_dir, cpus):
        """Call genes with Prodigal.

        Parameters
        ----------
        genome_files : list of str
            Nucleotide fasta files to call genes on.
        called_genes : boolean
            Flag indicating genes are already called.
        output_dir : str
            Directory to store called genes.
        cpus : int
            Number of cpus to use.
        """

        self.logger.info('  Identifying genes within genomes:')

        # populate producer queue with data to process
        producer_queue = mp.Queue()
        for gf in genome_files:
            producer_queue.put(gf)

        for _ in range(cpus):
            producer_queue.put(None)

        try:
            consumer_queue = mp.Queue()
            producer_proc = [mp.Process(target=self.__producer, args=(called_genes, output_dir, producer_queue, consumer_queue)) for _ in range(cpus)]
            consumer_proc = mp.Process(target=self.__consumer, args=(len(genome_files), consumer_queue))

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


class ProdigalGeneFeatureParser():
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename):
        """Initialization.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        check_file_exists(filename)

        self.genes = {}
        self.last_coding_base = {}

        self.__parseGFF(filename)

        self.coding_base_masks = {}
        for seq_id in self.genes:
            self.coding_base_masks[seq_id] = self.__build_coding_base_mask(seq_id)

    def __parseGFF(self, filename):
        """Parse genes from GFF file.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        bGetTranslationTable = True
        for line in open(filename):
            if bGetTranslationTable and line.startswith('# Model Data'):
                self.translationTable = line.split(';')[4]
                self.translationTable = int(self.translationTable[self.translationTable.find('=') + 1:])
                bGetTranslationTable = False

            if line[0] == '#':
                continue

            line_split = line.split('\t')
            seq_id = line_split[0]
            if seq_id not in self.genes:
                geneCounter = 0
                self.genes[seq_id] = {}
                self.last_coding_base[seq_id] = 0

            geneId = seq_id + '_' + str(geneCounter)
            geneCounter += 1

            start = int(line_split[3])
            end = int(line_split[4])

            self.genes[seq_id][geneId] = [start, end]
            self.last_coding_base[seq_id] = max(self.last_coding_base[seq_id], end)

    def __build_coding_base_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        """

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        coding_base_mask = np.zeros(self.last_coding_base[seq_id])
        for pos in self.genes[seq_id].values():
            coding_base_mask[pos[0]:pos[1] + 1] = 1

        return coding_base_mask

    def coding_bases(self, seq_id, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end).

        To process the entire sequence set start to 0, and
        end to None.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        start : int
            Start calculation at this position in sequence.
        end : int
            End calculation just before this position in the sequence.
        """

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end == None:
            end = self.last_coding_base[seq_id]

        return np.sum(self.coding_base_masks[seq_id][start:end])
