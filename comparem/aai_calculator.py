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
import logging
import multiprocessing as mp

from numpy import mean, std

from biolib.seq_io import SeqIO
from biolib.common import make_sure_path_exists

"""
*****************************************************************************
To do:
 - need to take into account homologs within a genome when identifying
   reciprocal best blast hits.
*****************************************************************************
"""


class AAICalculator(object):
    """Calculate AAI between all pairs of genomes."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self.shared_genes = 'shared_genes'

    def _blast_hits(self, blast_table, per_identity_threshold, per_aln_len_threshold):
        """Identify homologs among BLAST hits.

        Determines the best hit for each query sequence which
        satisfies the conditions of being a homologous gene.

        Parameters
        ----------
        blast_table : str
            Name of table containing BLAST hits.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.

        Returns
        -------
        dict
           Parameters of top hit to a homolgous gene for each query sequence.
        """

        hits = {}
        for line in open(blast_table):
            line_split = line.split('\t')

            query_seq_id = line_split[0]
            query_len = int(line_split[1])

            sub_seq_id = line_split[2]
            # subject_len = int(line_split[3])

            aln_len = int(line_split[4])
            per_ident = float(line_split[5])
            evalue = float(line_split[6])
            bitscore = float(line_split[7])

            if per_ident >= per_identity_threshold:
                per_aln_len = aln_len * 100.0 / query_len

                if per_aln_len >= per_aln_len_threshold:
                    if query_seq_id not in hits:  # take first hit passing criteria
                        hits[query_seq_id] = [sub_seq_id, per_ident, per_aln_len, evalue, bitscore]

        return hits

    def __producer(self, gene_dir, blast_dir, gene_ext, output_dir, per_identity_threshold, per_aln_len_threshold, producer_queue, consumer_queue):
        """Identify reciprocal best blast hits between pairs of genomes.

        Parameters
        ----------
        gene_dir : str
            Directory with amino acid genes in fasta format.
        blast_dir : str
            Directory with reciprocal blast between genome pairs.
        gene_ext : str
            Extension of fasta files containing genes.
        output_dir : str
            Directory to store AAI results.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        producer_queue : queue
            Queue containing pairs of genomes to process.
        consumer_queue : queue
            Queue to store completion of AAI calculations.
        """

        if gene_ext[0] != '.':
            gene_ext = '.' + gene_ext

        shared_genes_dir = os.path.join(output_dir, self.shared_genes)
        make_sure_path_exists(shared_genes_dir)

        while True:
            genome_idA, genome_idB = producer_queue.get(block=True, timeout=None)
            if genome_idA == None:
                break

            seqIO = SeqIO()

            # count number of genes in each genome
            genes_in_genomeA = seqIO.read_fasta(os.path.join(gene_dir, genome_idA + gene_ext))
            genes_in_genomeB = seqIO.read_fasta(os.path.join(gene_dir, genome_idB + gene_ext))

            # find blast hits between genome A and B
            output_fileAB = os.path.join(blast_dir, genome_idA + '-' + genome_idB + '.blastp.tsv')
            blast_hitsAB = self._blast_hits(output_fileAB, per_identity_threshold, per_aln_len_threshold)

            # find blast hits between genomes B and A
            output_fileBA = os.path.join(blast_dir, genome_idB + '-' + genome_idA + '.blastp.tsv')
            blast_hitsBA = self._blast_hits(output_fileBA, per_identity_threshold, per_aln_len_threshold)

            # find reciprocal best blast hits
            fout_seqs = open(os.path.join(shared_genes_dir, genome_idA + '-' + genome_idB + '.shared_genes.faa'), 'w')

            fout_stats = open(os.path.join(shared_genes_dir, genome_idA + '-' + genome_idB + '.rbb_hits.tsv'), 'w')
            fout_stats.write(genome_idA + '\t' + genome_idB + '\tPercent Identity\tPercent Alignment Length\te-value\tbitscore\n')

            perIdentityHits = []
            for querySeqId, stats in blast_hitsAB.iteritems():
                subSeqId, perIdent, perAlnLen, evalue, bitscore = stats
                if subSeqId in blast_hitsBA and querySeqId == blast_hitsBA[subSeqId][0]:
                    fout_stats.write('%s\t%s\t%.2f\t%.2f\t%.2g\t%.2f\n' % (querySeqId, subSeqId, perIdent, perAlnLen, evalue, bitscore))

                    # take average of percent identity in both blast directions as
                    # the results will be similar, but not identical
                    avgPerIdentity = 0.5 * (perIdent + blast_hitsBA[subSeqId][1])
                    perIdentityHits.append(avgPerIdentity)

                    # write out shared genes
                    fout_seqs.write('>' + querySeqId + '_' + genome_idA + '\n')
                    fout_seqs.write(genes_in_genomeA[querySeqId] + '\n')

                    fout_seqs.write('>' + subSeqId + '_' + genome_idB + '\n')
                    fout_seqs.write(genes_in_genomeB[subSeqId] + '\n')

            fout_seqs.close()
            fout_stats.close()

            mean_per_identity_hits = 0
            if len(perIdentityHits) > 0:
                mean_per_identity_hits = mean(perIdentityHits)

            std_per_identity_hits = 0
            if len(perIdentityHits) >= 2:
                std_per_identity_hits = std(perIdentityHits)

            consumer_queue.put((genome_idA, genome_idB, len(genes_in_genomeA), len(genes_in_genomeB), len(perIdentityHits), mean_per_identity_hits, std_per_identity_hits))

    def __consumer(self, num_genome_pairs, output_dir, consumer_queue):
        """Track completion of reciprocal best blast runs.

        Parameters
        ----------
        num_genome_pairs : int
            Number of genome pairs to processed.
        output_dir : str
            Directory to store AAI results.
        consumer_queue : queue
            Queue used to indicate results of AAI calculation.
        """

        output_file = os.path.join(output_dir, 'aai_summary.tsv')
        fout = open(output_file, 'w')
        fout.write('Genome Id A\tGenes in A\tGenome Id B\tGenes in B\tOrthologous Genes\tMean AAI\tStd AAI\n')

        processed_items = 0
        while True:
            genome_idA, genome_idB, genes_in_genomeA, genes_in_genomeB, rbb_hits, mean_per_identity_hits, std_per_identity_hits = consumer_queue.get(block=True, timeout=None)
            if genome_idA == None:
                break

            processed_items += 1
            status = '    Finished processing %d of %d (%.2f%%) genome pairs.' % (processed_items, num_genome_pairs, float(processed_items) * 100 / num_genome_pairs)
            sys.stdout.write('%s\r' % status)
            sys.stdout.flush()

            fout.write('%s\t%d\t%s\t%d\t%d\t%.2f\t%.2f\n' % (genome_idA, genes_in_genomeA, genome_idB, genes_in_genomeB, rbb_hits, mean_per_identity_hits, std_per_identity_hits))

        sys.stdout.write('\n')

        self.logger.info('')
        self.logger.info('  Summary of AAI between genomes: %s' % output_file)

        fout.close()

    def run(self, genome_ids, gene_dir, blast_dir, gene_ext, per_iden_threshold, per_aln_len_threshold, output_dir, cpus):
        """Calculate amino acid identity (AAI) between pairs of genomes.

        Parameters
        ----------
        genome_ids : list of str
            Unique ids of genomes to process.
        gene_dir : str
            Directory with amino acid genes in fasta format.
        blast_dir : str
            Directory with reciprocal blast between genome pairs.
        gene_ext : str
            Extension of fasta files containing genes.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        output_dir : str
            Directory to store AAI results.
        cpus : int
            Number of cpus to use.
        """

        self.logger.info('  Calculating amino acid identity between all pairs of genomes:')

        # populate producer queue with data to process
        producer_queue = mp.Queue()
        genome_ids.sort(key=str.lower)

        num_pairs = 0
        for i in xrange(0, len(genome_ids)):
            for j in xrange(i + 1, len(genome_ids)):
                producer_queue.put((genome_ids[i], genome_ids[j]))
                num_pairs += 1

        for _ in range(cpus):
            producer_queue.put((None, None))

        try:
            consumer_queue = mp.Queue()

            producer_proc = [mp.Process(target=self.__producer, args=(gene_dir,
                                                                        blast_dir,
                                                                        gene_ext,
                                                                        output_dir,
                                                                        per_iden_threshold,
                                                                        per_aln_len_threshold,
                                                                        producer_queue,
                                                                        consumer_queue)) for _ in range(cpus)]
            write_proc = mp.Process(target=self.__consumer, args=(num_pairs, output_dir, consumer_queue))

            write_proc.start()

            for p in producer_proc:
                p.start()

            for p in producer_proc:
                p.join()

            consumer_queue.put((None, None, None, None, None, None, None))
            write_proc.join()
        except:
            for p in producer_proc:
                p.terminate()

            write_proc.terminate()
