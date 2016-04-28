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
from collections import defaultdict

from numpy import mean, std

import biolib.seq_io as seq_io
from biolib.parallel import Parallel
from biolib.common import make_sure_path_exists, remove_extension

"""
*****************************************************************************
To do:
 - need to take into account homologs within a genome when identifying
   reciprocal best blast hits.
*****************************************************************************
"""


class AAICalculator(object):
    """Calculate AAI between all pairs of genomes."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        self.shared_genes = 'shared_genes'
        self.blast_table_file = 'all_hits.tsv'
        self.gene_file = 'all_genes.faa'

    def _genome_offsets(self, blast_table):
        """Read blast table to determine byte offsets of hits for each genome.

        Parameters
        ----------
        blast_table : str
            File containing blast hits.

        Returns
        -------
        dict : d[genome_id] -> (start_pos, end_pos)
           Start and end byte offsets of hits for each genome in blast table.
        """

        offset_table = {}
        with open(blast_table, 'rb', 512 * (10 ** 6)) as f:
            cur_query_genome = None
            start_pos = 0
            end_pos = 0
            for hit in f:
                query = hit[0:hit.find('\t')]
                query_genome = query[query.rfind('~') + 1:]

                if query_genome != cur_query_genome:
                    if cur_query_genome:
                        offset_table[cur_query_genome] = (start_pos, end_pos)

                    cur_query_genome = query_genome
                    start_pos = end_pos

                end_pos += len(hit)

            offset_table[cur_query_genome] = (start_pos, end_pos)

        return offset_table

    def _valid_hits(self, blast_stream, offset_table,
                            per_identity_threshold, per_aln_len_threshold,
                            query_genome_id, subject_genome_id):
        """Identify best hits from a genome meeting the specified criteria.

        Hits from genes within query genome are identified which
        satisfy the percent identity threshold, percent alignment
        length threshold, and are to the specified subject genome.
        For each gene, the hit with the highest bitscore is identified.

        Parameters
        ----------
        blast_stream : stream
            Stream to table with blast hits.
        offset_table : d[genome_id] -> (start_pos, end_pos)
           Start and end byte offsets of hits for each genome in blast table.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        query_genome_id : str
            Unique id of genome to obtained hits for.
        subject_genome_id : str
            Unique id of genome to considered hits to.

        Returns
        -------
        dict : d[query_id] -> list with blast hit information
           Hits from query genome meeting specified criteria.
        """

        # get valid hits for genome
        hits = {}
        start_pos, end_pos = offset_table[query_genome_id]
        blast_stream.seek(start_pos)
        while blast_stream.tell() < end_pos:
            hit = blast_stream.readline().split('\t')

            perc_iden = float(hit[2])
            if perc_iden < per_identity_threshold:
                continue

            query_id = hit[0]
            query_coverage = int(hit[7]) - int(hit[6])
            per_aln_len = query_coverage * 100.0 / self.gene_lengths[query_id]
            
            if per_aln_len < per_aln_len_threshold:
                continue

            subject_id = hit[1]
            subject_genome = subject_id[subject_id.rfind('~') + 1:]
            if subject_genome != subject_genome_id:
                continue

            evalue = float(hit[10])
            bitscore = float(hit[11])

            prev_hit = hits.get(query_id, None)
            if not prev_hit:
                hits[query_id] = [subject_id, perc_iden, per_aln_len, evalue, bitscore]
            elif prev_hit[4] < bitscore:
                # for each gene, keep the hit with the highest bitscore
                hits[query_id] = [subject_id, perc_iden, per_aln_len, evalue, bitscore]

        return hits

    def _producer(self, genome_info_pairs):
        """Identify reciprocal best blast hits between pairs of genomes.

        Parameters
        ----------
        genome_info_pairs : ((genome_idA, # genes), (genome_idB, # genes))
            Identifier of genomes to process.
        """

        blast_stream = open(self.blast_table, 'rb', 32 * (10 ** 6))

        # get genome ID and number of genes in genomes to process
        genome_infoA, genome_infoB = genome_info_pairs
        genome_idA, genes_in_genomeA = genome_infoA
        genome_idB, genes_in_genomeB = genome_infoB

        # find blast hits between genome A and B, and vice versa
        hitsAB = self._valid_hits(blast_stream, self.offset_table,
                                    self.per_identity_threshold, self.per_aln_len_threshold,
                                    genome_idA, genome_idB)
        hitsBA = self._valid_hits(blast_stream, self.offset_table,
                                    self.per_identity_threshold, self.per_aln_len_threshold,
                                    genome_idB, genome_idA)

        # report reciprocal best blast hits
        fout_stats = open(os.path.join(self.shared_genes_dir, genome_idA + '-' + genome_idB + '.rbb_hits.tsv'), 'w')
        fout_stats.write(genome_idA + '\t' + genome_idB + '\tPercent Identity\tPercent Alignment Length\te-value\tbitscore\n')

        per_identity_hits = []
        for query_id, hit_stats in hitsAB.iteritems():
            subject_id, per_identA, per_aln_lenA, evalueA, bitscoreA = hit_stats
            if subject_id in hitsBA and query_id == hitsBA[subject_id][0]:
                _subject_id, per_identB, per_aln_lenB, evalueB, bitscoreB = hitsBA[subject_id]

                # take average of statistics in both blast directions as
                # the results will be similar, but not identical
                per_ident = 0.5 * (per_identA + per_identB)
                per_identity_hits.append(per_ident)

                per_aln_len = 0.5 * (per_aln_lenA + per_aln_lenB)
                evalue = 0.5 * (evalueA + evalueB)
                bitscore = 0.5 * (bitscoreA + bitscoreB)

                fout_stats.write('%s\t%s\t%.2f\t%.2f\t%.2g\t%.2f\n' % (query_id, subject_id, per_ident, per_aln_len, evalue, bitscore))

        fout_stats.close()

        mean_per_identity_hits = 0
        if len(per_identity_hits) > 0:
            mean_per_identity_hits = mean(per_identity_hits)

        std_per_identity_hits = 0
        if len(per_identity_hits) >= 2:
            std_per_identity_hits = std(per_identity_hits)

        return (genome_idA,
                    genes_in_genomeA,
                    genome_idB,
                    genes_in_genomeB,
                    len(per_identity_hits),
                    mean_per_identity_hits,
                    std_per_identity_hits)

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : tuple
            Summary statistics for a genome pair.
        consumer_data : list
            Summary statistics of amino acid identity between genome pairs.

        Returns
        -------
        consumer_data
            Summary statistics of amino acid identity between genome pairs.
        """

        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = []

        consumer_data.append(produced_data)

        return consumer_data

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

    def run(self, blast_dir, per_iden_threshold, per_aln_len_threshold, output_dir):
        """Calculate amino acid identity (AAI) between pairs of genomes.

        Parameters
        ----------
        blast_dir : str
            Directory with reciprocal blast between genome pairs.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        output_dir : str
            Directory to store AAI results.
        """

        self.blast_dir = blast_dir

        self.per_identity_threshold = per_iden_threshold
        self.per_aln_len_threshold = per_aln_len_threshold
        self.output_dir = output_dir

        shared_genes_dir = os.path.join(output_dir, self.shared_genes)
        make_sure_path_exists(shared_genes_dir)
        self.shared_genes_dir = shared_genes_dir

        # calculate length of genes and number of genes each genome
        self.logger.info('Calculating length of genes.')
        self.gene_lengths  = {}
        genes_in_genomes = defaultdict(int)
        for seq_id, seq in seq_io.read_fasta_seq(os.path.join(self.blast_dir, self.gene_file)):
            if seq[-1] == '*':
                self.gene_lengths[seq_id] = len(seq) - 1
            else:
                self.gene_lengths[seq_id] = len(seq)
                
            genome_id = seq_id[seq_id.find('~')+1:]
            genes_in_genomes[genome_id] += 1

        # get byte offset of hits from each genome
        self.logger.info('Indexing blast hits.')
        self.blast_table = os.path.join(self.blast_dir, self.blast_table_file)
        self.offset_table = self._genome_offsets(self.blast_table)

        # calculate AAI between each pair of genomes in parallel
        ng = len(genes_in_genomes)
        num_pairs = (ng*ng - ng) / 2
        self.logger.info('Calculating amino acid identity between all %d pairs of genomes:' % num_pairs)

        genome_info_pairs = []
        genome_ids = genes_in_genomes.keys()
        for i in xrange(0, len(genome_ids)):
            genome_idI = genome_ids[i]
            genome_infoI = (genome_idI, genes_in_genomes[genome_idI])
            for j in xrange(i + 1, len(genome_ids)):
                genome_idJ = genome_ids[j]
                genome_infoJ = (genome_idJ, genes_in_genomes[genome_idJ])
                genome_info_pairs.append((genome_infoI, genome_infoJ))

        if len(genome_info_pairs) == 0:
            self.logger.warning('No genome pairs identified.')
            return

        parallel = Parallel(self.cpus)
        consumer_data = parallel.run(self._producer, self._consumer, genome_info_pairs, self._progress)

        # write results for each genome pair
        aai_summay_file = os.path.join(output_dir, 'aai_summary.tsv')
        fout = open(aai_summay_file, 'w')
        fout.write('Genome Id A\tGenes in A\tGenome Id B\tGenes in B\t# orthologous genes\tMean AAI\tStd AAI\n')

        for data in consumer_data:
            fout.write('%s\t%d\t%s\t%d\t%d\t%.2f\t%.2f\n' % data)

        fout.close()

        self.logger.info('Summary of AAI between genomes: %s' % aai_summay_file)
