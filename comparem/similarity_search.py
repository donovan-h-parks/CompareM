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
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import subprocess
import logging
import tempfile

import biolib.seq_io as seq_io
from biolib.common import concatenate_files, remove_extension, make_sure_path_exists
from biolib.parallel import Parallel
from biolib.external.diamond import Diamond
from biolib.external.blast import Blast


class SimilaritySearch(object):
    """Performs similarity search of gene sequences."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus
        
    def _prefix_gene_identifiers(self, gene_files, keep_headers, gene_out_dir):
        """Prefix all gene IDs with genome IDs: <genome_id>~<gene_id>.
        
        Parameters
        ----------
        gene_files : list of str
            Genes in fasta files to modify.
        keep_headers : boolean
            If True, indicates FASTA headers already have the format <genome_id>~<gene_id>.
        gene_out_dir : str
            Directory to store modified gene files.
        Returns
        -------
        list
            List of strings indicating files with modified gene IDs.
        """
        
        modified_gene_files = []
        for gf in gene_files:
            if keep_headers:
                modified_gene_files.append(gf)
            else:                
                genome_id = remove_extension(gf)

                aa_file = os.path.join(gene_out_dir, genome_id + '.faa')
                fout = open(aa_file, 'w')
                for seq_id, seq, annotation in seq_io.read_fasta_seq(gf, keep_annotation=True):
                    fout.write('>' + genome_id + '~' + seq_id  + ' ' + annotation + '\n')
                    fout.write(seq + '\n')
                fout.close()

                modified_gene_files.append(aa_file)
                
        return modified_gene_files
        
    def _run_blastp(self, query_gene_file, 
                            target_gene_file, 
                            evalue, 
                            per_identity, 
                            per_aln_len,
                            max_hits,
                            tmp_dir,
                            output_dir):
        """Perform similarity search of query genes against target genes.

        Parameters
        ----------
        query_gene_file : str
            File with all query sequences.
        target_gene_files : str
            File with all subject sequences
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        max_hits : int
            Maximum number of hits to report per query sequences.
        tmp_dir : str
            Directory to store temporary files.
        output_dir : str
            Directory to store blast results.
        """
        
        # concatenate all gene files and create a single diamond database
        self.logger.info('Creating diamond database (be patient!).')
        
        blast = Blast(self.cpus, silent=True)
        blast.create_blastp_db(target_gene_file)
        
        # create temporary hits table
        if tmp_dir:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            if os.path.isdir('/dev/shm'):
                tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir='/dev/shm', delete=False)
            else:
                tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_hits_table.close()

        # blast all genes against the database
        self.logger.info('Identifying hits between query and target genomes (be patient!).')
        hits_daa_file = os.path.join(output_dir, 'hits')
        blast.blastp(query_gene_file, target_gene_file, tmp_hits_table.name, evalue, max_hits, task='blastp-fast')
        
        # sort hit table
        self.logger.info('Sorting table with hits (be patient!).')
        hits_table_file = os.path.join(output_dir, 'hits_sorted.tsv')
        os.system("LC_ALL=C sed -i 's/~/\t/g' %s" % tmp_hits_table.name)
        os.system("LC_ALL=C sort --parallel=8 -o %s -k1,1 -k3,3 %s" % (tmp_hits_table.name, tmp_hits_table.name))
        os.system('mv %s %s' % (tmp_hits_table.name, hits_table_file))
        
    def _run_diamond(self, query_gene_file, 
                            target_gene_file, 
                            evalue, 
                            per_identity, 
                            per_aln_len,
                            max_hits,
                            high_mem,
                            tmp_dir,
                            output_dir):
        """Perform similarity search of query genes against target genes.

        Parameters
        ----------
        query_gene_file : str
            File with all query sequences.
        target_gene_files : str
            File with all subject sequences
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        max_hits : int
            Maximum number of hits to report per query sequences.
        tmp_dir : str
            Directory to store temporary files.
        output_dir : str
            Directory to store blast results.
        """
        
        # concatenate all gene files and create a single diamond database
        self.logger.info('Creating diamond database (be patient!).')
        
        diamond_db = os.path.join(output_dir, 'target_genes')

        diamond = Diamond(self.cpus)
        if high_mem:
            diamond.make_database(target_gene_file, diamond_db, block_size=8)
        else:
            diamond.make_database(target_gene_file, diamond_db)

        # blast all genes against the database
        self.logger.info('Identifying hits between query and target genomes (be patient!).')
        hits_daa_file = os.path.join(output_dir, 'hits')
        
        if high_mem:
            diamond.blastp(query_gene_file, diamond_db, evalue, per_identity, per_aln_len, max_hits, hits_daa_file, tmp_dir, chunk_size=1)
        else:
            diamond.blastp(query_gene_file, diamond_db, evalue, per_identity, per_aln_len, max_hits, hits_daa_file, tmp_dir)

        # create flat hits table
        if tmp_dir:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            if os.path.isdir('/dev/shm'):
                tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir='/dev/shm', delete=False)
            else:
                tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_hits_table.close()
                
        self.logger.info('Creating table with hits.')
        diamond.view(hits_daa_file + '.daa', tmp_hits_table.name)
        
        # sort hit table
        self.logger.info('Sorting table with hits (be patient!).')
        hits_table_file = os.path.join(output_dir, 'hits_sorted.tsv')
        os.system("LC_ALL=C sed -i 's/~/\t/g' %s" % tmp_hits_table.name)
        os.system("LC_ALL=C sort --parallel=8 -o %s -k1,1 -k3,3 %s" % (tmp_hits_table.name, tmp_hits_table.name))
        os.system('mv %s %s' % (tmp_hits_table.name, hits_table_file))

    def run(self, query_gene_files, 
                    target_gene_files,
                    evalue, 
                    per_identity, 
                    per_aln_len,
                    high_mem,
                    tmp_dir,
                    blastp,
                    keep_headers,
                    output_dir):
        """Perform similarity search of query genes against target genes.

        Parameters
        ----------
        query_gene_files : list of str
            Query genes in fasta files to process.
        target_gene_files : list of str
            Query genes in fasta files to process.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        tmp_dir : str
            Directory to store temporary files.
        blastp : boolean
            If True blasp-fast is used instead of DIAMOND.
        keep_headers : boolean
            If True, indicates FASTA headers already have the format <genome_id>~<gene_id>.
        output_dir : str
            Directory to store blast results.
        """
           
        # modify gene ids to include genome ids in order to ensure
        # all gene identifiers are unique across the set of genomes
        self.logger.info('Appending genome identifiers to all gene identifiers.')
        query_out_dir = os.path.join(output_dir, 'query_genes')
        make_sure_path_exists(query_out_dir)
        modified_query_gene_files = self._prefix_gene_identifiers(query_gene_files, 
                                                                    keep_headers, 
                                                                    query_out_dir)
        
        query_gene_file = os.path.join(output_dir, 'query_genes.faa')
        concatenate_files(modified_query_gene_files, query_gene_file)
            
        if query_gene_files == target_gene_files:
            target_gene_file = query_gene_file
        else:
            target_out_dir = os.path.join(output_dir, 'target_genes')
            make_sure_path_exists(target_out_dir)
            modified_target_gene_files = self._prefix_gene_identifiers(target_gene_files, 
                                                                        keep_headers, 
                                                                        target_out_dir)

            target_gene_file = os.path.join(output_dir, 'target_genes.faa')
            concatenate_files(modified_target_gene_files, target_gene_file)

        if blastp:
            self._run_blastp(query_gene_file, 
                                target_gene_file, 
                                evalue, 
                                per_identity, 
                                per_aln_len,
                                len(target_gene_files) * 10,
                                tmp_dir,
                                output_dir)
        else:
            self._run_diamond(query_gene_file, 
                                target_gene_file, 
                                evalue, 
                                per_identity, 
                                per_aln_len,
                                len(target_gene_files) * 10,
                                high_mem,
                                tmp_dir,
                                output_dir)
