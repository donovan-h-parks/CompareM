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
import platform
import shutil
import itertools

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

    def _producer(self, cmd):
        """Run command."""
    
        process = subprocess.Popen(
                        ["bash", "-c", cmd],
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        stdin = None)
        stdout, stderr = process.communicate(None)

        if process.returncode != 0:
            self.logger.error('Failed to execute:')
            self.logger.error(cmd)
            self.logger.error('Program returned code: %d' % process.returncode)
            
            if stderr.strip():
                self.logger.error('Stderr:')
                self.logger.error(stderr.strip().decode('utf-8'))
            
            if stdout.strip():
                self.logger.error('Stdout:')
                self.logger.error(stdout.strip().decode('utf-8'))
                
            sys.exit(-1)
            
        return True
            
    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes."""

        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = []

        consumer_data.append(produced_data)

        return consumer_data

    def _progress_genomes(self, processed_items, total_items):
        """Report progress of consumer processes."""

        return '  Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, 
                                                                        total_items, 
                                                                        float(processed_items) * 100 / total_items)
                                                                        
    def _progress_comparisons(self, processed_items, total_items):
        """Report progress of consumer processes."""

        return '  Finished processing %d of %d (%.2f%%) comparisons.' % (processed_items, 
                                                                        total_items, 
                                                                        float(processed_items) * 100 / total_items)
                                                                    
    def _create_db_cmd(self, similarity_method, seq_file, db_prefix):
        """Get command for creating sequence database."""
        
        if similarity_method == 'AAI_DIAMOND':
            cmd = 'diamond makedb --quiet -p %d --in %s -d %s' % (1, seq_file, db_prefix + '.dmnd')
        elif similarity_method == 'AAI_BLASTP-FAST':
            cmd = 'makeblastdb -dbtype prot -in %s -out %s  > /dev/null' % (seq_file, db_prefix + '.db')
        elif similarity_method in ['ANI_BLASTN', 'ANI_DC-MEGABLAST', 'ANI_MEGABLAST']:
            cmd = 'makeblastdb -dbtype nucl -in %s -out %s  > /dev/null' % (seq_file, db_prefix + '.db')
        else:
            self.logger.error('Unknown similarity search method: %s' % similarity_method)
            sys.exit(-1)
        
        return cmd
        
    def _search_cmd(self, similarity_method,
                            query_file,
                            db_prefix,
                            hit_table,
                            evalue,
                            per_identity, 
                            per_aln_len, 
                            sensitive):
        """Get command for searching sequence database."""

        if similarity_method == 'AAI_DIAMOND':
            args = ''
            if sensitive:
                args += ' --sensitive'

            cmd = "diamond blastp --quiet -p %d -q %s -d %s -e %g --id %f --query-cover %f -k %d -o %s -f %s %s" % (1,
                                                                                                                    query_file,
                                                                                                                    db_prefix,
                                                                                                                    evalue,
                                                                                                                    per_identity,
                                                                                                                    per_aln_len,
                                                                                                                    1,
                                                                                                                    hit_table,
                                                                                                                    '6',
                                                                                                                    args)
        elif similarity_method == 'AAI_BLASTP-FAST':
            cmd = "blastp -task blastp-fast"
            cmd += " -num_threads %d -query %s -db %s -out %s -evalue %g" % (1, query_file, db_prefix + '.db', hit_table, evalue)
            cmd += " -max_target_seqs 1 -max_hsps 1"
            cmd += " -outfmt '%s'" % '6 qseqid qlen sseqid slen length mismatch gaps pident bitscore evalue'
        elif similarity_method in ['ANI_BLASTN', 'ANI_DC-MEGABLAST', 'ANI_MEGABLAST']:
            if similarity_method == 'ANI_BLASTN':
                cmd = "blastn -task blastn"
            elif similarity_method == 'ANI_DC-MEGABLAST':
                cmd = "blastn -task dc-megablast"
            elif similarity_method == 'ANI_MEGABLAST':
                cmd = "blastn -task megablast"
            
            cmd += " -xdrop_gap_final 150 -dust no"
            cmd += " -num_threads %d -query %s -db %s -out %s -evalue %g" % (1, query_file, db_prefix + '.db', hit_table, evalue)
            cmd += " -max_target_seqs 1 -max_hsps 1"
            cmd += " -outfmt '%s'" % '6 qseqid qlen sseqid slen length mismatch gaps pident bitscore evalue'
        else:
            self.logger.error('Unknown similarity search method: %s' % similarity_method)
            sys.exit(-1)
        
        return cmd
        
    def run(self, 
                similarity_method,
                query_files,
                target_files, 
                evalue, 
                per_identity, 
                per_aln_len,
                sensitive,
                self_search,
                output_dir):
        """Perform similarity search between query and target files.

        Parameters
        ----------
        similarity_method : str
            Similarity search method to use.
        query_files : list
            FASTA files with query sequences.
        target_files : list
            FASTA files with target sequences.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        self_search : boolean
            Allow searching between genomes with same name.
        output_dir : str
            Directory to store blast results.
        """
        
        assert(similarity_method in ['AAI_DIAMOND', 'AAI_BLASTP-FAST', 'ANI_BLASTN', 'ANI_DC-MEGABLAST', 'ANI_MEGABLAST'])
    
        parallel = Parallel(self.cpus)

        tmp_db_dir = os.path.join(output_dir, 'tmp_db')
        make_sure_path_exists(tmp_db_dir)

        # build databases
        self.logger.info('Creating sequence database for each genome.')
        cmds = []
        for target_file in target_files:
            stem_name = os.path.splitext(os.path.split(target_file)[-1])[0]
            db_prefix = os.path.join(tmp_db_dir, stem_name)
            cmd = self._create_db_cmd(similarity_method, target_file, db_prefix)
            cmds.append(cmd)

        rtn = parallel.run(self._producer, 
                            self._consumer, 
                            cmds, 
                            self._progress_genomes if not self.logger.is_silent else None)

        # perform pairwise similarity search
        self.logger.info('Performing similarity search between genomes.')
        cmds = []
        for query_file in query_files:
            query_name = os.path.splitext(os.path.split(query_file)[-1])[0]
            for target_file in target_files:
                target_name = os.path.splitext(os.path.split(target_file)[-1])[0]
                
                if (query_name == target_name) and not self_search:
                    continue
        
                hit_table = os.path.join(output_dir, '%s_vs_%s.tsv' % (query_name, target_name))
                cmd = self._search_cmd(similarity_method,
                                        query_file,
                                        os.path.join(tmp_db_dir, target_name),
                                        hit_table,
                                        evalue,
                                        per_identity, 
                                        per_aln_len, 
                                        sensitive)

                cmds.append(cmd)

        rtn = parallel.run(self._producer, 
                            self._consumer, 
                            cmds, 
                            self._progress_comparisons if not self.logger.is_silent else None)
                            
        # clear up temporary files
        if tmp_db_dir:
            shutil.rmtree(tmp_db_dir)
