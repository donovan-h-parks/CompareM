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

import os
import sys
import logging
from collections import defaultdict

from comparem.reciprocal_blast import ReciprocalBlast
from comparem.reciprocal_diamond import ReciprocalDiamond
from comparem.aai_calculator import AAICalculator
from comparem.codon_usage import CodonUsage
from comparem.amino_acid_usage import AminoAcidUsage
from comparem.kmer_usage import KmerUsage
from comparem.lgt_dinucleotide import LgtDinucleotide
from comparem.lgt_codon import LgtCodon
from comparem.PCoA import PCoA
from comparem.plots.heatmap import Heatmap

import biolib.seq_io as seq_io
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.prodigal import Prodigal
from biolib.common import (remove_extension,
                             make_sure_path_exists,
                             check_dir_exists,
                             concatenate_files)


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def _genome_files(self, genome_dir, genome_ext):
        """Identify genomes files.

        Parameters
        ----------
        genome_dir : str
            Directory containing genomes of interest.
        genome_ext : str
            Extension of genome files.

        Returns
        -------
        list
            Name of genome files in directory.
        """

        check_dir_exists(genome_dir)

        genome_files = []
        for f in os.listdir(genome_dir):
            if f.endswith(genome_ext):
                genome_files.append(os.path.join(genome_dir, f))

        if not genome_files:
            self.logger.warning('  [Warning] No genomes found. Check the --genome_ext flag used to identify genomes.')
            sys.exit()

        return genome_files

    def _write_usage_profile(self, genome_usage, feature_set, output_file):
        """Write out occurrence of specified features for each genome.

        Parameters
        ----------
        genome_usage : d[genome_id][feature] -> count
            Occurrence of genomic feature in genome
        feature_set : iterable
            All genomic features.
        output_file : str
            File to produce.
        """

        sorted_feature_set = sorted(feature_set)

        fout = open(output_file, 'w')
        fout.write('Genome ID')
        for feature in sorted_feature_set:
            fout.write('\t' + feature)
        fout.write('\n')

        totals = defaultdict(int)
        for genome_id, features in genome_usage.iteritems():
            for feature in sorted_feature_set:
                totals[genome_id] += features.get(feature, 0)

        for genome_id, features in genome_usage.iteritems():
            fout.write(genome_id)

            for feature in sorted_feature_set:
                fout.write('\t%.2f%%' % (features.get(feature, 0) * 100.0 / totals[genome_id]))
            fout.write('\n')

    def ani(self, options):
        """ANI command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - ani] Calculating the ANI between genome pairs.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)

        self.logger.info('')
        self.logger.info('  Average nucleotide identity information written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - call_genes] Identifying genes within genomes.')
        self.logger.info('*******************************************************************************')

        make_sure_path_exists(options.output_dir)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)
        if not genome_files:
            self.logger.warning('  [Warning] No genome files found. Check the --genome_ext flag used to identify genomes.')
            sys.exit()

        prodigal = Prodigal(options.cpus)
        summary_stats = prodigal.run(genome_files, False, options.force_table, False, options.output_dir)

        # write gene calling summary
        fout = open(os.path.join(options.output_dir, 'call_genes.summary.tsv'), 'w')
        fout.write('Genome Id\tSelected translation table\tTable 4 coding density\tTable 11 coding density\n')
        for genome_id, stats in summary_stats.iteritems():
            fout.write('%s\t%d\t%.2f%%\t%.2f%%\n' % (genome_id,
                                                     stats.best_translation_table,
                                                     stats.coding_density_4,
                                                     stats.coding_density_11))
        fout.close()

        self.logger.info('')
        self.logger.info('  Identified genes written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def rblast(self, options):
        """Reciprocal blast command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - rblast] Performing reciprocal blast between genomes.')
        self.logger.info('*******************************************************************************')

        check_dir_exists(options.protein_dir)
        make_sure_path_exists(options.output_dir)

        aa_gene_files = []
        for f in os.listdir(options.protein_dir):
            if f.endswith(options.protein_ext):
                aa_gene_files.append(os.path.join(options.protein_dir, f))

        if not aa_gene_files:
            self.logger.warning('  [Warning] No gene files found. Check the --protein_ext flag used to identify gene files.')
            sys.exit()

        # modify gene ids to include genome ids in order to ensure
        # all gene identifiers are unique across the set of genomes,
        # also removes the trailing asterisk used to identify the stop
        # codon
        self.logger.info('')
        self.logger.info('  Appending genome identifiers to all gene identifiers.')
        gene_out_dir = os.path.join(options.output_dir, 'genes')
        make_sure_path_exists(gene_out_dir)
        modified_aa_gene_files = []
        for gf in aa_gene_files:
            genome_id = remove_extension(gf)

            aa_file = os.path.join(gene_out_dir, genome_id + '.faa')
            fout = open(aa_file, 'w')
            for seq_id, seq, annotation in seq_io.read_fasta_seq(gf, keep_annotation=True):
                fout.write('>' + seq_id + '~' + genome_id + ' ' + annotation + '\n')
                if seq[-1] == '*':
                    seq = seq[0:-1]
                fout.write(seq + '\n')
            fout.close()

            modified_aa_gene_files.append(aa_file)

        # perform the reciprocal blast with blastp or diamond
        self.logger.info('')
        if options.blastp:
            rblast = ReciprocalBlast(options.cpus)
            rblast.run(modified_aa_gene_files, options.evalue, options.output_dir)

            # concatenate all blast tables to mimic output of diamond, all hits
            # for a given genome MUST be in consecutive order to fully mimic
            # the expected results from diamond
            self.logger.info('')
            self.logger.info('  Creating single file with all blast hits (be patient!).')
            blast_files = sorted([f for f in os.listdir(options.output_dir) if f.endswith('.blastp.tsv')])
            hit_tables = [os.path.join(options.output_dir, f) for f in blast_files]
            concatenate_files(hit_tables, os.path.join(options.output_dir, 'all_hits.tsv'))
        else:
            rdiamond = ReciprocalDiamond(options.cpus)
            rdiamond.run(modified_aa_gene_files, options.evalue, options.per_identity, options.output_dir)

        self.logger.info('')
        self.logger.info('  Reciprocal blast hits written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def aai(self, options):
        """AAI command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - aai] Calculating the AAI between homologs in genome pairs.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.rblast_dir)
        make_sure_path_exists(options.output_dir)

        genome_ids = []
        protein_dir = os.path.join(options.rblast_dir, 'genes')
        for f in os.listdir(protein_dir):
            if f.endswith('.faa'):
                genome_id = remove_extension(f, '.faa')
                genome_ids.append(genome_id)

        if not genome_ids:
            self.logger.warning('  [Warning] No gene files found. Check the --protein_ext flag used to identify gene files.')
            sys.exit()

        aai_calculator = AAICalculator(options.cpus)
        aai_calculator.run(genome_ids,
                            protein_dir,
                            options.rblast_dir,
                            options.per_identity,
                            options.per_aln_len,
                            options.write_shared_genes,
                            options.output_dir)

        shared_genes_dir = os.path.join(options.output_dir, aai_calculator.shared_genes)
        self.logger.info('')
        self.logger.info('  Identified homologs between genome pairs written to: %s' % shared_genes_dir)

        self.time_keeper.print_time_stamp()

    def aa_usage(self, options):
        """Amino acid usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - aa_usage] Calculating amino acid usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.protein_dir)

        # get list of files with called genes
        gene_files = []
        files = os.listdir(options.protein_dir)
        for f in files:
            if f.endswith(options.protein_ext):
                gene_files.append(os.path.join(options.protein_dir, f))

        # warn use if no files were found
        if len(gene_files) == 0:
            self.logger.warning('  [Warning] No gene files found. Check the --protein_ext flag used to identify gene files.')
            return

        # calculate amino acid usage
        amino_acid_usage = AminoAcidUsage(options.cpus)
        genome_aa_usage, aa_set = amino_acid_usage.run(gene_files)

        # write out results
        self._write_usage_profile(genome_aa_usage, aa_set, options.output_file)

        self.logger.info('')
        self.logger.info('  Amino acid usage written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def codon_usage(self, options):
        """Codon usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - codon_usage] Calculating codon usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)

        # get list of files with called genes
        gene_files = []
        files = os.listdir(options.gene_dir)
        for f in files:
            if f.endswith(options.gene_ext):
                gene_files.append(os.path.join(options.gene_dir, f))

        # warn use if no files were found
        if len(gene_files) == 0:
            self.logger.warning('  [Warning] No gene files found. Check the --gene_ext flag used to identify gene files.')
            return

        # calculate amino acid usage
        codon_usage = CodonUsage(options.cpus, options.keep_ambiguous)
        genome_codon_usage, codon_set, _mean_length = codon_usage.run(gene_files)

        # write out results
        self._write_usage_profile(genome_codon_usage, codon_set, options.output_file)

        self.logger.info('')
        self.logger.info('  Codon usage written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def stop_usage(self, options):
        """Stop codon usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - stop_usage] Calculating stop codon usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)

        # get list of files with called genes
        gene_files = []
        files = os.listdir(options.gene_dir)
        for f in files:
            if f.endswith(options.gene_ext):
                gene_files.append(os.path.join(options.gene_dir, f))

        # warn use if no files were found
        if len(gene_files) == 0:
            self.logger.warning('  [Warning] No gene files found. Check the --gene_ext flag used to identify gene files.')
            return

        # calculate amino acid usage
        codon_usage = CodonUsage(options.cpus, keep_ambiguous=False, stop_codon_only=True)
        genome_codon_usage, codon_set, mean_gene_length = codon_usage.run(gene_files)

        # write out results
        fout = open(options.output_file, 'w')
        for codon in codon_set:
            fout.write('\t' + codon)
            if mean_gene_length:
                fout.write('\t' + codon + ': avg. seq. length')
        fout.write('\n')

        for genome_id, codons in genome_codon_usage.iteritems():
            fout.write(genome_id)

            for codon in codon_set:
                fout.write('\t%d' % codons.get(codon, 0))

                if mean_gene_length:
                    mean_len = mean_gene_length[genome_id].get(codon, None)
                    if mean_len:
                        fout.write('\t%.1f' % mean_len)
                    else:
                        fout.write('\tna')
            fout.write('\n')

        self.logger.info('')
        self.logger.info('  Stop codon usage written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def kmer_usage(self, options):
        """Kmer usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - kmer_usage] Calculating kmer usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        if options.k > 10 or options.k <= 0:
            self.logger.warning('[Warning] CompareM only support kmers with k <= 10.')
            sys.exit(0)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)

        # calculate amino acid usage
        kmer_usage = KmerUsage(options.k, options.cpus)
        genome_kmer_usage, kmer_set = kmer_usage.run(genome_files)

        # write out results
        self.logger.info('')
        self.logger.info('  Writing kmer profile to file (be patient!).')
        self._write_usage_profile(genome_kmer_usage, kmer_set, options.output_file)

        self.logger.info('')
        self.logger.info('  Kmer usage written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def lgt_di(self, options):
        """LGT dinucleotide usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - lgt_di] Calculating dinuceotide (3rd,1st) usage of genes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)

        # get list of files with called genes
        gene_files = []
        files = os.listdir(options.gene_dir)
        for f in files:
            if f.endswith(options.gene_ext):
                gene_files.append(os.path.join(options.gene_dir, f))

        # warn use if no files were found
        if len(gene_files) == 0:
            self.logger.warning('  [Warning] No gene files found. Check the --gene_ext flag used to identify gene files.')
            return

        lgt_dinucleotide = LgtDinucleotide(options.cpus)
        lgt_dinucleotide.run(gene_files, options.crit_value, options.output_dir)

        self.logger.info('')
        self.logger.info('  Dinucleotide usage written to directory: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def lgt_codon(self, options):
        """LGT dinucleotide usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - lgt_codon] Calculating codon usage of genes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)

        # get list of files with called genes
        gene_files = []
        files = os.listdir(options.gene_dir)
        for f in files:
            if f.endswith(options.gene_ext):
                gene_files.append(os.path.join(options.gene_dir, f))

        # warn use if no files were found
        if len(gene_files) == 0:
            self.logger.warning('  [Warning] No gene files found. Check the --gene_ext flag used to identify gene files.')
            return

        lgt_codon = LgtCodon(options.cpus)
        lgt_codon.run(gene_files, options.output_dir)

        self.logger.info('')
        self.logger.info('  Codon usage written to directory: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - unique] Identifying genes present in a single genome.')
        self.logger.info('*******************************************************************************')

        self.time_keeper.print_time_stamp()

    def pcoa_plot(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - pcoa_plot] Generating PCoA plot showing relative similarity of genomes.')
        self.logger.info('*******************************************************************************')

        self.logger.info('')
        self.logger.info('  Performing PCoA.')
        pcoa = PCoA()
        pcoa.plot(options.aai_summary_file)

        self.time_keeper.print_time_stamp()

    def heatmap(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - heatmap] Generating heatmap showing relative similarity of genomes.')
        self.logger.info('*******************************************************************************')

        self.logger.info('')
        self.logger.info('  Making heatmap.')
        heatmapper = Heatmap(options.aai_summary_file, options.output_file)
        heatmapper.plot(options.cluster, options.method, options.metric)

        self.time_keeper.print_time_stamp()

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        try:
            if options.bVerbose:
                logging.basicConfig(format='', level=logging.DEBUG)
            elif options.bQuiet:
                logging.basicConfig(format='', level=logging.ERROR)
            else:
                logging.basicConfig(format='', level=logging.INFO)
        except:
            logging.basicConfig(format='', level=logging.INFO)

        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'call_genes'):
            self.call_genes(options)
        elif(options.subparser_name == 'rblast'):
            self.rblast(options)
        elif(options.subparser_name == 'aai'):
            self.aai(options)
        elif(options.subparser_name == 'aai_wf'):
            root_dir = options.output_dir
            make_sure_path_exists(root_dir)

            options.output_dir = os.path.join(root_dir, 'genes')
            self.call_genes(options)

            options.protein_ext = 'faa'
            options.protein_dir = os.path.join(root_dir, 'genes')
            options.output_dir = os.path.join(root_dir, 'rblast')
            self.rblast(options)

            options.output_dir = root_dir
            options.rblast_dir = os.path.join(root_dir, 'rblast')
            self.aai(options)
        elif(options.subparser_name == 'aa_usage'):
            self.aa_usage(options)
        elif(options.subparser_name == 'codon_usage'):
            self.codon_usage(options)
        elif(options.subparser_name == 'kmer_usage'):
            self.kmer_usage(options)
        elif(options.subparser_name == 'stop_usage'):
            self.stop_usage(options)
        elif(options.subparser_name == 'lgt_di'):
            self.lgt_di(options)
        elif(options.subparser_name == 'lgt_codon'):
            self.lgt_codon(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'pcoa_plot'):
            self.pcoa_plot(options)
        elif(options.subparser_name == 'heatmap'):
            self.heatmap(options)
        else:
            self.logger.error('  [Error] Unknown CompareM command: "' + options.subparser_name + '"\n')
            sys.exit()

        return 0
