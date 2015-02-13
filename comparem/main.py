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

from comparem.time_keeper import TimeKeeper
from comparem.prodigal import Prodigal
from comparem.blast import Blast
from comparem.diamond import Diamond
from comparem.aai_calculator import AAICalculator
from comparem.codon_usage import CodonUsage
from comparem.amino_acid_usage import AminoAcidUsage
from comparem.dinucleotide_usage import DinucleotideUsage
from comparem.PCoA import PCoA
from comparem.common import (remove_extension,
                             make_sure_path_exists,
                             check_dir_exists)


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - call_genes] Identifying genes within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)

        genome_files = []
        for f in os.listdir(options.genome_dir):
            if f.endswith(options.genome_ext):
                genome_files.append(os.path.join(options.genome_dir, f))

        if not genome_files:
            self.logger.warning('  [Warning] No genomes found. Check the --genome_ext flag used to identify genomes.')
            sys.exit()

        prodigal = Prodigal(options.cpus, options.genes, options.output_dir)
        prodigal.run(genome_files)

        self.logger.info('')
        self.logger.info('  Identified genes written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def blast(self, options):
        """Blast command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - blast] Performing reciprocal blast between genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)
        make_sure_path_exists(options.output_dir)

        aa_gene_files = []
        for f in os.listdir(options.gene_dir):
            if f.endswith(options.gene_ext):
                aa_gene_files.append(os.path.join(options.gene_dir, f))

        if not aa_gene_files:
            self.logger.warning('  [Warning] No gene files found. Check the --gene_ext flag used to identify gene files.')
            sys.exit()

        if options.diamond:
            diamond = Diamond(options.cpus, options.evalue, options.gene_ext, options.output_dir)
            diamond.run(aa_gene_files)
        else:
            blast = Blast(options.cpus, options.evalue, options.gene_ext, options.output_dir)
            blast.run(aa_gene_files)

        self.logger.info('')
        self.logger.info('  Reciprocal blast hits written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def aai(self, options):
        """AAI command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - aai] Identifying the AAI between homologs in genome pairs.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.gene_dir)
        check_dir_exists(options.blast_dir)
        make_sure_path_exists(options.output_dir)

        genome_ids = []
        for f in os.listdir(options.gene_dir):
            if f.endswith(options.gene_ext):
                genome_id = remove_extension(f, options.gene_ext)
                genome_ids.append(genome_id)

        if not genome_ids:
            self.logger.warning('  [Warning] No gene files found. Check the extension (-x) used to identify gene files.')
            sys.exit()

        aai_calculator = AAICalculator()
        aai_calculator.run(genome_ids,
                            options.gene_dir,
                            options.blast_dir,
                            options.gene_ext,
                            options.per_identity,
                            options.per_aln_len,
                            options.output_dir,
                            options.cpus)

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
        amino_acid_usage = AminoAcidUsage()
        genome_aa_usage, aa_set = amino_acid_usage.run(gene_files, options.cpus)

        # write out results
        fout = open(options.output_file, 'w')
        for aa in aa_set:
            fout.write('\t' + aa)
        fout.write('\n')

        for genome_id, codons in genome_aa_usage.iteritems():
            fout.write(genome_id)

            for aa in aa_set:
                fout.write('\t%d' % codons.get(aa, 0))
            fout.write('\n')

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
        genome_codon_usage, codon_set = codon_usage.run(gene_files)

        # write out results
        fout = open(options.output_file, 'w')
        for codon in codon_set:
            fout.write('\t' + codon)
        fout.write('\n')

        for genome_id, codons in genome_codon_usage.iteritems():
            fout.write(genome_id)

            for codon in codon_set:
                fout.write('\t%d' % codons.get(codon, 0))
            fout.write('\n')

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

    def lgt_usage(self, options):
        """LGT dinucleotide usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - lgt_usage] Calculating dinuceotide (3rd,1st) usage.')
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
        dinucleotide_usage = DinucleotideUsage(options.output_dir, options.cpus, options.keep_ambiguous)
        dinucleotide_usage.run(gene_files)

        self.logger.info('')
        self.logger.info('  Dinucleotide usage written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - unique] Identifying genes present in a single genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        self.time_keeper.print_time_stamp()

    def pcoa_plot(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CompareM - pcoa_plot] Generating PCoA plot showing relative similarity of genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        self.logger.info('  Performing PCoA.\n')
        pcoa = PCoA()
        pcoa.plot(options.aai_summary_file)

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
        elif(options.subparser_name == 'blast'):
            self.blast(options)
        elif(options.subparser_name == 'aai'):
            self.aai(options)
        elif(options.subparser_name == 'aai_wf'):
            root_dir = options.output_dir
            make_sure_path_exists(root_dir)

            options.output_dir = os.path.join(root_dir, 'genes')
            self.call_genes(options)

            options.output_dir = os.path.join(root_dir, 'blast')
            options.gene_dir = os.path.join(root_dir, 'genes')
            options.gene_ext = '.genes.faa'
            self.blast(options)

            options.output_dir = root_dir
            options.blast_dir = os.path.join(root_dir, 'blast')
            self.aai(options)
        elif(options.subparser_name == 'aa_usage'):
            self.aa_usage(options)
        elif(options.subparser_name == 'codon_usage'):
            self.codon_usage(options)
        elif(options.subparser_name == 'stop_usage'):
            self.stop_usage(options)
        elif(options.subparser_name == 'lgt_usage'):
            self.lgt_usage(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'pcoa_plot'):
            self.pcoa_plot(options)
        else:
            self.logger.error('  [Error] Unknown AAI command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
