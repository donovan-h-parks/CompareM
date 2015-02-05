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

from aai.time_keeper import TimeKeeper

from aai.prodigal import Prodigal
from aai.blast import Blast
from aai.aai_calculator import AAICalculator
from aai.codon_usage import CodonUsage
from aai.amino_acid_usage import AminoAcidUsage
from aai.PCoA import PCoA


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.timeKeeper = TimeKeeper()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - call_genes] Identifying genes within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        if not os.path.exists(options.output_dir):
            os.mkdir(options.output_dir)

        output_dir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        genome_files = []
        for f in os.listdir(options.genome_dir):
            if f.endswith(options.extension):
                genome_files.append(os.path.join(options.genome_dir, f))

        if not genome_files:
            print '  [Error] No genomes found. Check the extension (-x) used to identify genomes.'
            sys.exit()

        prodigal = Prodigal()
        prodigal.run(genome_files, options.bCalledGenes, output_dir, options.cpus)

        self.timeKeeper.printTimeStamp()

    def ortholog(self, options):
        """Ortholog command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - ortholog] Identifying orthologs within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        gene_dir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(gene_dir):
            print "  [Error] You must specify the output directory used with the 'call_genes' command"
            sys.exit(-1)

        output_dir = os.path.join(options.output_dir, 'blastp')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        aaGeneFiles = []
        for f in os.listdir(gene_dir):
            if f.endswith('faa'):
                aaGeneFiles.append(os.path.join(gene_dir, f))

        blast = Blast()
        blast.run(aaGeneFiles, options.evalue, output_dir, options.cpus)

        self.timeKeeper.printTimeStamp()

    def calculate(self, options, db=None):
        """Calculate command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - calculate] Calculating AAI between orthologs')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        gene_dir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(gene_dir):
            print "  [Error] You must specify the output directory used with the 'call_genes' command"
            sys.exit(-1)

        output_dir = os.path.join(options.output_dir, 'aai')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        genomeIds = []
        for f in os.listdir(gene_dir):
            if f.endswith('faa'):
                genomeId = f[0:f.rfind('.genes.faa')]
                genomeIds.append(genomeId)

        aai_calculator = AAICalculator()
        aai_calculator.run(genomeIds, options.per_identity, options.per_aln_len, options.output_dir, options.cpus)

        self.timeKeeper.printTimeStamp()

    def aa_usage(self, options, db=None):
        """Amino acid usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - aa_usage] Calculating amino acid usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        if not os.path.exists(options.input_dir):
            self.logger.error('  Input directory does not exists: %s' % options.input_dir)
            sys.exit()

        # get list of files with called genes
        geneFiles = []
        files = os.listdir(options.input_dir)
        for f in files:
            if f.endswith(options.extension):
                geneFiles.append(os.path.join(options.input_dir, f))

        # warn use if no files were found
        if len(geneFiles) == 0:
            self.logger.warning('  No files found. Perhaps you need to set the extension (-x) flag.')
            return

        # calculate amino acid usage
        aminoAcidUsage = AminoAcidUsage()
        aaInGenomes, aaSet = aminoAcidUsage.run(geneFiles, options.cpus)

        # write out results
        fout = open(options.output_file, 'w')
        for aa in aaSet:
            fout.write('\t' + aa)
        fout.write('\n')

        for genomeId, codons in aaInGenomes.iteritems():
            fout.write(genomeId)

            for aa in aaSet:
                fout.write('\t%d' % codons.get(aa, 0))
            fout.write('\n')

        self.logger.info('')
        self.logger.info('  Amino acid usage written to: %s' % options.output_file)

        self.timeKeeper.printTimeStamp()

    def codon_usage(self, options, db=None):
        """Codon usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - codon_usage] Calculating codon usage within each genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        if not os.path.exists(options.input_dir):
            self.logger.error('  Input directory does not exists: %s' % options.input_dir)
            sys.exit()

        # get list of files with called genes
        geneFiles = []
        files = os.listdir(options.input_dir)
        for f in files:
            if f.endswith(options.extension):
                geneFiles.append(os.path.join(options.input_dir, f))

        # warn use if no files were found
        if len(geneFiles) == 0:
            self.logger.warning('  No files found. Perhaps you need to set the extension (-x) flag.')
            return

        # calculate amino acid usage
        codonUsage = CodonUsage()
        codonsInGenomes, codonSet = codonUsage.run(geneFiles, options.keep_ambiguous)

        # write out results
        fout = open(options.output_file, 'w')
        for aa in codonSet:
            fout.write('\t' + aa)
        fout.write('\n')

        for genomeId, codons in codonsInGenomes.iteritems():
            fout.write(genomeId)

            for aa in codonSet:
                fout.write('\t%d' % codons.get(aa, 0))
            fout.write('\n')

        self.logger.info('')
        self.logger.info('  Codon usage written to: %s' % options.output_file)

        self.timeKeeper.printTimeStamp()

    def unique(self, options, db=None):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - unique] Identifying genes present in a single genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        self.timeKeeper.printTimeStamp()

    def pcoa_plot(self, options, db=None):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - pcoa_plot] Generating PCoA plot showing relative similarity of genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        self.logger.info('  Performing PCoA.\n')
        pcoa = PCoA()
        pcoa.plot(options.aai_summary_file)

        self.timeKeeper.printTimeStamp()

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
        elif(options.subparser_name == 'ortholog'):
            self.ortholog(options)
        elif(options.subparser_name == 'calculate'):
            self.calculate(options)
        elif(options.subparser_name == 'aai_wf'):
            self.call_genes(options)
            self.ortholog(options)
            self.calculate(options)
        elif(options.subparser_name == 'aa_usage'):
            self.aa_usage(options)
        elif(options.subparser_name == 'codon_usage'):
            self.codon_usage(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'pcoa_plot'):
            self.pcoa_plot(options)
        else:
            self.logger.error('  [Error] Unknown AAI command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
