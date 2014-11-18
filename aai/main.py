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

from aai.timeKeeper import TimeKeeper

from aai.prodigal import Prodigal
from aai.blast import Blast
from aai.calculateAAI import CalculateAAI

class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.timeKeeper = TimeKeeper()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - call_genes] Identify genes within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        if not os.path.exists(options.output_dir):
            os.mkdir(options.output_dir)
        
        outputDir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        
        genomeFiles = []
        for f in os.listdir(options.genome_dir):
            if f.endswith(options.extension):
                genomeFiles.append(os.path.join(options.genome_dir, f))
                
        if not genomeFiles:
            print '  [Error] No genomes found. Check the extension (-x) used to identify genomes.'
            sys.exit()
        
        prodigal = Prodigal()
        prodigal.run(genomeFiles, options.bCalledGenes, outputDir, options.threads)

        self.timeKeeper.printTimeStamp()

    def ortholog(self, options):
        """Ortholog command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - ortholog] Identify orthologs within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        geneDir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(geneDir):
            print "  [Error] You must specify the output directory used with the 'call_genes' command"
            sys.exit(-1)
        
        outputDir = os.path.join(options.output_dir, 'blastp')
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
            
        aaGeneFiles = []
        for f in os.listdir(geneDir):
            if f.endswith('faa'):
                aaGeneFiles.append(os.path.join(geneDir, f))
                
        blast = Blast()
        blast.run(aaGeneFiles, options.evalue, outputDir, options.threads)

        self.timeKeeper.printTimeStamp()

    def calculate(self, options, db=None):
        """Calculate command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [AAI - calculate] Calculate AAI between orthologs')
        self.logger.info('*******************************************************************************')
        self.logger.info('')
        
        geneDir = os.path.join(options.output_dir, 'called_genes')
        if not os.path.exists(geneDir):
            print "  [Error] You must specify the output directory used with the 'call_genes' command"
            sys.exit(-1)
            
        outputDir = os.path.join(options.output_dir, 'aai')
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
            
        genomeIds = []
        for f in os.listdir(geneDir):
            if f.endswith('faa'):
                genomeId = f[0:f.rfind('.genes.faa')]
                genomeIds.append(genomeId)
        
        calculateAAI = CalculateAAI()
        calculateAAI.run(genomeIds, options.per_identity, options.per_aln_len, options.output_dir, options.threads)

        self.timeKeeper.printTimeStamp()

    def parseOptions(self, options):
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
        else:
            self.logger.error('  [Error] Unknown AAI command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
