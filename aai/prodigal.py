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

from aai.seq_io import SeqIO
from aai.common import checkFileExists

import numpy as np


class Prodigal(object):
    """Wrapper for running Prodigal."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self._checkForProdigal()

    def _checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it."""
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.info("  Make sure prodigal is on your system path.")
            sys.exit()

    def __workerThread(self, bCalledGenes, outputDir, queueIn, queueOut):
        while True:
            genomeFile = queueIn.get(block=True, timeout=None)
            if genomeFile == None:
                break

            genomeId = ntpath.basename(genomeFile)
            genomeId = genomeId[0:genomeId.rfind('.')]

            aaGeneFile = os.path.join(outputDir, genomeId + '.genes.faa')
            ntGeneFile = os.path.join(outputDir, genomeId + '.genes.fna')
            gffFile = os.path.join(outputDir, genomeId + '.gff')

            if bCalledGenes:
                os.system('ln -s %s %s' % (os.path.abspath(genomeFile), aaGeneFile))
            else:
                tmpDir = tempfile.mkdtemp()

                seqIO = SeqIO()
                seqs = seqIO.readFasta(genomeFile)

                # determine number of bases
                totalBases = 0
                for seq in seqs.values():
                    totalBases += len(seq)

                # call genes under different translation tables
                tableCodingDensity = {}
                for translationTable in [4, 11]:
                    aaGeneFileTmp = os.path.join(tmpDir, 'table' + str(translationTable), genomeId + '.genes.faa')
                    ntGeneFileTmp = os.path.join(tmpDir, 'table' + str(translationTable), genomeId + '.genes.fna')
                    gffFileTmp = os.path.join(tmpDir, 'table' + str(translationTable), genomeId + '.gff')

                    # check if there is sufficient bases to calculate prodigal parameters
                    if totalBases < 100000:
                        procedureStr = 'meta'  # use best precalculated parameters
                    else:
                        procedureStr = 'single'  # estimate parameters from data

                    cmd = 'prodigal -p %s -q -f gff -g %d -a %s -d %s -i %s > %s 2> /dev/null' % (procedureStr, translationTable, aaGeneFileTmp, ntGeneFileTmp, genomeFile, gffFileTmp)
                    os.system(cmd)

                    # determine coding density
                    prodigalParser = ProdigalGeneFeatureParser(gffFile)

                    codingBases = 0
                    for seqId, seq in seqs.iteritems():
                        codingBases += prodigalParser.codingBases(seqId)

                    codingDensity = float(codingBases) / totalBases
                    tableCodingDensity[translationTable] = codingDensity

                # determine best translation table
                bestTranslationTable = 11
                if (tableCodingDensity[4] - tableCodingDensity[11] > 0.05) and tableCodingDensity[4] > 0.7:
                    bestTranslationTable = 4

                shutil.copyfile(os.path.join(tmpDir, 'table' + str(bestTranslationTable), genomeId + '.genes.faa'), aaGeneFile)
                shutil.copyfile(os.path.join(tmpDir, 'table' + str(translationTable), genomeId + '.genes.fna'), ntGeneFile)
                shutil.copyfile(os.path.join(tmpDir, 'table' + str(translationTable), genomeId + '.gff'), gffFile)

                # clean up temporary files
                shutil.rmtree(tmpDir)

            queueOut.put(genomeId)

    def __writerThread(self, numDataItems, writerQueue):
        processedItems = 0
        while True:
            genomeId = writerQueue.get(block=True, timeout=None)
            if genomeId == None:
                break

            processedItems += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genomes.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genomeFiles, bCalledGenes, outputDir, cpus):
        self.logger.info('  Identifying genes within genomes:')

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for gf in genomeFiles:
            workerQueue.put(gf)

        for _ in range(cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(bCalledGenes, outputDir, workerQueue, writerQueue)) for _ in range(cpus)]
            writeProc = mp.Process(target=self.__writerThread, args=(len(genomeFiles), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

                writeProc.terminate()


class ProdigalGeneFeatureParser():
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename):
        """Initialization."""
        checkFileExists(filename)

        self.genes = {}
        self.lastCodingBase = {}

        self.__parseGFF(filename)

        self.codingBaseMasks = {}
        for seqId in self.genes:
            self.codingBaseMasks[seqId] = self.__buildCodingBaseMask(seqId)

    def __parseGFF(self, filename):
        """Parse genes from GFF file."""
        bGetTranslationTable = True
        for line in open(filename):
            if bGetTranslationTable and line.startswith('# Model Data'):
                self.translationTable = line.split(';')[4]
                self.translationTable = int(self.translationTable[self.translationTable.find('=') + 1:])
                bGetTranslationTable = False

            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId not in self.genes:
                geneCounter = 0
                self.genes[seqId] = {}
                self.lastCodingBase[seqId] = 0

            geneId = seqId + '_' + str(geneCounter)
            geneCounter += 1

            start = int(lineSplit[3])
            end = int(lineSplit[4])

            self.genes[seqId][geneId] = [start, end]
            self.lastCodingBase[seqId] = max(self.lastCodingBase[seqId], end)

    def __buildCodingBaseMask(self, seqId):
        """Build mask indicating which bases in a sequences are coding."""

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        codingBaseMask = np.zeros(self.lastCodingBase[seqId])
        for pos in self.genes[seqId].values():
            codingBaseMask[pos[0]:pos[1] + 1] = 1

        return codingBaseMask

    def codingBases(self, seqId, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end)."""

        # check if sequence has any genes
        if seqId not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end == None:
            end = self.lastCodingBase[seqId]

        return np.sum(self.codingBaseMasks[seqId][start:end])
