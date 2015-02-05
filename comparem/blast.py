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
import multiprocessing as mp


class Blast(object):
    def __init__(self):
        self.logger = logging.getLogger()

        self._checkForBlast()

    def _checkForBlast(self):
        """Check to see if Blast is on the system before we try to run it."""
        try:
            subprocess.call(['blastp', '-help'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  Make sure blastp is on your system path.")
            sys.exit()

    def __workerThread(self, outputDir, evalue, queueIn, queueOut):
        while True:
            aaGeneFileA, aaGeneFileB = queueIn.get(block=True, timeout=None)
            if aaGeneFileA == None:
                break

            genomeIdA = ntpath.basename(aaGeneFileA)
            genomeIdA = genomeIdA[0:genomeIdA.rfind('.genes.faa')]

            genomeIdB = ntpath.basename(aaGeneFileB)
            genomeIdB = genomeIdB[0:genomeIdB.rfind('.genes.faa')]

            dbA = os.path.join(outputDir, genomeIdA + '.db')
            dbB = os.path.join(outputDir, genomeIdB + '.db')

            outputFileAB = os.path.join(outputDir, genomeIdA + '-' + genomeIdB + '.blastp.tsv')
            cmd = "blastp -query %s -db %s -out %s -max_target_seqs 1 -evalue %s -outfmt '6 qseqid qlen sseqid slen length pident evalue bitscore'" % (aaGeneFileA, dbB, outputFileAB, str(evalue))
            os.system(cmd)

            outputFileBA = os.path.join(outputDir, genomeIdB + '-' + genomeIdA + '.blastp.tsv')
            cmd = "blastp -query %s -db %s -out %s -max_target_seqs 1 -evalue %s -outfmt '6 qseqid qlen sseqid slen length pident evalue bitscore'" % (aaGeneFileB, dbA, outputFileBA, str(evalue))
            os.system(cmd)

            queueOut.put((aaGeneFileA, aaGeneFileB))

    def __writerThread(self, numDataItems, writerQueue):
        processedItems = 0
        while True:
            aaGeneFileA, _ = writerQueue.get(block=True, timeout=None)
            if aaGeneFileA == None:
                break

            processedItems += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genome pairs.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, aaGeneFiles, evalue, outputDir, cpus):
        # create the blast databases in serial
        self.logger.info('  Creating blast databases.')

        for aaGeneFile in aaGeneFiles:
            genomeId = ntpath.basename(aaGeneFile)
            genomeId = genomeId[0:genomeId.rfind('.genes.faa')]

            blastDB = os.path.join(outputDir, genomeId + '.db')
            logFile = os.path.join(outputDir, genomeId + '.log')
            if not os.path.exists(blastDB):
                cmd = 'makeblastdb -dbtype prot -in %s -out %s -logfile %s' % (aaGeneFile, blastDB, logFile)
                os.system(cmd)

        # perform blast in parallel
        self.logger.info('')
        self.logger.info('  Identifying blast hits between all pairs of genomes:')

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        numPairs = 0
        for i in xrange(0, len(aaGeneFiles)):
            for j in xrange(i + 1, len(aaGeneFiles)):
                workerQueue.put((aaGeneFiles[i], aaGeneFiles[j]))
                numPairs += 1

        for _ in range(cpus):
            workerQueue.put((None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(outputDir, evalue, workerQueue, writerQueue)) for _ in range(cpus)]
            writeProc = mp.Process(target=self.__writerThread, args=(numPairs, writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

                writeProc.terminate()
