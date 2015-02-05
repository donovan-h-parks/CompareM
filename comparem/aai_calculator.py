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

from comparem.seq_io import SeqIO


class AAICalculator(object):
    """Calculate AAI between all pairs of genomes."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def _blast_hits(self, blastTable, perIdentityThreshold, perAlnLenThreshold):
        """Identify homologous genes among BLAST hits.

        Determines the best hit for each query sequence which
        satisfies the conditions of being a homologous gene.

        Parameters
        ----------
        blastTable : str
            Name of table containing BLAST hits.
        perIdentityThreshold : float
            Percent identity threshold used to define a homologous gene.
        perAlnLenThreshold : float
            Alignment length threshold used to define a homologous gene.

        Returns
        -------
        dict
           Parameters of top hit to a homolgous gene for each query sequence.
        """

        hits = {}
        for line in open(blastTable):
            lineSplit = line.split('\t')

            querySeqId = lineSplit[0]
            queryLen = int(lineSplit[1])

            subSeqId = lineSplit[2]
            # subjectLen = int(lineSplit[3])

            alnLen = int(lineSplit[4])
            perIdent = float(lineSplit[5])
            evalue = float(lineSplit[6])
            bitscore = float(lineSplit[7])

            if perIdent >= perIdentityThreshold:
                perAlnLen = alnLen * 100.0 / queryLen

                if perAlnLen >= perAlnLenThreshold:
                    if querySeqId not in hits:  # take first hit passing criteria
                        hits[querySeqId] = [subSeqId, perIdent, perAlnLen, evalue, bitscore]

        return hits

    def __workerThread(self, outputDir, perIdentityThreshold, perAlnLenThreshold, queueIn, queueOut):

        genesOutputDir = os.path.join(outputDir, 'called_genes')
        blastOutputDir = os.path.join(outputDir, 'blastp')
        aaiOutputDir = os.path.join(outputDir, 'aai')

        while True:
            genomeIdA, genomeIdB = queueIn.get(block=True, timeout=None)
            if genomeIdA == None:
                break

            seqIO = SeqIO()

            # count number of genes in each genome
            genesInGenomeA = seqIO.readFasta(os.path.join(genesOutputDir, genomeIdA + '.genes.faa'))
            genesInGenomeB = seqIO.readFasta(os.path.join(genesOutputDir, genomeIdB + '.genes.faa'))

            # find blast hits between genome A and B
            outputFileAB = os.path.join(blastOutputDir, genomeIdA + '-' + genomeIdB + '.blastp.tsv')
            blastHitsAB = self._blastHits(outputFileAB, perIdentityThreshold, perAlnLenThreshold)

            # find blast hits between genomes B and A
            outputFileBA = os.path.join(blastOutputDir, genomeIdB + '-' + genomeIdA + '.blastp.tsv')
            blastHitsBA = self._blastHits(outputFileBA, perIdentityThreshold, perAlnLenThreshold)

            # find reciprocal best blast hits
            foutSeqs = open(os.path.join(aaiOutputDir, genomeIdA + '-' + genomeIdB + '.shared_genes.faa'), 'w')

            foutStats = open(os.path.join(aaiOutputDir, genomeIdA + '-' + genomeIdB + '.rbb_hits.tsv'), 'w')
            foutStats.write(genomeIdA + '\t' + genomeIdB + '\tPercent Identity\tPercent Alignment Length\te-value\tbitscore\n')

            perIdentityHits = []
            for querySeqId, stats in blastHitsAB.iteritems():
                subSeqId, perIdent, perAlnLen, evalue, bitscore = stats
                if subSeqId in blastHitsBA and querySeqId == blastHitsBA[subSeqId][0]:
                    foutStats.write('%s\t%s\t%.2f\t%.2f\t%.2g\t%.2f\n' % (querySeqId, subSeqId, perIdent, perAlnLen, evalue, bitscore))

                    # take average of percent identity in both blast directions as
                    # the results will be similar, but not identical
                    avgPerIdentity = 0.5 * (perIdent + blastHitsBA[subSeqId][1])
                    perIdentityHits.append(avgPerIdentity)

                    # write out shared genes
                    foutSeqs.write('>' + querySeqId + ' ' + genomeIdA + '\n')
                    foutSeqs.write(genesInGenomeA[querySeqId] + '\n')

                    foutSeqs.write('>' + subSeqId + ' ' + genomeIdB + '\n')
                    foutSeqs.write(genesInGenomeB[subSeqId] + '\n')

            foutSeqs.close()
            foutStats.close()

            meanPerIdentityHits = 0
            if len(perIdentityHits) > 0:
                meanPerIdentityHits = mean(perIdentityHits)

            stdPerIdentityHits = 0
            if len(perIdentityHits) >= 2:
                stdPerIdentityHits = std(perIdentityHits)

            queueOut.put((genomeIdA, genomeIdB, len(genesInGenomeA), len(genesInGenomeB), len(perIdentityHits), meanPerIdentityHits, stdPerIdentityHits))

    def __writerThread(self, numDataItems, outputDir, writerQueue):
        outputFile = os.path.join(outputDir, 'aai_summary.tsv')
        fout = open(outputFile, 'w')
        fout.write('Genome Id A\tGenes in A\tGenome Id B\tGenes in B\tOrthologous Genes\tMean AAI\tStd AAI\n')

        processedItems = 0
        while True:
            genomeIdA, genomeIdB, genesInGenomeA, genesInGenomeB, rbbHits, meanPerIdentityHits, stdPerIdentityHits = writerQueue.get(block=True, timeout=None)
            if genomeIdA == None:
                break

            processedItems += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genome pairs.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            fout.write('%s\t%d\t%s\t%d\t%d\t%.2f\t%.2f\n' % (genomeIdA, genesInGenomeA, genomeIdB, genesInGenomeB, rbbHits, meanPerIdentityHits, stdPerIdentityHits))

        sys.stdout.write('\n')

        self.logger.info('')
        self.logger.info('  Summary of AAI between genomes: %s' % outputFile)

        fout.close()

    def run(self, genomeIds, perIdentityThreshold, perAlnLenThreshold, outputDir, cpus):
        self.logger.info('  Calculating amino acid identity between all pairs of genomes:')

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        genomeIds.sort(key=str.lower)

        numPairs = 0
        for i in xrange(0, len(genomeIds)):
            for j in xrange(i + 1, len(genomeIds)):
                workerQueue.put((genomeIds[i], genomeIds[j]))
                numPairs += 1

        for _ in range(cpus):
            workerQueue.put((None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(outputDir, perIdentityThreshold, perAlnLenThreshold, workerQueue, writerQueue)) for _ in range(cpus)]
            writeProc = mp.Process(target=self.__writerThread, args=(numPairs, outputDir, writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None, None, None, None, None, None))
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

                writeProc.terminate()
