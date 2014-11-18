#!/usr/bin/env python

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
import multiprocessing as mp

from numpy import mean, std

class CalculateAAI(object):
    def __init__(self):
        pass
        
    def __countGenes(self, geneFile):
        geneCount = 0
        for line in open(geneFile):
            if line[0] == '>':
                geneCount += 1
                
        return geneCount
        
    def __blastHits(self, blastTable, perIdentityThreshold, perAlnLenThreshold):
        hits = {}
        for line in open(blastTable):
            lineSplit = line.split('\t')
            
            querySeqId = lineSplit[0]
            queryLen = int(lineSplit[1])
            
            subSeqId = lineSplit[2]
            subjectLen = int(lineSplit[3])
            
            alnLen = int(lineSplit[4])
            perIdent = float(lineSplit[5])
            evalue = float(lineSplit[6])
            bitscore = float(lineSplit [7])
            
            if perIdent >= perIdentityThreshold:
                perAlnLen = alnLen * 100.0 / queryLen
                
                if perAlnLen >= perAlnLenThreshold:
                    if querySeqId not in hits: # take first hit passing criteria
                        hits[querySeqId] = [subSeqId, perIdent, perAlnLen, evalue, bitscore]
            
        return hits
          
    def __workerThread(self, outputDir, perIdentity, perAlnLen, queueIn, queueOut):
        
        genesOutputDir = os.path.join(outputDir, 'called_genes')
        blastOutputDir = os.path.join(outputDir, 'blastp')
        aaiOutputDir = os.path.join(outputDir, 'aai')
        
        while True:
            genomeIdA, genomeIdB = queueIn.get(block=True, timeout=None)
            if genomeIdA == None:
                break
            
            # count number of genes in each genome
            genesInGenomeA = self.__countGenes(os.path.join(genesOutputDir, genomeIdA + '.genes.faa'))
            genesInGenomeB = self.__countGenes(os.path.join(genesOutputDir, genomeIdB + '.genes.faa'))
            
            # find blast hits between genome A and B
            outputFileAB = os.path.join(blastOutputDir, genomeIdA + '-' + genomeIdB + '.blastp.tsv')
            blastHitsAB = self.__blastHits(outputFileAB, perIdentity, perAlnLen)
            
            # find blast hits between genomes B and A
            outputFileBA = os.path.join(blastOutputDir, genomeIdB + '-' + genomeIdA + '.blastp.tsv')
            blastHitsBA = self.__blastHits(outputFileBA, perIdentity, perAlnLen)

            # find reciprocal best blast hits
            fout = open(os.path.join(aaiOutputDir, genomeIdB + '-' + genomeIdA + '.rbb_hits.tsv'), 'w')
            fout.write(genomeIdA + '\t' + genomeIdB + '\tPercent Identity\tPercent Alignment Length\te-value\tbitscore\n')
            
            rbbHits = 0
            perIdentityHits = []
            for querySeqId, stats in blastHitsAB.iteritems():
                subSeqId, perIdent, perAlnLen, evalue, bitscore = stats
                if subSeqId in blastHitsBA and querySeqId == blastHitsBA[subSeqId][0]:
                    fout.write('%s\t%s\t%.2f\t%.2f\t%f\t%.2f\n' % (querySeqId, subSeqId, perIdent, perAlnLen, evalue, bitscore))
                    rbbHits += 1
                    perIdentityHits.append(perIdent)
            fout.close()
    
            meanPerIdentityHits = 0
            if len(perIdentityHits) > 0:
                meanPerIdentityHits = mean(perIdentityHits)
            
            stdPerIdentityHits = 0
            if len(perIdentityHits) >= 2:
                stdPerIdentityHits = std(perIdentityHits)
            queueOut.put((genomeIdA, genomeIdB, genesInGenomeA, genesInGenomeB, rbbHits, meanPerIdentityHits, stdPerIdentityHits))

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
            statusStr = '    Finished processing %d of %d (%.2f%%) genome pairs.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()
            
            fout.write('%s\t%d\t%s\t%d\t%d\t%.2f\t%.2f\n' % (genomeIdA, genesInGenomeA, genomeIdB, genesInGenomeB, rbbHits, meanPerIdentityHits, stdPerIdentityHits))

        sys.stdout.write('\n')
        
        print ''
        print '  Summary of AAI between genomes: %s' % outputFile
        
        fout.close()

    def run(self, genomeIds, perIdentity, perAlnLen, outputDir, threads):   
        print '  Calculating amino acid identity between all pairs of genomes:'
             
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        numPairs = 0
        for i in xrange(0, len(genomeIds)):            
            for j in xrange(i+1, len(genomeIds)):
                workerQueue.put((genomeIds[i], genomeIds[j]))
                numPairs += 1
        
        for _ in range(threads):
            workerQueue.put((None, None))
        
        try:
            workerProc = [mp.Process(target = self.__workerThread, args = (outputDir, perIdentity, perAlnLen, workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target = self.__writerThread, args = (numPairs, outputDir, writerQueue))
            
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

