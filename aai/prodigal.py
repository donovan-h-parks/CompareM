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
import subprocess
import ntpath
import multiprocessing as mp

class Prodigal(object):
    def __init__(self):
        self.__checkForProdigal()
    
    def __checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it."""
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "  [Error] Make sure prodigal is on your system path."
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
            
            if not bCalledGenes:
                cmd = ('prodigal -p single -q -m -g 11 -a %s -d %s -i %s > %s 2> /dev/null' % (aaGeneFile, ntGeneFile, genomeFile, gffFile))
                os.system(cmd)
            else:
                os.system('ln -s %s %s' % (os.path.abspath(genomeFile), aaGeneFile))
    
            queueOut.put(genomeId)

    def __writerThread(self, numDataItems, writerQueue):
        processedItems = 0
        while True:
            genomeId = writerQueue.get(block=True, timeout=None)
            if genomeId == None:
                break
            
            processedItems += 1
            statusStr = '    Finished processing %d of %d (%.2f%%) genomes.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genomeFiles, bCalledGenes, outputDir, threads):
        print '  Identifying genes within genomes:'
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        for gf in genomeFiles:
            workerQueue.put(gf)
        
        for _ in range(threads):
            workerQueue.put(None)
        
        try:
            workerProc = [mp.Process(target = self.__workerThread, args = (bCalledGenes, outputDir, workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target = self.__writerThread, args = (len(genomeFiles), writerQueue))
            
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

