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

__prog_name__ = 'createMatrix'
__prog_desc__ = 'create matrix of AAI results'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
from collections import defaultdict

class MyClass(object):
    def __init__(self):
      pass

    def run(self, aaiSummaryFile, aaiMatrix, orthologyMatrix):
      # read pairwise comparisons
      matrix = defaultdict(dict)
      with open(aaiSummaryFile) as f:
        f.readline()
        for line in f:
          lineSplit = line.split('\t')
          genomeA = lineSplit[0]
          genomeB = lineSplit[2]
          orthologusGenes = int(lineSplit[4])
          aai = float(lineSplit[5])

          matrix[genomeA][genomeB] = [orthologusGenes, aai]
          matrix[genomeB][genomeA] = [orthologusGenes, aai]

      sampleIds = sorted(matrix.keys())

      # write out AAI matrix
      fout = open(aaiMatrix, 'w')

      fout.write('\t' + '\t'.join(sampleIds) + '\n')
      for i in xrange(0, len(sampleIds)):
        sampleIdI = sampleIds[i]
        row = matrix[sampleIdI]

        fout.write(sampleIdI)
        for j in xrange(0, len(sampleIds)):
          sampleIdJ = sampleIds[j]
          if i == j:
            fout.write('\t0')
          else:
            fout.write('\t%.2f' % row[sampleIdJ][1])
        fout.write('\n')
      fout.close()

      # write out number of orthologus protein matrix
      fout = open(orthologyMatrix, 'w')

      fout.write('\t' + '\t'.join(sampleIds) + '\n')
      for i in xrange(0, len(sampleIds)):
        sampleIdI = sampleIds[i]
        row = matrix[sampleIdI]

        fout.write(sampleIdI)
        for j in xrange(0, len(sampleIds)):
          sampleIdJ = sampleIds[j]
          if i == j:
            fout.write('\t1')
          else:
            fout.write('\t%d' % row[sampleIdJ][0])
        fout.write('\n')
      fout.close()

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('aai_summary_file', help='file indicating AAI between all pairs of genomes')
    parser.add_argument('aai_matrix', help='output AAI matrix')
    parser.add_argument('orthology_matrix', help='orthologous gene count matrix')

    args = parser.parse_args()

    try:
        myClass = MyClass()
        myClass.run(args.aai_summary_file, args.aai_matrix, args.orthology_matrix)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
