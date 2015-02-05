###############################################################################
#
# PCA.py - Principal component analysis
#
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

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
# from sklearn import manifold


class PCoA:
    """http://stackoverflow.com/questions/1730600/principal-component-analysis-in-python"""
    def __init__(self):
        pass

    def plot(self, aaiSummaryFile):
        # create matrix from pairwise comparisons
        matrix = defaultdict(dict)
        with open(aaiSummaryFile) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                genomeA = lineSplit[0]
                genomeB = lineSplit[2]
                aai = float(lineSplit[5])

                matrix[genomeA][genomeB] = aai
                matrix[genomeB][genomeA] = aai

        data = np.array([])
        sampleIds = matrix.keys()
        for i, sampleIdI in enumerate(sampleIds):
            row = []
            for j, sampleIdJ in enumerate(sampleIds):
                if i == j:
                    row.append(0)
                else:
                    row.append(1.0 - matrix[sampleIdI][sampleIdJ] / 100.0)

            data = np.append(data, row)

        data = np.reshape(data, (len(matrix), len(matrix)))

        mds = manifold.MDS(n_components=2, dissimilarity="precomputed")
        coords = mds.fit(data).embedding_

        # lle = manifold.LocallyLinearEmbedding(n_neighbors=5)
        # coords = lle.fit(data).embedding_
        # print lle.reconstruction_error_

        print '  Stress of metric MDS embedding: %.2f' % mds.stress_

        plt.subplots_adjust(bottom=0.1)
        plt.scatter(
            coords[:, 0], coords[:, 1], marker='o'
            )

        for label, x, y in zip(sampleIds, coords[:, 0], coords[:, 1]):
            plt.annotate(
                label,
                xy=(x, y), xytext=(-20, 20),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        plt.show()
