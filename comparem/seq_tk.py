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

"""
*****************************************************************************
To do:
 - this class still needs a lot of work
 - anything for calculating sequence statistics should be put in here
*****************************************************************************
"""


class SeqToolkit(object):
    def __init__(self):
        """Initialization."""
        pass

    def gc(self, seq):
        """Calculate GC content of a sequence.

        GC is calculated as (G+C)/(A+C+G+T), where
        each of these terms represents the number
        of nucleotides within the sequence. Ambiguous
        and degenerate bases are ignored. Uracil (U)
        is treated as a thymine (T).

        Parameters
        ----------
        seq : str
            Nucleotide sequence.

        Returns
        -------
        float
            GC content of sequence.
        """

        s = seq.upper()
        a = s.count('A')
        c = s.count('C')
        g = s.count('G')
        t = s.count('T') + s.count('U')

        return float(g + c) / (a + c + g + t)

    def ambiguous_nucleotides(self, seq):
        """Count ambiguous or degenerate nucleotides in a sequence.

        Any base that is not a A, C, G, or T/U is considered
        to be ambiguous or degenerate.

        Parameters
        ----------
        seq : str
            Nucleotide sequence.

        Returns
        -------
        int
            Number of ambiguous and degenerate bases.
        """

        s = seq.upper()
        a = s.count('A')
        c = s.count('C')
        g = s.count('G')
        t = s.count('T') + s.count('U')

        return len(seq) - (a + c + g + t)

    def N50(self, seqs):
        """Calculate N50 for a set of sequences.

         N50 is defined as the length of the longest
         sequence, L, for which 50% of the total bases
         are present in sequences of length >= L.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence ids.

        Returns
        -------
        int
            N50 for the set of sequences.
        """

        seq_lens = [len(x) for x in seqs.values()]
        threshold = sum(seq_lens) / 2.0

        seq_lens.sort(reverse=True)

        current_sum = 0
        for seq_len in seq_lens:
            current_sum += seq_len
            if current_sum >= threshold:
                N50 = seq_len
                break

        return N50

    def mean_length(self, seqs):
        """Calculate mean length of sequences.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence ids.

        Returns
        -------
        int
            Mean length of sequences.
        """

        """***TO DO***"""
        pass

    def max_length(self, seqs):
        """Identify longest sequence.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence ids.

        Returns
        -------
        int
            Length of longest sequence.
        """

        """***TO DO***"""
        pass

    def coding_density(self):
        pass

    def extract_contigs(self, seq, num_ambiguous_bases=10):
        pass
