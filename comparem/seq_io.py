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

import sys
import gzip


class SeqIO(object):
    def __init__(self):
        """Initialization."""
        pass

    def read_fasta(self, fasta_file):
        """Read sequences from fasta file.

        Parameters
        ----------
        fasta_file : str
            Name of fasta file to read.

        Returns
        -------
        dict : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        """
        try:
            if fasta_file.endswith('.gz'):
                open_file = gzip.open
            else:
                open_file = open

            seqs = {}
            for line in open_file(fasta_file):
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    seq_id = line[1:].split(None, 1)[0]
                    seqs[seq_id] = []
                else:
                    seqs[seq_id].append(line[0:-1])

            for seq_id, seq in seqs.iteritems():
                seqs[seq_id] = ''.join(seq)
        except:
            print  "[Error] Failed to process sequence file: " + fasta_file
            sys.exit()

        return seqs

    def read_seq(self, fasta_file):
        """Generator function to read sequences from fasta file.

        This function is intended to be used as a generator
        in order to avoid having to have large sequence files
        in memory.

        Example:
        seq_io = SeqIO()
        for seq_id, seq in seq_io.read_seq(fasta_file):
            print seq_id
            print seq

        Parameters
        ----------
        fasta_file : str
            Name of fasta file to read.

        Yields
        ------
        list : [seq_id, seq]
            Unique id of the sequence followed by the sequence itself.
        """
        try:
            if fasta_file.endswith('.gz'):
                open_file = gzip.open
            else:
                open_file = open

            seq_id = None
            seq = None
            for line in open_file(fasta_file):
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if seq_id != None:
                        yield seq_id, ''.join(seq)

                    seq_id = line[1:].split(None, 1)[0]
                    seq = []
                else:
                    seq.append(line[0:-1])

            # report last sequence
            yield seq_id, ''.join(seq)
        except:
            print  "[Error] Failed to process sequence file: " + fasta_file
            sys.exit()

    def extract_seqs(self, fasta_file, seqs_to_extract):
        """Extract specific sequences from fasta file.

        Parameters
        ----------
        fasta_file : str
            Fasta file containing sequences.
        seqs_to_extract : set
            Ids of sequences to extract.

        Returns
        -------
        dict : dict[seq_id] -> seq
            Dictionary of sequences indexed by sequence id.
        """

        seqs = {}

        for line in open(fasta_file):
            if line[0] == '>':
                seq_id = line[1:].partition(' ')[0]

                seq_of_interest = False
                if seq_id in seqs_to_extract:
                    seqs[seq_id] = []
                    seq_of_interest = True
            elif seq_of_interest:
                seqs[seq_id].append(line[0:-1])

        for seq_id, seq in seqs.iteritems():
            seqs[seq_id] = ''.join(seq)

        return seqs

    def write_fasta(self, seqs, output_file):
        """Write sequences to fasta file.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        output_file : str
            Name of fasta file to produce.
        """

        if output_file.endswith('.gz'):
            fout = gzip.open(output_file, 'wb')
        else:
            fout = open(output_file, 'w')

        for seq_id, seq in seqs.iteritems():
            fout.write('>' + seq_id + '\n')
            fout.write(seq + '\n')
        fout.close()
