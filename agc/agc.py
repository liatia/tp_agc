#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Lara Herrmann"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lara Herrmann"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lara Herrmann"
__email__ = "lara.herrma@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, 
                        required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """ Lecture du fichier fastq en entree
        Retourne : générateur de séquences de longueur l >= minseqlen: yield sequence
    """
    if amplicon_file.endswith(".gz"):
        with gzip.open(amplicon_file, "rb") as filin:
            sequence = b''
            for line in filin:
                if line.startswith(b">"):
                    if (len(sequence) >= minseqlen):
                        yield sequence.decode('ascii')
                    sequence = b''
                else:
                    sequence += line.strip()
        yield sequence.decode('ascii')
    else:
        with open(amplicon_file, "r") as filin:
            sequence = ''
            for line in filin:
                if line.startswith(">"):
                    if (len(sequence) >= minseqlen):
                        yield sequence
                    sequence = ''
                else:
                    sequence += line.strip()
        yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    occu_dict = {}
    del_list = []
    for sequence in list(read_fasta(amplicon_file, minseqlen)):
        if not sequence in occu_dict:
            occu_dict[sequence] = 0
        occu_dict[sequence] += 1
    for sequence in occu_dict:
        if occu_dict[sequence] < mincount:
            del_list.append(sequence)
    for sequence in sorted(occu_dict, key=occu_dict.get, reverse=True):
        if sequence in del_list:
            del occu_dict[sequence]
        else:
            yield(sequence, occu_dict[sequence])


def get_chunks(sequence, chunk_size):
    print("nouvel appel")
    chunk_list = []
    start = 0
    end = chunk_size
    while end < len(sequence):
        chunk_list.append(sequence[start:end])
        print(len(sequence[start:end]))
        start += chunk_size
        end += chunk_size
    if len(chunk_list) < 4 :
        raise(ValueError("Less than 4 chunks"))
    return chunk_list


def cut_kmer(sequence, kmer_size):
    """ Coupe les sequences en kmers
        Retourne : tous les kmers trouves dans les sequences
    """
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i + kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for kmer in list(cut_kmer(sequence, kmer_size)):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    indentic_nucl = 0
    for nt in zip(alignment_list[0], alignment_list[1]):
        if nt[0] == nt[1]:
            indentic_nucl += 1
    alignemnt_length = len(alignment_list[0])
    return indentic_nucl/alignemnt_length * 100


def detect_chimera(perc_identity_matrix):
    pass


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


def write_OTU(OTU_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
