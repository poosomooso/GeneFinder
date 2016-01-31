# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Serena Chen

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
stopCodons = ["TAA", "TGA", "TAG"]
startCodon = 'ATG'

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'

    Added additional tests so that all nucleotide cases are accounted for.
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    raise ValueError("get_complement(nucleotide) only supports a single character of either 'A', 'C', 'G', or 'T'.")


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    These cases go through all the nucleotides and are sufficiently un-symmetrical.
    """
    newSeq = ""
    for i in range(len(dna)-1, -1, -1):
        newSeq+=get_complement(dna[i])
    return newSeq


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGGGTTATGATTAGATGAAATAG")
    'ATGGGTTATGAT'

    Added additional test to make sure that it only checks for the first orf.
    """
    index = 0
    while index<len(dna) and not (dna[index:index+3] in stopCodons):
        index+=3
    if(index>len(dna)):
        index = len(dna)
    return dna[0:index]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TAGTAAATGTAGATGCATGAATGTAGATAGATGTGCCC")
    ['ATG', 'ATGCATGAATGTAGA', 'ATGTGCCC']

    Added additional test to make sure that it works even if the sequence doesn't start at 'ATG', 
    and also that it doesn't check for stop codons before a start codon has been found
    """
    index = 0
    orfs = []
    while index<len(dna):
        if(dna[index:index+3]==startCodon):
            o = rest_of_ORF(dna[index:]) 
            orfs.append(o)
            index+=len(o)
        else:
            index+=3
    return orfs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("A")
    []

    We haven't yet checked what happens when there are no results
    """
    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    I think these tests are sufficient
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


# def longest_ORF(dna):
#     """ Finds the longest ORF on both strands of the specified DNA and returns it
#         as a string
#     >>> longest_ORF("ATGCGAATGTAGCATCAAA")
#     'ATGCTACATTCGCAT'
#     """
#     # TODO: implement this
#     pass


# def longest_ORF_noncoding(dna, num_trials):
#     """ Computes the maximum length of the longest ORF over num_trials shuffles
#         of the specfied DNA sequence

#         dna: a DNA sequence
#         num_trials: the number of random shuffles
#         returns: the maximum length longest ORF """
#     # TODO: implement this
#     pass


# def coding_strand_to_AA(dna):
#     """ Computes the Protein encoded by a sequence of DNA.  This function
#         does not check for start and stop codons (it assumes that the input
#         DNA sequence represents an protein coding region).

#         dna: a DNA sequence represented as a string
#         returns: a string containing the sequence of amino acids encoded by the
#                  the input DNA fragment

#         >>> coding_strand_to_AA("ATGCGA")
#         'MR'
#         >>> coding_strand_to_AA("ATGCCCGCTTT")
#         'MPA'
#     """
#     # TODO: implement this
#     pass


# def gene_finder(dna):
#     """ Returns the amino acid sequences that are likely coded by the specified dna

#         dna: a DNA sequence
#         returns: a list of all amino acid sequences coded by the sequence dna.
#     """
#     # TODO: implement this
#     pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
