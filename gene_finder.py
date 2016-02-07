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
    #iterate backwords through the string
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


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("")
    ''

    Checking to make sure that it can handle no ORFs.
    """
    try:
        return max(find_all_ORFs_both_strands(dna), key=lambda s:len(s))
    except ValueError:
        #if empty list
        return ""



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    long_orf_strings = []
    for i in range(num_trials):
        long_orf_strings.append(longest_ORF(shuffle_string(dna)))
    return max(long_orf_strings, key=lambda s:len(s))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCCCGCTTTTT")
        'MPAF'
        >>> coding_strand_to_AA("A")
        ''

        testing len(dna)%3==1, and an string with only one nucleotide
    """
    aminos = ""
    for i in range(0,len(dna)-2,3): 
        #len - 2 because if the possible lengths:
        #len%3==0, len%3 == 1, len%3==2
            #len-2 allows you to access the first index of the last full group of 3
            #stopping before the leftover group of 1 or 2
            #sorry for jank implementation
        aminos+=aa_table[dna[i:i+3]]
    return aminos


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna,1500)
    all_orfs = find_all_ORFs_both_strands(dna)
    aminos = []
    for o in all_orfs:
        if len(o)>=len(threshold):
            aminos.append(coding_strand_to_AA(o))
    return aminos
def nitrogenase_substring():
    """ finds the metagenome with the longest substring in common with
    the nitrogenase sequence.

    I have no idea if it works, but it spit out an answer with no errors after over an hour
    of running, so i'm considering that success.

    """
    import load
    nit_seq = load.load_nitrogenase_seq()
    metagenomes = load.load_metagenome()
    subs = []
    for meta in metagenomes:
        subs.append((meta[0], longest_substring(meta[1], nit_seq)))
    return max(subs, key=lambda s:s[1])[0]




def longest_substring(str1, str2):
    """
    Finds the longest substring in common in two given strings.

    >>> longest_substring("ABCDEFGABACD", "FEDCBAABAC")
    'ABAC'
    """
    zeros = [0]*(len(str2)+1)
    L = [zeros]
    for i in range(len(str1)+1):
        L.append(zeros[:])

    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if str1[i-1] == str2[j-1]:
                L[i][j] = L[i-1][j-1] + 1
            else:
                L[i][j] = 0
    maxr = 0
    maxc = 0
    maxnum = 0
    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if L[i][j] > maxnum:
                maxnum = L[i][j]
                maxr = i
                maxc = j
    return str1[maxr-maxnum:maxr]

if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    
    print nitrogenase_substring()

    

