# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Joseph Lee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


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
    >>> get_complement('E')
    'Not a valid nucleotide'
    """
    #I put all of the explanation for the unit tests in here
    #because otherwise it thought they were part of the unit
    #tests and would return failed

    #This unit test verifies that when the function is given 'G'
    # it returns the inverse 'C'

    #This unit test verifies that when the function is given 'T'
    #it returns the inverse 'A'

    #This unit test is to verify that if something besides
    #A, T,C, or G is given as the nucleotide, the function
    #recognizes it as invalid

    if nucleotide=='A':
    	return'T'
    elif nucleotide=='C':
    	return'G'
    elif nucleotide=='T':
    	return'A'
    elif nucleotide=='G':
    	return'C'
    else:
    	return'Not a valid nucleotide'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("ATGCCCTGATCTAHAACT")
    'Not a valid dna string'
    """
    #added a unit test to protect against invalid dna strings.
    index=-1#				Here index is set to negative 1 so that it pulls the last letter in the string first
    reverse_complement=[]#	Initializing an empty list - later the complementary dna strand will be written to this
    while index>=-len(dna):
    	nucleotide=dna[index]#    This holds the character in a variable that is fed into get_complement
    	complementary_codon=get_complement(nucleotide)#		This gets the complement of the current codon in index
    	if complementary_codon=='Not a valid nucleotide':
    		return 'Not a valid dna string'
        reverse_complement.append(complementary_codon)#		This adds that codon to the reverse_complement list
        index=index-1#										This moves the index down by one (i.e. closer to the beginning)
    return ''.join(reverse_complement)#						This returns the reverse complement (.join converts list to string)
        




def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.  I added a unit test to ensure that if a stop codon
        is not detected, then the function returns the full dna strand

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGGAGATGTTAAAGTGATC")
    'ATGGAGATGTTAAAG'
    >>> rest_of_ORF("ATGCATGAATGTAGATAGATGTGCCC")
    'ATGCATGAATGTAGA'
    >>> rest_of_ORF("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
    'AAAAAAAAAAAAAAAAAAAAAAAAAAA'
    """
    stop_codon1='TAG'#		Just defining the stop codons here
    stop_codon2='TGA'
    stop_codon3='TAA'
    start_frame=0#			Setting initial values for the start and end of the first frame
    end_frame=3
    ORF_length=0#			initial value for the length of the ORF

    while ORF_length<=len(dna):
    	if dna[start_frame:end_frame]==stop_codon1 or dna[start_frame:end_frame]==stop_codon2 or dna[start_frame:end_frame]==stop_codon3:
#												If the 3 letters from start frame to end frame are one of the stop codons
    		return dna[:start_frame]#			then return the dna string up to but not including the stop codon

    	end_frame=end_frame+3#					this is outside the if loop and just indexes by 3 to look at the next 3 characters	
    	start_frame=start_frame+3#				updating the start frame by indexing it by 3
    	ORF_length=ORF_length+3#				adding 3 to the tracked ORF length for purposes of the while loop
    return dna#									this line ensures that if an in frame stop codon is never hit, it returns the whole dna strand






    #								THIS WAS THE CODE I FIRST WROTE FOR REST_OF_ORF FUNCTION, BUT BECAUSE I USED .FIND, IT WOULD LOOK FOR
    #								A STOP CODON REGARDLESS OF WHETHER IT WAS IN THE RIGHT FRAME OR NOT.  I FIXED THIS WITH THE IMPLEMENTATION
    #								ABOVE, BUT I AM KEEPING THIS CODE IN COMMENTS JUST FOR DOCUMENTATION


    #						Added 2 unit tests - the first one I added verifies that the program will recognize the TAA stop codon
    #						The second unit test I added ensures that it will recognize all stop codons in a sequence, but stop at the first
    #



  # stop_codon1='TAG'#						Here all of the stop codons have been defined
 #  stop_codon2='TGA'#						Order is arbitrary and does not impact the code
#   stop_codon3='TAA'#
    
  # stop_loc1=dna.find(stop_codon1)#		stop location 1 is the index at which stop codon 1 is found.  If not found this is set to -1
 #  stop_loc2=dna.find(stop_codon2)#		stop location 2 is the index at which stop codon 1 is found.  If not found this is set to -1
#   stop_loc3=dna.find(stop_codon3)#		stop location 3 is the index at which stop codon 1 is found.  If not found this is set to -1

    
   #if stop_codon1==stop_codon2==stop_codon3==-1:#	If none of the stop codons are found, return the entire string
    #   return dna[:]
   #else:#											Otherwise:
   #    templist=[stop_loc1,stop_loc2,stop_loc3]#	Put the locations into a list
  #     stop_loc=min(templist)#						Set the stop location to be the minumum (this way if there are multiple stop codons in a dna string it will stop at the first)
 #      if stop_loc+1%3==0:
#       	return dna[:stop_loc-3]#					Return the dna string up to but not including the stop codon (the -3 takes off the stop codon)




def find_all_ORFs_oneframe(dna,frame_number):
#
#	NOTE: I MODIFIED THIS FUNCTION TO ACCEPT FRAME NUMBER AS AN INPUT. ESSENTIALLY WHAT THIS DOES
#		  IS IT ALLOWS ME TO SPECIFY WHERE THE START OF THE ORF IS LOCATED.  FOR PURPOSES OF TESTING
#		  find_all_ORFs_oneframe I JUST SPECIFIED THE FRAME NUMBER TO BE 1 - I.E. THE FIRST FRAME
#		  WHICH STARTS AT INDEX 0
#
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC",1)
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGAAAAAAAAAAAAAAAATGTGCCC",1)
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    #the additional unit test checks to make sure the function is looking for the start codon ATG properly
    all_ORF_oneframe=[]
    start_ORF=frame_number-1#		again this is setting the start of the ORF to be the frame number minus 1...i.e. 0 for frame 1.
    length_ORFs=0

    while length_ORFs<=len(dna):
    	if dna[start_ORF:start_ORF+3]=='ATG':#				If the first 3 letters in the strand are a start codon
    		current_ORF=rest_of_ORF(dna[start_ORF: ])#		the current ORF is the result of rest_of_ORF from the start of ORF to the end of the strand
#															internal to this function is that it stops if it reaches an in frame stop codon
    		all_ORF_oneframe.append(current_ORF)#			add this ORF to the list of all the ORFs in the frame
    		start_ORF=start_ORF+len(current_ORF)#			reset the starting index to be the original start + the length of the just recorded ORF
    		length_ORFs=length_ORFs+len(current_ORF)+3#		update the length of all the ORFs combined.  The +3 is there to account for the stop codon
    	start_ORF=start_ORF+3#								this is outside the if statement...if there is no start codon, index by 3
    	length_ORFs=length_ORFs+3#							add 3 to the total length of ORFs since this length is what is being tracked by the while loop    
    return all_ORF_oneframe




    


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
    """
    first_ORF=find_all_ORFs_oneframe(dna,1)#	THIS IS WHERE MY MODIFIED find_all_ORFs_oneframe FUNCTION REALLY SHINES!!!
    second_ORF=find_all_ORFs_oneframe(dna,2)#	essentially to find all the ORFs I can literally just shift
    third_ORF=find_all_ORFs_oneframe(dna,3)#	the frame by 1
    all_ORFs=first_ORF+second_ORF+third_ORF#	and then add all of them together in one list
    return all_ORFs#							and have my function return that list.


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    complementary_strand=get_reverse_complement(dna)
    all_ORFs_both_strands=find_all_ORFs(dna)+find_all_ORFs(complementary_strand)
    return all_ORFs_both_strands

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs_both_strands=find_all_ORFs_both_strands(dna)
    longest_ORF=max(all_ORFs_both_strands,key=len)
    return longest_ORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
        >>> from load import load_seq
		>>> dna = load_seq("./data/X73525.fa")
		>>> longest_ORF_noncoding(dna,700)
		700

    """
    longest_ORFs=[]
    for i in range(num_trials):
    	random_dna=shuffle_string(dna)
    	longest_ORFs.append(longest_ORF(random_dna))
    length_longest_ORF_noncoding=len(max(longest_ORFs,key=len))
    return length_longest_ORF_noncoding


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
    """
    start_slice=0
    end_slice=3
    amino_acids=[]
    for i in range(len(dna)/3):
    	amino_acid=aa_table[dna[start_slice:end_slice]]
    	amino_acids.append(amino_acid)
    	start_slice=start_slice+3
    	end_slice=end_slice+3
    amino_acids=''.join(amino_acids)
    return amino_acids

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        >>> from load import load_seq
		>>> dna = load_seq("./data/X73525.fa")
		>>> gene_finder(dna)
		 

    """
    threshold = longest_ORF_noncoding(dna, 1500)
    all_ORFs_both_strands=find_all_ORFs_both_strands(dna)
    amino_acids_dna=[]
    for ORF in all_ORFs_both_strands:
    	if len(ORF)>threshold:
    		amino_acids_dna.append(coding_strand_to_AA(ORF))
    return amino_acids_dna


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(gene_finder, globals())
