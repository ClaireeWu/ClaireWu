def is_base_pair(ch1, ch2):
    """(str, str) -> bool
    
    Precondition: ch1 and ch2 are both one single character representing a 
    base('A', 'T', 'C', or 'G')
    
    Return True if and only if ch1 and ch2 form a base pair
    
    >>> is_base_pair('A', 'C')
    False
    >>> is_base_pair('T', 'T')
    False
    >>> is_base_pair('A', 'T')
    True
    """
    
    pair1 = 'AT'
    pair2 = 'CG'
    
    if ch1 != ch2:
        if ch1 in pair1 and ch2 in pair1:
            return True
        elif ch1 in pair2 and ch2 in pair2:
            return True
    return False
        
def is_dna(strand1, strand2):
    """(str, str) -> bool
    
    Precondition: strand1 and strand2 are of equal length and only contain the 
    four characters that represent bases
    
    Return True if and only if strand1 and strand2 form a properly base-paired 
    DNA molecule.
    
    >>> is_dna('AT', 'CT')
    False
    >>> is_dna('AT', 'TA')
    True
    >>> is_dna('ATCGG', 'TAGCC')
    True
    """
    
    for i in range(len(strand1)):
        if not is_base_pair(strand1[i],strand2[i]):
            return False
    return True

def is_dna_palindrome(strand1, strand2):
    """(str, str) -> bool
    
    Return True if and only if strand1 and strand2 form a DNA palindrome
    
    >>> is_dna_palindrome('AT', 'CT')
    False
    >>> is_dna_palindrome('ATCG', 'ATCG')
    False
    >>> is_dna_palindrome('AT', 'TA')
    True
    """
    
    total = strand1 + strand2
    
    for i in range(len(total) // 2):
        if total[i] != total[-1 - i]:
            return False
    return True    
    
def restriction_sites(strand, recognition_sequence):
    """(str, str) -> list of int
        
    Return a list of all the indices where the sequence appears in the strand
        
    >>> restriction_sites('ACCTTTTC', 'GAATTC')
    []
    >>> restriction_sites('GCGGCCGCTTTC', 'GCGGCCGC')
    [0]
    >>> restriction_sites('AGCTAGCTAGCT', 'AGCT')
    [0, 4, 8]
    """
    
    sites = []
    
    if strand.find(recognition_sequence) == -1:
        return sites
    sites.append(strand.find(recognition_sequence))
    i = strand.find(recognition_sequence)
    while strand.find(recognition_sequence, i + 1) != -1:
        sites.append(strand.find(recognition_sequence, i + 1))
        i = strand.find(recognition_sequence, i + 1)
    return sites

def match_enzymes(strand, enzymes, recognition_sequences):
    """(str, list of str, list of str) ->list of two-item[str, list of int]lists
    
    Return a list of two-item lists where the first item of each two-item list 
    is the enzymes and the second item is the recognition_sequences of the 
    strand that the enzymes cuts
    
    >>> match_enzymes('CCCCATC', ['EcoRV'], ['GATATC'])
    [['EcoRV', []]]
    >>> match_enzymes('TCTAGATCTAGATCTAGA', ['XbaI'], ['TCTAGA'])
    [['XbaI', [0, 6, 12]]]
    >>> match_enzymes('TCTAGATCTAGATCTAGA', ['AluI', 'TaqI'], ['AGCT', 'TCGA'])
    [['AluI', []], ['TaqI', []]]
    """
    
    matched = []
    
    for i in range(len(recognition_sequences)):
        current_enzymes = []
        current_sites = restriction_sites(strand, recognition_sequences[i])
        current_enzymes.append(enzymes[i])
        current_enzymes.append(current_sites)
        matched.append(current_enzymes)
    return matched
          
def one_cutters(strand, enzymes, recognition_sequences):
    """(str, list of str, list of str) ->list of two-item[str, int] lists
    
    Return a list of two-item lists representing the 1-cutters for the DNA 
    strand.The first item of each two-item list is enzymes and the second item 
    is the recognition_sequences.
    
    >>> one_cutters('TCTAGATCTAGATCTAGA', ['XbaI', 'HI'], ['TCTAGA', 'GAT'])
    []
    >>> one_cutters('AAAGATTCTAGAT', ['XbaI', 'HI'], ['TCTAGA', 'GAT'])
    [['XbaI', 6]]
    >>> one_cutters('AGCTTCGATHGACGC', ['TaqI', 'AluI', 'HgaI'], ['TCGA', 
    'AGCT', 'GACGC'])
    [['TaqI', 4], ['AluI', 0], ['HgaI', 10]]
    """
    
    all_cutters = match_enzymes(strand, enzymes, recognition_sequences)
    single_cutters = []
    
    for i in all_cutters:
        if len(i[1]) == 1:
            current = []
            current.append(i[0])
            current.append(i[1][0])
            single_cutters.append(current)
    return single_cutters

def correct_mutations(mutated_strands, clean, enzyme_names, sequences):
    """(list of str, str, list of str, list of str) -> NoneType
    
    Modifies mutated_strands that share a 1-cutter with clean by replacing all 
    bases starting at the 1-cutter in mutated_strand with all bases starting at 
    the 1-cutter in clean, up to and including the end of the strand.
    
    >>> mutated_strands = ['ACGTGGCCTAGCT', 'CAGCTGATCG']
    >>> clean = 'ACGGCCTT'
    >>> enzyme_names = ['HaeIII', 'HgaI', 'AluI']
    >>> sequences = ['GGCC', 'GACGC', 'AGCT']
    >>> correct_mutations(mutated_strands, clean, enzyme_names, sequences)
    ['ACGTGGCCTT', 'CAGCTGATCG']
    >>> mutated_strands = ['ACGTGGCCTAGCT', 'CAGCTGACGCTTAA']
    >>> clean = 'ATTTGACGCTTTT'
    >>> enzyme_names = ['HaeIII', 'HgaI', 'AluI']
    >>> sequences = ['GGCC', 'GACGC', 'AGCT']
    >>> correct_mutations(mutated_strands, clean, enzyme_names, sequences)
    ['ACGTGGCCTAGCT', 'CAGCTGACGCTTTT']
    >>> mutated_strands = ['ACGTGTCGCGCTA', 'ACAGCTCTAGACTTAA']
    >>> clean = 'ATTTGACGCTT'
    >>> enzyme_names = ['SalI', 'HgaI', 'XbaI']
    >>> sequences = ['GTCGAC', 'GACGC', 'TCTAGA']
    >>> correct_mutations(mutated_strands, clean, enzyme_names, sequences)
    ['ACGTGTCGCGCTA', 'ACAGCTCTAGACTTAA']
    """

    cutter = one_cutters(clean, enzyme_names, sequences)[0]
    
    for i in range(len(enzyme_names)):
        if enzyme_names[i] == cutter[0]:
            index = i
            
    replacing_sequences = clean[clean.find(sequences[index]):]
    
    for i in range(len(mutated_strands)):
        if sequences[index] in mutated_strands[i]:
            cut_index = mutated_strands[i].find(sequences[index])
            mutated_strands[i] = mutated_strands[i][:cut_index] + \
            replacing_sequences
            

    
