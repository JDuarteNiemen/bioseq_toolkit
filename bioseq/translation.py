def translate(dnaseq):
    '''
    Translates DNA into its corresponding amino acids in all possible reading frames.
    One argument which must be a string of DNA.
    '*' stand for stop codons. 'X' represent unknown amino acids.
    '''

    # ValueError raised if input is not a DNA string
    if not isinstance(dnaseq, str):
        raise ValueError('Argument must be a DNA sequence of type string')

    # Making sure everything is upper case
    dnaseq = dnaseq.upper()

    # Converting input string into a list
    dnaseq = list(dnaseq)

    # Creating all possible forward reading frames
    f1 = ''.join(dnaseq[0:])
    f2 = ''.join(dnaseq[1:])
    f3 = ''.join(dnaseq[2:])

    # Creating copy of sequence and reversing it
    rv_dnaseq = dnaseq.copy()
    rv_dnaseq.reverse()

    # Creating all possible reverse reading frames
    r1 = ''.join(rv_dnaseq[0:])
    r2 = ''.join(rv_dnaseq[1:])
    r3 = ''.join(rv_dnaseq[2:])

    # Initialize codon dictionary
    codon_dict = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y',
        'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
        'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
        'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S',
        'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
        'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    # Create list of frames
    frame_list = [f1, f2, f3, r1, r2, r3]
    frame_list_names = ['f1', 'f2', 'f3', 'r1', 'r2', 'r3']

    # Translating dna sequence to amino acid sequences
    translated_list = []  # Creating list for translated sequences
    for seq in frame_list:  # looping through all DNA sequences
        aaseq = ''  # Creating string for amino acid sequences
        for n in range(0, len(seq), 3):  # Moves through sequence 3 at a time
            rf = seq[n:n + 3]  # takes one codon at a time
            aaseq += codon_dict.get(rf, 'x')  # Getting the amino acid based on codon
        translated_list.append(aaseq)  # adding amino acid to list

    # zip together list of frame list names and the translated sequence
    possible_translations = dict(zip(frame_list_names, translated_list))

    return (possible_translations)