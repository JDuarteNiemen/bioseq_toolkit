import re
from bioseq.translation import translate
from bioseq.fasta import readFASTAseq, writeFASTA

def openReadingFrame(aaseq, n=0):
    '''
    Returns the amino acid sequence from the first Methionine (start codon) and the first stop codon reached.
    Function uses a regular expression thus requires 're' to be imported.
    2 arguments, first is the sequence the second argument is which match is returned 1st, 2nd etc...
    Default is set to the first match found.
    If no match is found an empty string is produced
    '''

    # Making sequence uppercase
    aaseq = aaseq.upper()

    # creating regular expression to find amino acid sequence between Methionine up until stop codon
    pattern = r"M[A-Z]*(?=\*)"

    # Searching sequence for my regular expression
    match = re.findall(pattern, aaseq)

    # Returns the match number specified in second argument.
    if match:
        return (match[n])
    else:
        return ('')  # returns empty string if no match is found.


def candidateProtein(dnaseq):
    '''
    Returns the longest Open reading frame in a DNA sequence.
    The sequence is converted into its Amino acid sequence and the longest string between a Methionine and a stop codon is returned.
    '''
    # translating DNA to Amino acid sequence
    t1 = translate(dnaseq)
    aaseq = t1['f1']  # taking the first forward reading frame

    # creating loop to return all possible orfs
    i = 0  # Creating index to cycle through
    all_orf = []
    while True:  # loop through all possible indexes of orfs
        try:
            orf = openReadingFrame(aaseq, i)
            if not orf:
                break
            all_orf.append(orf)
            i += 1  # increment index by one
        except IndexError:
            break

    # Loop through all orf sequences and keep the longest
    longest = all_orf[0]
    for seq in all_orf:
        if len(seq) > len(longest):
            longest = seq

    if longest:
        return longest
    else:
        return ''


def maximalORF(inputfile, outputfile, proteinname):
    '''
    Takes an input of a file containing a string of DNA and outputs a fasta file in of the longest Open reading frame in the sequence.
    The function requires the use of other function found in this module: 'readFASTAseq', 'translate', 'candidateProtein', and 'writeFASTA'.

    '''

    # Reading input file
    dnaseq = readFASTAseq(inputfile)

    # Getting the longest amino acid sequence from dna sequence extracted from file
    candidate = ''
    try:
        candidate = candidateProtein(dnaseq)
    except IndexError:
        print("Error: No ORFs found (IndexError).")

    # This is the string to write to output file

    # Writing to outputfile
    writeFASTA(candidate, proteinname, outputfile)

