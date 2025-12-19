from bioseq.fasta import readFASTAseq

def AAtypes(aaseq):
    '''
    Calculates the proportion of polar, small and hydrophobic amino acids found in the input string.
    Output is a tuple in polar, small, hydrophobic order.
    '''
    # Ensuring all characters are upper case
    aaseq = aaseq.upper()

    # Dictionary of
    aa_classes = {'polar': ['W', 'H', 'K', 'R', 'Y', 'T', 'C', 'S', 'N', 'Q', 'D', 'E'],
                  'small': ['D', 'N', 'T', 'S', 'C', 'A', 'G', 'P', 'V'],
                  'hydrophobic': ['I', 'V', 'L', 'M', 'F', 'Y', 'W', 'H', 'K', 'T', 'C', 'A']}

    try:
        counts = {}
        for k in aa_classes:
            count = 0
            for a in aaseq:
                if a in aa_classes[k]:
                    count = count + 1
            proportion = round(count / len(aaseq), 3)
            counts[k] = proportion

        stats = tuple(counts.values())
        return stats
    except ZeroDivisionError as e:
        print('Error:', e, '| Empty input, please enter a sequence.')


def AAtypetable(filelist, outputfile):
    '''
    Takes a list of fasta files containing amino acid sequence data. Sequence is read and the proportion of polar, small and hydrophobic residues are calculated.
    Results are outputted to a file in CSV format.
    Requires the AAtypes() function.
    '''

    # Creating header format
    header = '#' + 'Filename' + ',' + 'Polar' + ',' + 'Small' + ',' + 'Hydrophobic' + '\n'

    # Creating output file and opening it for writing
    with open(outputfile, 'wt') as OUTF:
        OUTF.write(header)

        # Looping through files to get stats and writing tuples to the output file
        for file in filelist:  # Looping through file names in the list
            try:
                seq = readFASTAseq(file)
                stats = AAtypes(seq)  # Getting proportions of amino acid types
                for i in range(0, len(stats)):
                    data = ','.join([str(i) for i in stats])
                line_format = file + ',' + data + '\n'
                OUTF.write(line_format)

            except FileNotFoundError as e:
                print('Error:', e, '|', 'File does not exist in this directory')



