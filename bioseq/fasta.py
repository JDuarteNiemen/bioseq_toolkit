def readFASTAseq(fastafile, header='>', case='UPPER', strip='STRIP'):
    '''
    Reads a fasta file and returns the sequence found in the fasta file as a string.

    Second argument specifies the starting character of the header, the default is '>' typically used in fasta files.
    This function is compatible for various file types to exclude headers. Any value can be put which will exclude any line beginning with iput value.

    Third argument gives user the option to output everything as upper case, lower case or to keep it in the original casing. Default is set to output as upper case.
    There are three options: 'UPPER', 'LOWER' and 'original'

    Fourth arugument keeps special characters that specify line breaks or tab breaks etc. Default is set to remove them
    There are 2 options: 'STRIP' and 'KEEP'.
    '''
    # Raising error if arguments do not match options avaliable
    if case not in ['UPPER', 'LOWER', 'ORIGINAL']:
        raise ValueError
    if strip not in ['STRIP', 'KEEP']:
        raise ValueError

    # Opening file for reading
    with open(fastafile, 'rt') as INFILE:
        seq = ''
        for line in INFILE:
            if not line.startswith(header):  # skipping headers
                if strip == 'STRIP':
                    seq += line.rstrip()  # removes special characters
                else:
                    seq += line

        if case == 'UPPER':  # making sequence upper case
            seq = seq.upper()

        elif case == 'LOWER':  # making sequence lower case
            seq = seq.lower()

        if case == 'original':  # keeping sequence in original format
            seq = seq

    return seq


def writeFASTA(sequence, description, filename):
    '''
    Writes a description line as a header, then sequence in the body of the file, breaking every 60 lines to a specified filename.
    '''

    # Setting header format
    header = '>' + description + '\n'
    # Setting sequence format
    linewidth = 60
    pos = 0

    # Opening file in writing format
    with open(filename, 'wt') as OUTF:
        # Writing header to file
        OUTF.write(header)

        # Writing sequence to file using format
        while pos < len(sequence):
            OUTF.write(sequence[pos:pos + linewidth] + '\n')
            pos = pos + linewidth
