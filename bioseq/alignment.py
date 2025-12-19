import re

def simple_align(sequences, scoring, gap=-1):
    '''
    Performs an alignment between 2 sequences using a one step look ahead algorithm.
    Sequences are aligned character by character, being scored against each other and with the next character in the opposite sequence.
    The match returning the highest value plus the gap in indel cases is chosen for the alignment.
    Input is a tuple containing two sequences to be aligned in the form of a string
    Scoring object is required for alignment it is not created inside function. ALignment requires the Scoring class to be loaded
    Gap penalty is set to a default of -1.
    Input sequences are converted to uppercase for usage in scoring but original casing is returned in displayed alignment
    '''

    # Make input sequences upper case due to upper case nature of matrices. Original casing is returned in alignment
    upper_sequences = ()
    for seq in sequences:
        seq = seq.upper()
        upper_sequences += (seq,)

    # Set gap score penalty
    gap = gap

    score = scoring

    # grabs the sequence and then grabs an amino acid in the sequence to do scoring on
    alignment_score = 0  # initialise alignment score

    # Track position in sequences
    i_pos = 0
    j_pos = 0

    # Initialising strings for aligned sequences
    aln1 = ''
    aln2 = ''

    # Initialising boundaries
    i_boundary = len(sequences[0])
    j_boundary = len(sequences[1])

    while i_pos < i_boundary and j_pos < j_boundary:  # Ensuring position is within sequences boundaries
        # Calculating score for match
        c1 = score.match(upper_sequences[0][i_pos], upper_sequences[1][j_pos])

        # Calculating score for deletion
        if i_pos + 1 < i_boundary and j_pos + 1 < j_boundary:  # Ensuring that look ahead doesnt go past seq length
            c2 = score.match(upper_sequences[0][i_pos],
                             upper_sequences[1][j_pos + 1]) + gap  # Computing look ahead score
        else:
            c2 = gap + gap  # Look ahead score if look ahead has exceeded length of seq

        # Calculating score for insertion
        if i_pos + 1 < i_boundary and j_pos + 1 < j_boundary:  # Ensuring look ahead doesnt go past seq length
            c3 = score.match(upper_sequences[0][i_pos + 1],
                             upper_sequences[1][j_pos]) + gap  # Computing look ahead score
        else:
            c3 = gap + gap  # Look ahead score if look ahead has exceeded length of seq

        # Deciding which choice should be made based on look ahead score
        # Deciding that matches win ties with insertions and deletions. Insertions will win ties over deletions.
        # If matching current positions yields highest lookahead score
        if c1 >= c2 and c1 >= c3:
            aln1 = ''.join([aln1, sequences[0][i_pos]])  # Add current symbol to aligned sequence
            aln2 = ''.join([aln2, sequences[1][j_pos]])  # Add current symbol to aligned sequence
            alignment_score += c1  # Add match score to alignment score
            i_pos += 1  # Increase first seq position
            j_pos += 1  # Increase second seq position


        # If deletion choice yields highest lookahead score
        elif c2 > c1 and c2 > c3:
            aln1 = ''.join([aln1, '-'])  # Add - to represent deletion in aligned sequence
            aln2 = ''.join([aln2, sequences[1][j_pos]])  # Add current symbol to aligned sequence
            alignment_score += gap  # Add gap penalty to alignment score
            i_pos += 0  # Maintain first seq current position
            j_pos += 1  # Increase second seq position


        # If insertion choice yields highest lookahead score
        elif c3 > c1 and c3 >= c2:
            aln1 = ''.join([aln1, sequences[0][i_pos]])  # Add current symbol to aligned sequence
            aln2 = ''.join([aln2, '-'])  # Add - to represent insertion in aligned sequence
            alignment_score += gap  # Add gap penalty to alignment score
            i_pos += 1  # Increase first seq position
            j_pos += 0  # Maintain first seq current position

    # Finish aligning leftover positions
    while i_pos < i_boundary:
        aln1 = ''.join([aln1, sequences[0][i_pos]])
        aln2 = ''.join([aln2, '-'])
        alignment_score += gap
        i_pos += 1

    while j_pos < j_boundary:
        aln1 = ''.join([aln1, '-'])
        aln2 = ''.join([aln2, sequences[1][j_pos]])
        alignment_score += gap
        j_pos += 1

    # Create tuple of aligned sequence and alignment score
    result = (aln1, aln2, alignment_score)

    return result


def seeded_simple_align(sequences, scoring, seed, gap=-1):
    '''
    Performs a seeded alignment between 2 sequences using a one step look ahead algorithm.
    Seeds must be the same length and when gaps are removed are prefix matches of their corresponding sequences.
    Seeds must be a tuple containing sequences as type string.
    Sequences are aligned character by character, being scored against each other and with the next character in the opposite sequence.
    The match returning the highest value plus the gap in indel cases is chosen for the alignment.
    Input is a tuple containing two sequences to be aligned in the form of a string
    Scoring object is required for alignment it is not created inside function. ALignment requires the Scoring class to be loaded
    Gap penalty is set to a default of -1.
    Input sequences are converted to uppercase for usage in scoring but original casing is returned in displayed alignment
    '''

    # Setting variable name for scoring object
    score = scoring

    # Removing gaps from seed sequences
    prefixes = ()
    for seq in seed:
        seq = re.sub('-', '', seq)
        prefixes = prefixes + (seq,)

    # Checking the prefixes match the beginning of sequences.
    if prefixes[0] != sequences[0][0:len(prefixes[0])] or prefixes[1] != sequences[1][0:len(prefixes[1])]:
        raise ValueError('Seeds are not prefixes of corresponding sequences')

    # Checking seeds are the same length
    if len(seed[0]) != len(seed[1]):
        raise ValueError('Seeds have different lengths')

    # Finding length of prefixes to find start position for alignment
    start_pos_i = len(prefixes[0])
    start_pos_j = len(prefixes[1])

    # Extracting sequences from the end of the prefixes to be used in alignment
    s1 = sequences[0][start_pos_i:]
    s2 = sequences[1][start_pos_j:]
    align = (s1, s2)

    # Running simple_align on unaligned sequences to get alignment from gap
    result = simple_align(align, scoring, gap)

    # Find alignment score of seeds
    seed_score = 0
    for i in range(0, len(seed[0])):
        if seed[0][i] == '-' and seed[1][i] == '-':
            seed_score += 0

        elif seed[0][i] == '-' or seed[1][i] == '-':
            seed_score += gap

        else:
            seed_score += score.match(seed[0][i], seed[1][i])

    # Append aligned sequence to seed, and add seed score to total alignment score.
    final_result = seed[0] + result[0], seed[1] + result[1], (result[2] + seed_score)

    return final_result
