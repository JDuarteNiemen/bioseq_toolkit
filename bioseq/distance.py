import math
from bioseq.fasta import readFASTAseq


# Creating my exception to raise if tuples are different length
class DimensionalityException(Exception):
    pass


def distance(veca, vecb, metric):
    '''
    Computes the distance between 2 vectors. 2 different methods, Euclidean distance and Manhattan distance.
    Metric argument determines the method used.
    'l1'=Euclidean distance. 'l2'=Manhattan distance.
    Input requires 2 tuples with equal length containing float values.
    '''

    # Check tuples are the same length
    if len(veca) != len(vecb):
        raise DimensionalityException

    # Check tuples are not empty
    if len(veca) == 0 or len(vecb) == 0:
        raise DimensionalityException

    if metric == 'l1':
        # Calculating euclidean distance
        roll = 0
        for i in range(0, len(veca)):
            euclid = (veca[i] - vecb[i]) ** 2
            roll += euclid
        euclid_distance = math.sqrt(roll)
        return round(euclid_distance, 3)

    elif metric == 'l2':
        # calculating manhattan distance
        manhattan_distance = 0
        for i in range(0, len(veca)):
            manhattan = abs((veca[i] - vecb[i]))
            manhattan_distance += manhattan
        return round(manhattan_distance, 3)

    else:
        raise ValueError


def readTable(filename):
    '''
    Takes the filename of a CSV file containing proportions of small, polar and hydrophobic amino acids and stores them in a dictionary.
    Filename stored as keys and the values are a tuple containing the floats from calculations.
    '''

    out = readFASTAseq(filename, '#', 'ORIGINAL', 'KEEP')  # Using readFASTAseq function with required arguments

    data = {}  # initialising dictionary for csv to be stored
    values = []  # initialising list for values to be stored
    lines = out.strip().split('\n')  # stripping special characters
    for line in lines:
        values.append(line.split(','))  # adding lines as list to

    # Convert values in tuple to float
    for i in range(0, len(values)):  # moving through list of lines
        for x in range(1, 4):  # moving through indexes in tuple
            values[i][x] = float(values[i][x])
        data[values[i][0]] = tuple(values[i][1:4])  # Adding file names as keys and calculations as values

    return data


def distanceMatrix(inputfile, outputfile, metric):
    '''
    Function takes a text file in csv format and outputs a tsv distance matrix file (.dmf)
    Calculates the distance between the proportions of amino acids for each protein using the given metric.
    'l1'=Euclidean distance, 'l2'=Manhattan distance.
    Requires the distance() function
    '''

    data = readTable(inputfile)

    list_keys = list(data.keys())
    list_values = list(data.values())

    # Create header format
    files = '\t'.join(list_keys)
    header = '# filename\t' + files + '\n'

    with open(outputfile + '.dmf', 'wt') as OUTF:
        # writing header
        OUTF.write(header)

    # Calculating values based on provided metric and then writing rows to file
    for x in range(0, len(list_values)):  # looping through each tuple
        line = ''
        for y in range(0, len(list_values)):  # looping through each tuple again to create matrix
            answer = distance(list_values[x], list_values[y], metric)  # Calculating distance using distance function
            line += '\t' + f"{answer:.3f}" + '\t'
        with open(outputfile + '.dmf', 'at') as OUTF:
            OUTF.write(str(list_keys[x]) + line + '\n')


