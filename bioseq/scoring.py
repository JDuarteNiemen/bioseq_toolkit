from Bio.Align import substitution_matrices

class Scoring:
    '''
    Class loads a scoring matrix from the Biopython subsitituiton_matrices, name must be provided in the argument as a string.
    It also accepts a custom substition matrix in type dictionary.
    Custom matrices must be in the format of {('res_1', 'res_2'): <score>, ('res_3', 'res_4'): <score>}. Where res_n stands for a 1-letter aminoacid code.
    Empty items in dictionary will default to a identity matrix. match=1, mismatch=0
    One attribute is matrix described above.
    One method, 'match('res1', 'res2')'. This returns the match score of the two residues using the given matrix
    '''

    def __init__(self, matrix):
        if isinstance(matrix, str):
            # Getting names of avalible matrices
            avaliable = substitution_matrices.load()

            #
            if matrix in avaliable:
                self.matrix = substitution_matrices.load(matrix)
            else:
                raise ValueError(f"Matrix '{matrix}' not available in Biopython. Please select one of {avaliable}")

        elif isinstance(matrix, dict):
            self.matrix = matrix

    def match(self, res1, res2):

        query = (res1, res2)
        mirror = (res2, res1)

        # Score for query pair
        if query in self.matrix:
            return self.matrix[query]

        # Find scoring for mirror pair
        if mirror in self.matrix:
            return self.matrix[mirror]

        # Residue pair not in matrix
        if query not in self.matrix:
            if res1 == res2:  # Residue pair is the same
                return 1
            return 0