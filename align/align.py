# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        # Didn't use these
        # self._back_A = None
        # self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        #initialize matrices
        self._matrix_dims = (len(self._seqA)+1, len(self._seqB)+1)
        self._align_matrix = np.full(self._matrix_dims, -np.inf)
        self._gapA_matrix = np.full(self._matrix_dims, -np.inf)
        self._gapB_matrix = np.full(self._matrix_dims, -np.inf)
        self._back = np.full(self._matrix_dims, -np.inf)

        self._align_matrix[0,0] = 0
        self._back[0,0] = 0
        #top row can move left (coded as 1)
        self._back[1:,0] = 1
        #left column can move up (coded as 2)
        self._back[0,1:] = 2 

        for i in range(self._matrix_dims[0]): 
            if i == 0:
                self._gapA_matrix[i, 0] = self.gap_open
            else:
                self._gapA_matrix[i, 0] = self.gap_extend+self._gapA_matrix[i-1,0]
            
        for j in range(self._matrix_dims[1]):   
            if j == 0:
                self._gapB_matrix[0, j] = self.gap_open
            else:
                self._gapB_matrix[0, j] = self.gap_extend+self._gapB_matrix[0,j-1]   
                        

        # TODO: Implement global alignment here
        
        for i in range(1, self._matrix_dims[0]):
            for j in range(1, self._matrix_dims[1]):
                
                #calculate scores for moves
                subst_score = self.sub_dict[(self._seqA[i-1], self._seqB[j-1])]
                options = [subst_score + self._align_matrix[i-1, j-1]]
                options += [subst_score + self._gapA_matrix[i-1,j-1]]
                options += [subst_score + self._gapB_matrix[i-1,j-1]]
                self._align_matrix[i,j] = max(options)
                
                self._gapA_matrix[i, j] = max(self._align_matrix[i-1,j]+self.gap_open+self.gap_extend, 
                                        self._gapA_matrix[i-1, j]+self.gap_extend)
                self._gapB_matrix[i, j] = max(self._align_matrix[i,j-1]+self.gap_open+self.gap_extend, 
                                        self._gapB_matrix[i, j-1]+self.gap_extend)
                options = [self._align_matrix[i,j],self._gapB_matrix[i, j],self._gapA_matrix[i, j],]
                
                #add best move to traceback matrix
                self._back[i,j] = options.index(max(options))    
        		    
        return self._backtrace()


    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        self.seqA_align = ''
        self.seqB_align = ''

        _seqA_pos = self._matrix_dims[0]-1
        _seqB_pos = self._matrix_dims[1]-1

        #extract alignment score from bottom right cell
        self.alignment_score = self._align_matrix[_seqA_pos, _seqB_pos]

        while True:

            _seqA_prev_pos = _seqA_pos
            _seqB_prev_pos = _seqB_pos

            _seqA_pos, _seqB_pos = self._get_next_cell(_seqA_prev_pos, _seqB_prev_pos)

            if _seqA_pos == _seqA_prev_pos:
                self.seqA_align += '-'
            else:
                self.seqA_align += self._seqA[_seqA_pos]
            
            if _seqB_pos == _seqB_prev_pos:
                self.seqB_align += '-'
            else:
                self.seqB_align += self._seqB[_seqB_pos]
                
            if _seqA_pos == 0 and _seqB_pos == 0:
                break

        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]

        return (self.alignment_score, self.seqA_align, self.seqB_align)

    def _get_next_cell(self, i:int,j:int):
        
        """Helper for _backtrace: Returns coordinates of next cell to go to in 
        traceback according to move encoded in traceback matrix"""
        
        move = self._back[i,j]
        #move diagonally
        if move == 0:
            return (i-1, j-1)
        #move left
        if move == 1:
            return (i, j-1)
        #move up
        if move == 2:
            return (i-1, j)

def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
