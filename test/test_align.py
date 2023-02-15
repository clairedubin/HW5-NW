# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    t = NeedlemanWunsch(sub_matrix_file='substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    score, seq1_algn, seq2_algn = t.align(seq1, seq2)

    correct_align_matrix =  np.array([[  0., -np.inf, -np.inf, -np.inf],
                    [-np.inf,   5., -11., -13.],
                    [-np.inf, -12.,   4.,  -8.],
                    [-np.inf, -12.,  -1.,   5.],
                    [-np.inf, -14.,  -6.,   4.]])

    assert np.array_equal(t._align_matrix, correct_align_matrix), "Align_matrix not as expected"
    assert seq1_algn == "MYQR", "Seq1 does not match expected sequence"
    assert seq2_algn == "M-QR", "Seq2 does not match expected sequence"

    assert t._gapA_matrix[0,0] == t.gap_open, "[0,0] cell of _gapA_matrix does not equal gap open score"
    assert t._gapB_matrix[0,0] == t.gap_open, "[0,0] cell of _gapB_matrix does not equal gap open score"


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    t = NeedlemanWunsch(sub_matrix_file='substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    score, seq3_algn, seq4_algn = t.align(seq3, seq4)

    assert score == 17, "Score does not match correct score of 17"
    assert seq3_algn == "MAVHQLIRRP", "Aligned seq3 does not match correct sequence"
    assert seq4_algn == "M---QLIRHP", "Aligned seq4 does not match correct sequence"



