import sys
sys.path.insert(1,'/Users/patrick/Projects/BMI-203/Project1')

import pytest
import numpy as np
import align.util
import align.algs

data_path = '/Users/patrick/Projects/BMI-203/Project1/test_data/'

def test_fasta_io():
    # Setup
    fa_file = 'test_fasta_io.fa'
    
    desired = 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPE' \
              'MAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSY' \
              'GGDEGAWTAVAGALMGEIEPDM'

    # Exercise
    actual = align.util.fa2str(data_path + fa_file)

    # Verify
    assert actual == desired

def test_scoring_matrix_io():
    # Setup
    matrix_file = 'BLOSUM62.mat'
    
    desired_dict = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 
                   'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 
                   'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 
                   'V': 19, 'B': 20, 'Z': 21, 'X': 22, '*': 23}
    desired_array_5x5 = np.array([[ 4, -1, -2, -2,  0],
                                  [-1,  5,  0, -2, -3],
                                  [-2,  0,  6,  1, -3],
                                  [-2, -2,  1,  6, -3],
                                  [ 0, -3, -3, -3,  9]])

    # Exercise
    actual_dict = align.util.mat2array(data_path + matrix_file)[0]
    actual_array = align.util.mat2array(data_path + matrix_file)[1]
    actual_array_5x5 = actual_array[0:5, 0:5]

    # Verify
    assert actual_dict == desired_dict
    np.testing.assert_array_equal(actual_array_5x5, desired_array_5x5)

def test_identical():
    # Setup
    seq_a = 'test_identical_1.fa'
    seq_b = 'test_identical_2.fa'
    gap = 11, 1
    sim_matrix = 'BLOSUM62.mat'
    
    desired = [382]

    # Exercise
    sw = align.algs.SmithWaterman(data_path + sim_matrix, gap)
    sw.align(data_path + seq_a, data_path +  seq_b)
    sw_actual = sw.score_

    nw = align.algs.NeedlemanWunsch(data_path + sim_matrix, gap)
    nw.align(data_path + seq_a, data_path +  seq_b)
    nw_actual = nw.score_

    # Verify
    assert sw_actual == desired
    assert nw_actual == desired

def test_alignment_score():
    # Setup
    
    # Exercise
    
    # Verify
	assert True
