import numpy as np
from align.util import argmaxdup, fa2str, mat2array


class PairwiseAligner(object):
    '''
    Base class for allignment algorithms.

    Warning: This class should not be used directly, use child classes
    NeedlemanWunsch or SmithWaterman instead.
    '''
    def __init__(self, sim_matrix, af_gap):
        # import .mat file into amino acid-interger dict and numpy array
        self.sim_matrix_ = mat2array(sim_matrix)
        self.af_gap_ = np.abs(af_gap)
        self.score_matrix_ = np.array([])
        self.trace_matrix_ = np.array([])
        self.score_ = []
        self.align_list_ = []

    def align(self, seq_a, seq_b):
        '''
        Calculates pairwise alignments and their scores. 
        
        Parameters
        ----------
        seq_a : fasta (.fa) file or str
            An amino acid or nucleotide sequence.
        seq_b : fasta (.fa) file or str
            A second amino acid or nucleotide sequence to be aligned to the
            first sequence.
            
        Returns
        ----------
        score : lists of numpy.float64
            The calculated alignment score.
        alignment : list of strings
            A visualization of the alignment.
        '''

        # Begine Score and Traceback Procedure #
        
        # Convert fasta files to strings
        seq_a = fa2str(seq_a)
        seq_b = fa2str(seq_b)
        
        # Set amino acid-interger dict and numpy array
        am = self.sim_matrix_[0]
        sim_m = self.sim_matrix_[1]

        # Initialize score matrix and traceback matrix
        len_m, len_n = len(seq_b)+1, len(seq_a)+1
        score_matrix = np.zeros([len_m, len_n])
        trace_matrix = np.zeros([len_m, len_n])

        # Needleman-Wunsch requires first row and column to be filled with
        # gap penalties
        if self.alg_ == 'nw':
            # Filling out the gap row and gap column
            af = -self.af_gap_[0]
            for i in range(1, len_m):
                score_matrix[i, 0] = af
                trace_matrix[i, 0] = 10
                af += -self.af_gap_[1]

            af = -self.af_gap_[0]
            for j in range(1, len_n):
                score_matrix[0, j] = af
                trace_matrix[0, j] = 100
                af += -self.af_gap_[1]

        # Developing the score and traceback matrices
        
        # Encoding system for traceback. 100 = gap in sequence a, 
        # 10 = gap in sequence b, 1 = match
        enc = [1, 10, 100]
        for m, i in enumerate(seq_b):
            for n, j in enumerate(seq_a):

                # Calculate match score
                match = sim_m[am[i], am[j]]

                # Set gap scores to gap initialization penalty
                gap_t = self.af_gap_[0]
                gap_l = self.af_gap_[0]

                if trace_matrix[m, n+1] == 10:
                    gap_t = self.af_gap_[1]

                if trace_matrix[m+1, n] == 100:
                    gap_l = self.af_gap_[1]

                d = score_matrix[m, n] + match
                t = score_matrix[m, n+1] - gap_t
                l = score_matrix[m+1, n] - gap_l

                # Determine maximum candidate
                vals = [d, t, l]
                input_value = max(vals)

                # Smith-Waterman sets maximum candidate score to zero if
                # the score is less than zero.
                if self.alg_ == 'sw':
                    if input_value < 0:
                        input_value = 0

                # Select trace to be recorded
                max_index = argmaxdup(vals)
                input_trace = sum([enc[i] for i in max_index])

                # Record scores and traces in their respective matrices
                score_matrix[m+1, n+1] = input_value
                trace_matrix[m+1, n+1] = input_trace

        # Update class parameters with matrices
        self.score_matrix_ = score_matrix
        self.trace_matrix_ = trace_matrix


        # Begin Traceback Procedure #

        max_dup = 1 # count of duplicate solutions for Smith-Waterman
        fork_count = 1 # count of duplicate solutions for NW and SW

        # The traceback procedure will continue until duplicate
        # solutions are exhausted
        while fork_count > 0 or max_dup > 0:
            # Initialize allignment strings
            align_a, align_b, pair = '', '', ''

            # Traceback starting point is dependent on algorithm type
            if self.alg_ == 'nw':
                # Global alignment starts at the bottom right hand corner
                m, n = len(seq_b), len(seq_a)

            else:
                # Local alignment starts at the location of the maximum score
                val, ind, counts = np.unique(score_matrix,
                                             return_index=True,
                                             return_counts=True)
                max_dup = counts[-1]
                np.unravel_index(ind[-1], score_matrix.shape)
                m, n = np.unravel_index(ind[-1], trace_matrix.shape)
                
            self.score_.append(self.score_matrix_[m, n])

            # Traceback
            while m > 0 or n > 0 or score_matrix[m,n] > 0:
                # Loop exit for local alignment
                if self.alg_ == 'sw':
                    if score_matrix[m, n] == 0:
                        break
                if trace_matrix[m, n] - 100 >= 0:
                    if trace_matrix[m, n] - 100 != 0 and fork_count == 1:
                        fork_count += 1
                        trace_matrix[m, n] = trace_matrix[m, n] - 100
                    align_a += seq_a[n-1]
                    align_b += '-'
                    pair += ' '
                    n = n-1
                elif trace_matrix[m, n] - 10 >= 0:
                    if trace_matrix[m, n] - 10 != 0 and fork_count == 1:
                        fork_count += 1
                        trace_matrix[m, n] = trace_matrix[m, n] - 10
                    align_a += '-'
                    align_b += seq_b[m-1]
                    pair += ' '
                    m = m-1
                else:
                    align_a += seq_a[n-1]
                    align_b += seq_b[m-1]
                    pair += '|'
                    m, n = m-1, n-1
            
            # Reverse alignment strings
            align_a = align_a[::-1]
            align_b = align_b[::-1]
            pair = pair[::-1]
            
            # Create fancy alignment string
            alignment = align_a + '\n' + pair + '\n' + align_b
            self.align_list_.append(alignment)

            # For local alignment we'll remove the traceback starting point so
            # that other identical alignments will be discoverable
            if self.alg_ == 'sw':
                score_matrix[np.unravel_index(ind[-1], trace_matrix.shape)] = 0

            # Update duplicate counts
            max_dup -= 1
            fork_count -= 1
            max_dup = max(0, max_dup)
            fork_count = max(0, fork_count)
        
        return self.score_, self.align_list_


class NeedlemanWunsch(PairwiseAligner):
    '''
    Calculates global pairwise alignments and their scores between two 
    sequences using the Needleman-Wunsch algorithm. Returns the alignment 
    score and the visual representation of the alignment as a string.
    Returns global alignments with the maximum score, and multiple global
    alignments if multiple alignments have same the maximum score.

    Parameters
    ----------
    sim_matrix : .mat file (str)
        Matrix of simularity scores between nucleotides or ammino acids.
    af_gap : tuple
        Affine gap penalty constants. The first value represents the gap
        opening penalty, and the second value represents the gap extension
        penalty. Values may be input as postive or negative.

    Attributes
    ----------
    sim_matrix_ : tuple (dict, np.array)
        Matrix of simularity scores between nucleotides or ammino acids.
    af_gap_ : tuple
        Affine gap penalty constants. The first value represents the gap
        opening penalty, and the second value represents the gap extension
        penalty.
    score_matrix_ : np.array
        Matrix of allignment scores generated for the input sequences.
    trace_matrix_ : np.array
        Traceback matrix used to generate the allignment.
    score_ : lists of numpy.float64
        The calculated alignment score.
    align_list_ : list of strings
        A visualization of the alignment.
    alg_ : str
        String defining alignment algorithm type

    Notes
    -----
    .mat files may be generated using matblas from blosum50.iij

    References
    ----------
    1. Needleman, Saul B., and Christian D. Wunsch. "A general method 
    applicable to the search for similarities in the amino acid sequence of 
    two proteins." Journal of molecular biology 48.3 (1970): 443-453.

    Examples
    --------
    >>> seq_a = 'TGTTACGG'
    >>> seq_b = 'GGTTGACTA'
    >>> gap = 11, 1
    >>> sim_matrix = 'BLOSUM62.mat'
    >>> nw = NeedlemanWunsch(sim_matrix, gap)
    >>> nw.align(seq_a, seq_b)
    >>> print(nw.score_)
    [14.0]
    >>> print(nw.align_list_[0])
    TGTT-ACGG
    |||| ||||
    GGTTGACTA
    '''
    def __init__(self, sim_matrix, af_gap):
        super(NeedlemanWunsch, self).__init__(sim_matrix, af_gap)
        self.alg_ = 'nw'

class SmithWaterman(PairwiseAligner):
    '''
    Calculates local pairwise alignments and their scores between two 
    sequences using the Smith-Waterman algorithm. Returns the alignment 
    score and the visual representation of the alignment as a string. 
    Returns local alignments with the maximum score, and multiple local
    alignments if multiple alignments have the same maximum score.

    Parameters
    ----------
    sim_matrix : .mat file (str)
        Matrix of simularity scores between nucleotides or ammino acids.
    af_gap : tuple
        Affine gap penalty constants. The first value represents the gap
        opening penalty, and the second value represents the gap extension
        penalty. Values may be input as postive or negative.

    Attributes
    ----------
    sim_matrix_ : np.array
        Matrix of simularity scores between nucleotides or ammino acids.
    af_gap_ : tuple
        Affine gap penalty constants. The first value represents the gap
        opening penalty, and the second value represents the gap extension
        penalty.
    score_matrix_ : np.array
        Matrix of allignment scores generated for the input sequences.
    trace_matrix_ : np.array
        Traceback matrix used to generate the allignment.
    score_ : lists of numpy.float64
        The calculated alignment score.
    align_list_ : list of strings
        A visualization of the alignment.
    alg_ : str
        String defining alignment algorithm type

    Notes
    -----
    .mat files may be generated using matblas from blosum50.iij
    
    References
    ----------
    1. Smith, Temple F., and Michael S. Waterman. 
    "Identification of common molecular subsequences." 
    Journal of molecular biology 147.1 (1981): 195-197.

    Examples
    --------
    >>> seq_a = 'TGTTACGG'
    >>> seq_b = 'GGTTGACTA'
    >>> gap = 11, 1
    >>> sim_matrix = 'BLOSUM62.mat'
    >>> sw = SmithWaterman(sim_matrix, gap)
    >>> sw.align(seq_a, seq_b)
    >>> print(sw.score_)
    [18.0]
    >>> print(sw.align_list_[0])
    GTT-AC
    ||| ||
    GTTGAC
    '''
    def __init__(self, sim_matrix, af_gap):
        super(SmithWaterman, self).__init__(sim_matrix, af_gap)
        self.alg_ = 'sw'