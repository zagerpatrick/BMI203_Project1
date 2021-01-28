# Project 1 - Sequence Alignment
![BuildStatus](https://github.com/zagerpatrick/BMI203_Project1/workflows/HW1/badge.svg?event=push)

Code should be run from the top level of the directory, see Jupyter notebook and function docstrings for examples of usage.

# API Documentation

<a name="notebook_util"></a>
# notebook\_util

<a name="align"></a>
# align

<a name="align.util"></a>
# align.util

<a name="align.util.argmaxdup"></a>
#### argmaxdup

```python
argmaxdup(a)
```

Returns the indices of the maximum values, returing multiple indices if
there are duplicate max values.

Parameters
----------
a : list
    input list

Returns
----------
maxdup: list
    list of indices of maximum values.

<a name="align.util.fa2str"></a>
#### fa2str

```python
fa2str(file)
```

Returns string of the sequence of the fasta file.

Parameters
----------
file : fasta (.fa) file
    An amino acid or nucleotide sequence.

Returns
----------
amino: str
    An amino acid or nucleotide sequence.

<a name="align.util.mat2array"></a>
#### mat2array

```python
mat2array(file)
```

Returns dict and np.array of the .mat file.

Parameters
----------
file : (.mat) file
    Matrix of simularity scores between nucleotides or ammino acids.

Returns
----------
amino: tuple (dict, np.array)
    Matrix of simularity scores between nucleotides or ammino acids.

<a name="align.algs"></a>
# align.algs

<a name="align.algs.PairwiseAligner"></a>
## PairwiseAligner Objects

```python
class PairwiseAligner(object)
```
Base class for allignment algorithms.
Warning: This class should not be used directly, use child classes
NeedlemanWunsch or SmithWaterman instead.
<a name="align.algs.PairwiseAligner.align"></a>
#### align
```python
 | align(seq_a, seq_b)
```

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

<a name="align.algs.NeedlemanWunsch"></a>
## NeedlemanWunsch Objects

```python
class NeedlemanWunsch(PairwiseAligner)
```
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
<a name="align.algs.SmithWaterman"></a>
## SmithWaterman Objects

```python
class SmithWaterman(PairwiseAligner)
```
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
