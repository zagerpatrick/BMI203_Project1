import numpy as np


def argmaxdup(a):
    '''
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
    '''

    if len(a) == 0:
        return []
    maxdup = [0]
    max_ = a[0]
    for i in range(1, len(a)):
        if a[i] > max_:
            maxdup = [i]
            max_ = a[i]
        elif a[i] == max_:
            maxdup.append(i)
    return maxdup


def fa2str(file):
    '''
    Returns string of the sequence of the fasta file.

    Parameters
    ----------
    file : fasta (.fa) file
        An amino acid or nucleotide sequence.

    Returns
    ----------
    amino: str
        An amino acid or nucleotide sequence.
    '''

    amino = ''
    with open(file) as file:
        for n, line in enumerate(file):
            if n == 0:
                pass
            else:
                amino += line[:-1]
    return amino


def mat2array(file):
    '''
    Returns dict and np.array of the .mat file.

    Parameters
    ----------
    file : (.mat) file
        Matrix of simularity scores between nucleotides or ammino acids.

    Returns
    ----------
    amino: tuple (dict, np.array)
        Matrix of simularity scores between nucleotides or ammino acids.
    '''

    i = 0
    n = 0
    amino = {}
    sim = []
    with open(file) as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                if i == 0:
                    for j in line:
                        if j == ' ' or j == '\n':
                            pass
                        else:
                            amino[j] = n
                            n += 1
                    i += 1
                else:
                    sim_l = []
                    for j in line.split():
                        sim_l.append(int(j))
                    sim.append(sim_l)
        sim = np.array(sim)
    return amino, sim