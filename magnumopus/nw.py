import numpy as np
from typing import Tuple

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> Tuple[Tuple[str, str], int]:
    """Align two sequences using the Needleman-Wunsch algorithm

    Args:
        seq_a: The first sequence
        seq_b: The second sequence
        match: match score to use in alignment
        mismatch: mismatch penalty to use in alignment
        gap: gap penalty to use in alignment
    
    Return:
        (aligned_seq_a, aligned_seq_b): The optimal alignment of two sequences
        alignment_score: Score of the alignment
    """
    #Lengths of input sequences
    len_a = len(seq_a)
    len_b = len(seq_b)

    #Initialize the score matrix with gap penalties in the top row and left column
    score_matrix = np.zeros((len_a + 1, len_b + 1))
    score_matrix[:,0] = np.linspace(0, len_a * gap, len_a + 1)
    score_matrix[0,:] = np.linspace(0, len_b * gap, len_b +1)

    #Initialize the pointer matrix
    pointer_matrix = np.zeros((len_a + 1, len_b + 1))
    pointer_matrix[:,0] = 3
    pointer_matrix[0,:] = 4

    #Create temporary scores
    t = np.zeros(3)
    for i in range(len_a):
        for j in range(len_b):
            if seq_a[i] == seq_b[j]:
                t[0] = score_matrix[i,j] + match
            else:
                t[0] = score_matrix[i,j] + mismatch
            t[1] = score_matrix[i,j+1] + gap
            t[2] = score_matrix[i+1,j] + gap
            tmax = np.max(t)
            #Filling the matrix
            score_matrix[i+1,j+1] = tmax
            if t[0] == tmax:
                pointer_matrix[i+1,j+1] += 2 #Diagonal
            if t[1] == tmax:
                pointer_matrix[i+1,j+1] += 3 #Vertical
            if t[2] == tmax:
                pointer_matrix[i+1,j+1] +=4 #Horizontal
    
    #Traceback through optimal alignment
    aligned_seq_a = []
    aligned_seq_b = []
    i, j = len_a, len_b
    while i>0 or j>0:
        if pointer_matrix[i,j] in [2,5,6,9]: #Diagonal move
            aligned_seq_a.append(seq_a[i-1])
            aligned_seq_b.append(seq_b[j-1])
            i -= 1
            j -= 1
        elif pointer_matrix[i,j] in [3,5,7,9]: #Vertical move
            aligned_seq_a.append(seq_a[i-1])
            aligned_seq_b.append('-')
            i -= 1
        elif pointer_matrix[i,j] in [4,6,7,9]: #Horizontal move
            aligned_seq_a.append('-')
            aligned_seq_b.append(seq_b[j-1])
            j -= 1
    
    #Reverse the aligned sequences
    aligned_seq_a = ''.join(aligned_seq_a)[::-1]
    aligned_seq_b = ''.join(aligned_seq_b)[::-1]

    #Calculate the alignment score (botton-right cell of the score matrix)
    alignment_score = int(score_matrix[len_a][len_b])

    return (aligned_seq_a, aligned_seq_b), alignment_score
