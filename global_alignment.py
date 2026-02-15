import numpy as np

INF = 10**9


# ---------- FASTA READER ----------
def read_fasta(file):
    seq = ""
    with open(file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


# ---------- GLOBAL ALIGNMENT WITH AFFINE GAPS ----------
def affine_global_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):

    n, m = len(seq1), len(seq2)

    M = np.full((n+1, m+1), -INF)
    X = np.full((n+1, m+1), -INF)
    Y = np.full((n+1, m+1), -INF)

    tM = np.zeros((n+1, m+1), dtype=int)
    tX = np.zeros((n+1, m+1), dtype=int)
    tY = np.zeros((n+1, m+1), dtype=int)

    # ---------- INITIALIZATION ----------
    M[0,0] = 0

    for i in range(1, n+1):
        X[i,0] = gap_open + (i-1)*gap_extend
        tX[i,0] = 1

    for j in range(1, m+1):
        Y[0,j] = gap_open + (j-1)*gap_extend
        tY[0,j] = 2

    # ---------- FILL MATRICES ----------
    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            prev = [M[i-1,j-1], X[i-1,j-1], Y[i-1,j-1]]
            M[i,j] = max(prev) + score
            tM[i,j] = np.argmax(prev)

            open_gap = M[i-1,j] + gap_open + gap_extend
            extend_gap = X[i-1,j] + gap_extend

            if open_gap >= extend_gap:
                X[i,j] = open_gap
                tX[i,j] = 0
            else:
                X[i,j] = extend_gap
                tX[i,j] = 1

            open_gap = M[i,j-1] + gap_open + gap_extend
            extend_gap = Y[i,j-1] + gap_extend

            if open_gap >= extend_gap:
                Y[i,j] = open_gap
                tY[i,j] = 0
            else:
                Y[i,j] = extend_gap
                tY[i,j] = 2

    # ---------- FINAL SCORE ----------
    matrices = [M[n,m], X[n,m], Y[n,m]]
    matrix = np.argmax(matrices)
    score = matrices[matrix]

    # ---------- TRACEBACK ----------
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:

        if matrix == 0:
            prev = tM[i,j]
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1; j -= 1
            matrix = prev

        elif matrix == 1:
            prev = tX[i,j]
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
            matrix = prev

        else:
            prev = tY[i,j]
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
            matrix = prev

    return score, align1, align2


# ---------- MAIN PROGRAM ----------
file1 = input("Enter FASTA file 1: ")
file2 = input("Enter FASTA file 2: ")

seq1 = read_fasta(file1)
seq2 = read_fasta(file2)

match = int(input("Match score: "))
mismatch = int(input("Mismatch penalty: "))
gap_open = int(input("Gap opening penalty: "))
gap_extend = int(input("Gap extension penalty: "))

score, a1, a2 = affine_global_alignment(
    seq1, seq2, match, mismatch, gap_open, gap_extend
)

print("\n===== BEST GLOBAL ALIGNMENT =====")
print("Score:", score)
print(a1)
print(a2)
