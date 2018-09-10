from range_minimum_query import RangeMinimum
from itertools import combinations, chain
import suffix_array
import suffix_array_AlgorithmicAlley
import numpy as np
import time

def aligned(s_i, L, m):
    """
    Let l = si mod m, i.e., the starting position of the suffix si within a profile.
    Then this subroutine returns l/L if l is multiple of L, and −1 otherwise.
    """
    l = s_i % m
    return int(l / L) if l % L == 0 else -1

def HD(p_i, p_j, s, m):
    """
    Given two profiles pi and pj which share a substring of length L,
    starting at index lL, this subroutine computes
    the minimum of k and the Hamming distance between pi and pj.
    """
    # dist = 0
    # for t in range(0, m):
    #     dist += 1 if s[p_i * m + t] != s[p_j * m + t] else 0

    # return dist
    return np.count_nonzero(s[p_i * m: p_i * m + m] != s[p_j * m: p_j * m + m])

def fast_ham_distance(P, k):
    """
    An average‐case linear‐time algorithm to compute
    pairwise Hamming distances among a set of taxa,
    under a given Hamming distance threshold.

    Input: A set P of d profiles of length m each; an integer threshold 0 < k < m.
    Output: The set X of distinct pairs of profiles that are at Hamming distance at most k
    """
    elapsed_1 = 0
    elapsed_2 = 0
    elapsed_3 = 0

    # Initialization:
    d = len(P)
    m = len(P[0])
    s = np.concatenate(P)
    n = m * d
    L = int(m / (k+1))
    s_t1 = time.time()
    S, inv_S = suffix_array.suffix_array_best(s)
    # S = suffix_array_AlgorithmicAlley.suffix_array_ManberMyers(s)
    # inv_S = suffix_array.inverse_array(S)
    e_t1 = time.time()
    LCP = suffix_array.lcp(s, S, inv_S)
    # RMQ_lcp = RangeMinimum(LCP)
    HT = dict()
    elapsed_1 += (e_t1 - s_t1)

    # Candidate pairs enumeration:
    X = np.ones((d, d)) * k
    np.fill_diagonal(X, 0)
    l_p = -1
    C = [set() for i in range(0, int(m / L) + 1)]
    for i in range(1, n):
        l = LCP[i]
        if l >= L:
            p_i = int(S[i] / m)
            x = aligned(S[i], L, m)
            if x != -1:
                C[x].add(p_i)
            if l_p == -1:
                p_i_1 = int(S[i-1] / m)
                x = aligned(S[i-1], L, m)
                if x != 1:
                    C[x].add(p_i_1)
            l_p = l
        elif l_p != -1:
            # Pairs enumeration:
            s_t3 = time.time()
            for C_t in C:
                for p, q in combinations(C_t, 2):
                    if (p, q) not in HT:
                        HT[p, q] = 1
                        s_t2 = time.time()
                        delta = HD(p, q, s, m)
                        e_t2 = time.time()
                        elapsed_2 += (e_t2 - s_t2)
                        if delta <= k:
                            X[p][q] = delta
                            X[q][p] = delta
                C_t.clear()
            l_p = -1
            e_t3 = time.time()
            elapsed_3 += (e_t3 - s_t3)
            # C = [set() for i in range(0, int(m / L) + 1)]

    print('Suffix array construction:', elapsed_1)
    print('HD:', elapsed_2)
    print('Pairs enumeration:', elapsed_3)
    return X
