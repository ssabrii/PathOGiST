# https://louisabraham.github.io/notebooks/suffix_arrays.html

from itertools import zip_longest, islice

def inverse_array(l):
    n = len(l)
    ans = [0] * n
    for i in range(n):
        ans[l[i]] = i
    return ans

def to_int_keys_best(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    seen = set()
    ls = []
    for e in l:
        if not e in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def suffix_array_best(s):
    """
    suffix array of s
    O(n * log(n)^2)
    """
    n = len(s)
    k = 1
    line = to_int_keys_best(s)
    while max(line) < n - 1:
        line = to_int_keys_best(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)])
        k <<= 1
    return inverse_array(line), line

def lcp(s, sa, inv_sa):
    """
    Use Kasai algorithm to build LCP array
    http://www.mi.fu-berlin.de/wiki/pub/ABI/RnaSeqP4/suffix-array.pdf
    """
    n = len(inv_sa)
    lcp = [0] * n
    l = 0
    for i in range(n):
        if inv_sa[i] > 0:
            k = sa[inv_sa[i] - 1]
            while i + l < n and k + l < n and s[i + l] == s[k + l]:
                l += 1
            lcp[inv_sa[i]] = l
            if l > 0:
                l -= 1

    return lcp
