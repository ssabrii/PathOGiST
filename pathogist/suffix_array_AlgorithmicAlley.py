# https://github.com/benfulton/Algorithmic-Alley/blob/master/AlgorithmicAlley/SuffixArrays/sa.py

from collections import defaultdict, Counter

def get_suffix_array(s):
    return sorted(range(len(s)), key=lambda i: s[i:])

def sort_bucket(s, bucket, order=1):
    d = defaultdict(list)
    for i in bucket:
        key = tuple(s[i:i+order])
        d[key].append(i)
    result = []
    for k, v in sorted(d.items()):
        if len(v) > 1:
            result += sort_bucket(s, v, order*2)
        else:
            result.append(v[0])
    return result

def suffix_array_ManberMyers(s):
    return sort_bucket(s, (i for i in range(len(s))))
