#rlib/rrawarray.py! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#  Copyright © XYM
# CreateTime: 2017-01-13 15:34:46


#from math import sqrt, erf, erfc, pow
from math import sqrt, pow, erf, exp
import math
import itertools
from random import random
from rpython.rlib import rrandom
from rpython.rlib.rfloat import erfc
from array import array


ncdf = lambda x : erfc(-x / 1.4142135623730951) / 2


def sum(x):
    flag = 0
    #a = [[1] * 100 for elem in xrange(pow(10, 6))]
    for i in x:
        flag += i
        #a.append(i)
    return flag

# rank the sorted data
def rankdata(val):
    n = len(val)
    ts = []
    rank = range(n)
    start = end = 0
    for i in xrange(1, n):
        if val[start] < val[i]:
            if start < end:
                #ri = float(start + end) / (i - start)
                ri = (start + end) / 2.
                ts.append(i - start)
                #print 'hello', ri, i, start
                for j in xrange(start, i):
                    rank[j] = ri

            start = i
        else:
            end = i

    if start < end:
        ri = float(start + end) / 2.
        ts.append(n - start)
        #print 'hello', ri
        for j in xrange(start, n):
            rank[j] = ri

    # correct start from 0
    #for i in xrange(n):
    #    rank[i] += 1

    return rank, ts

# tiecorrect
#tiecorrect = lambda ts, n: ts and 1 - sum([t ** 3. - t for t in ts]) / (n ** 3 - n) or 1
tiecorrect = lambda ts, n: ts and 1 - sum([pow(t, 3) - t for t in ts]) / (pow(n, 3) - n) or 1

# Mann–Whitney U test
# return the z score of U
#def mannwhitneyu(x, y, use_continuity = True, alternative = 'two-sided'):
def mannwhitneyu(x, y, use_continuity):

    #n0, n1 = map(len, [x, y])
    n0, n1 = len(x), len(y)
    n = n0 + n1
    val = [0] * n
    lab = [0] * n
    p0 = p1 = p2 = 0
    while p0 < n0 and p1 < n1:
        if x[p0] < y[p1]:
            val[p2] = x[p0]
            lab[p2] = 0
            p0 += 1
            p2 += 1

        elif x[p0] > y[p1]:
            val[p2] = y[p1]
            lab[p2] = 1
            p1 += 1
            p2 += 1

        else:
            val[p2] = x[p0]
            lab[p2] = 0
            p0 += 1
            p2 += 1

            val[p2] = y[p1]
            lab[p2] = 1
            p1 += 1
            p2 += 1

    while p0 < n0:
        val[p2] = x[p0]
        lab[p2] = 0
        p0 += 1
        p2 += 1

    while p1 < n1:
        val[p2] = y[p1]
        lab[p2] = 1
        p1 += 1
        p2 += 1

    # correct the tie rank
    # tied rank
    rank, ts = rankdata(val)
    #rank0 = sum([a for a, b in izip(rank, lab) if b == 0])
    # add n0 to correct the start from 0
    #Rank0 = [a for a, b in izip(rank, lab) if b == 0]
    Rank0 = [rank[i] for i in xrange(n) if lab[i] == 0]

    rank0 = sum(Rank0) + n0


    u0 = n0 * n1 + n0 * (n0 + 1) / 2. - rank0
    #u1 = n1 + sum([a for a, b in zip(rank, lab) if b == 1]) - n1 * (n1 + 1) / 2
    u1 = n0 * n1 - u0

    u = min(u0, u1)
    mu = n0 * n1 / 2.
    var = n0 * n1 * (n + 1) * tiecorrect(ts, n) / 12.
    #sd = var ** .5
    sd = pow(var, .5)
    z = use_continuity and abs(u - mu) - .5 / sd or abs(u - mu) / sd

    #print sqrt(z)
    #print erf(z)
    return z

#def run(x, y, n):
def run(n):
    rg = rrandom.Random()
    a = [[1] * 100 for elem in xrange(n)]

    for i in xrange(n):
        #x = [0] * 32
        #y = [0] * 32
        #for j in xrange(32):
        #    x[j] = rg.random()
        #    y[j] = rg.random()
        #x = [rg.random() for elem in xrange(32)]
        x = [0] * 32
        y = [0] * 32
        #y = [rg.random() for elem in xrange(32)]
        z = mannwhitneyu(x, x, True)
    #print z
    return 1

def entry_point(argv):
    try:
        K = int(argv[1])
    except:
        K = 10

    #x = range(32)
    #y = range(32, 64)
    #run(x, y, K)
    run(K)
    return 0

def target(*args):
    return entry_point, None

if __name__ == "__main__":
    entry_point(sys.argv)
