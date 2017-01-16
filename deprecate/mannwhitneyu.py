#rlib/rrawarray.py! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#  Copyright © XYM
# CreateTime: 2017-01-13 15:34:46


from math import sqrt, erf, exp, pow
#from rpython.flowspace.operation import inplace_pow
import math
import itertools
from itertools import izip
from random import random
from rpython.rlib import rrandom
from rpython.rlib.rfloat import erfc
from array import array
from rpython.rtyper.lltypesystem.rffi import r_ushort #r_uint, r_int_fast64_t
from rpython.rlib.rarithmetic import intmask, r_uint32, r_uint
from time import time
import gc

# open a file
def fopen(name, chk = 128):
    #qry = args[0]
    f = open(name, 'r')
    c = ''
    #b = []
    while 1:
        a = f.read(1024 * 1024 * chk)
        if not a:
            break
        #print a

        b = a.split('\n')
        #yield b
        b[0] = c + b[0]
        for line in b[:-1]:
            yield line

        c = b[-1]

    f.close()


class Seq:
    def __init__(self, id = None, seq = None, description = None, qual = None ):
        self.id = id
        self.seq = seq
        self.description = description
        self.qual = qual
    def update(self, id = None, seq = None, description = None, qual = None ):
        self.id = id
        self.seq = seq
        self.description = description
        self.qual = qual


# parse the fastq or fasta
def parse(f, mode = 'fastq'):
    seq = Seq()
    output = []
    for i in f:
        print 'seq is', i
        if len(output) == 4:
            #seq.id, seq.seq, seq.description, seq.qual = output
            #seq.update(*output)
            a, b, c, d = output
            seq.update(a, b, c, d)
            yield seq
            output = [i[: -1]]
        else:
            output.append(i[: -1])

    if len(output) == 4:
        #seq.id, seq.seq, seq.description, seq.qual = output
        a, b, c, d = output
        seq.update(a, b, c, d)
        yield seq


code = [0] * 256
flag = 0
for i in 'ATGC':
    code[ord(i.lower())] = code[ord(i)] = flag
    flag += 1

scale = max(code) + 1

# fast convert kmer to number
def k2n(s, scale = 4, code = code):
    if scale == -1:
        #scale = code.max() + 1
        scale = max(code) + 1
    N = 0
    output = 0
    #for i in s[::-1]:
    for c in xrange(len(s) - 1, -1, -1):
        i = s[c]
        #print 'char', c, i, len(s), code[ord(i)]
        output += code[ord(i)] * int(pow(scale, N))
        #output += code[ord(i)] * pow(scale, N)
        N += 1

    return int(output)

s2n_next = lambda start, ksize, scale, code, char: start % int(pow(scale, (ksize - 1))) * scale + code[ord(char)]

# convert kmer to number from a seq
def seq2n(s, k = 15, scale = 4, code = code):
    if len(s) <= k:
        yield k2n(s, scale, code)
    else:
        start = k2n(s[: k], scale, code)
        yield start
        for i in s[k: ]:
            start = s2n_next(start, k, scale, code, i)
            yield start


# convert kmer to number from a seq
def seq2n(s, k = 15, scale = 4, code = code):
    if len(s) <= k:
        yield k2n(s, scale, code)
    else:
        start = k2n(s[: k], scale, code)
        yield start
        for i in s[k: ]:
            start = s2n_next(start, k, scale, code, i)
            yield start



# cdf of normal distribution
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
def mannwhitneyu(x, y, use_continuity = True):

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

    if z >= 0:
        p = 2 - 2 * ncdf(z)
    else:
        p = 2 * ncdf(z)

    return u, p


# the canopy algoithm
# for euc, t1 > t2
# for cor, t1 < t2
#def canopy(data, t1 = 2., t2 = 1.5, dist = pearson):
def canopy(data, t1 = 0, t2 = 1e-3, dist = mannwhitneyu):
    canopies = []
    #canopies = open('canopies.npy', 'w')
    idxs = range(len(data))
    #idxs =  array('i')
    #for i in xrange(len(data)):
    #    idxs.append(i)

    # shuffle the index
    #shuffle(idxs)

    # the can and keep array
    init0_time = time()
    init1_time = time()
    flag = 0
    while idxs:
        idx = idxs.pop()
        x = [intmask(val) for val in data[idx]]

        can = [idx]
        #can = array('i', [idx])
        #del can[:]
        #can.append(idx)

        keep = []
        #keep = array('i')
        #del keep[:]

        if dist == mannwhitneyu:
            print 'pearson'
            for elem in idxs:
            #while idxs:
            #    elem = idxs.pop()
                y = [intmask(val) for val in data[elem]]
                u, p = dist(x, y)
                if p > t1:
                    can.append(elem)
                if p < t2:
                    keep.append(elem)
                #print 'size', len(idxs)

        #else:
        #    for elem in idxs:
        #    #while idxs:
        #    #    elem = idxs.pop()
        #        y = [intmask(val) for val in data[elem]]
        #        p = dist(x, y)
        #        if p < t1:
        #            can.append(elem)
        #        if p > t2:
        #            keep.append(elem)

        #canopies.append(can)
        # use -1 as sep and save to disk
        #can.append(-1)
        #can.tofile(canopies)
        canopies.extend(can)
        print 'can size', len(can), len(keep), len(canopies)
        #del x, can, idxs
        idxs = keep
        #del keep
        gc.collect()
        print 'use time', time() - init1_time, 'total', time() - init0_time
        init1_time = time()

    #canopies.close()
    return canopies



#def run(x, y, n):
def run(n, qry):
    rg = rrandom.Random()
    a = [[r_ushort(int(rg.random() * pow(2, 15) - 1)) for elem0 in xrange(32)] for elem1 in xrange(n)]
    #print 'short add', a[0][0] + a[0][0]
    u = p = 0
    for i in xrange(1):
        #x = [0] * 32
        #y = [0] * 32
        #for j in xrange(32):
        #    x[j] = rg.random()
        #    y[j] = rg.random()
        #x = [rg.random() for elem in xrange(32)]
        #y = [rg.random() for elem in xrange(32)]
        x = [intmask(elem) for elem in [r_ushort(0)] * 32]
        y = [intmask(elem) for elem in xrange(32)]
        #x = y = [0] * 32
        u, p = mannwhitneyu(x, y, True)
    print 'p value is', u, p

    print k2n('ta' * 12), intmask(int('123'))
    test_seq = 'tgatcgctgtagctgatgctcatgctatgctatcgtagtcgtgctagctagcatcgatcgatcgctagaaacagctgcgtatctatctatatatatattaggagaatgtgagaga'
    test_n = seq2n(test_seq)
    canopy(a)

    print [r_uint(elem) for elem in test_n]
    f0 = fopen(qry, 64)
    f1 = fopen(qry, 64)
    seqs0 = parse(f0)
    seqs1 = parse(f1)
    while 1:
        try:
            seq0 = seqs0.next()
            seq1 = seqs1.next()
            print seq0.seq, seq1.seq
        except:
            break

    return 1

def entry_point(argv):
    try:
        K = int(argv[1])
    except:
        K = 10

    qry = argv[2]
    #x = range(32)
    #y = range(32, 64)
    #run(x, y, K)
    run(K, qry)
    return 0

def target(*args):
    return entry_point, None

if __name__ == "__main__":
    import sys
    entry_point(sys.argv)
