#rlib/rrawarray.py! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#  Copyright © XYM
# CreateTime: 2017-01-13 15:34:46

#from array import array
from math import sqrt, erf, exp, pow
import math
import itertools
#from itertools import izip
from random import random
from rpython.rlib import rrandom
from rpython.rlib.rfloat import erfc
from rpython.rtyper.lltypesystem.rffi import r_ushort #r_uint, r_int_fast64_t
from rpython.rlib.rarithmetic import intmask, r_uint32, r_uint
from rpython.rlib import rfile
from rpython.rlib import rmmap
from rpython.rlib.listsort import TimSort
#from struct import pack
from time import time
import gc

# quicksort
def qsort(x, y, k = 10):
    n = len(x)
    for i in xrange(n):
        for j in xrange(i + 1, n):
            xi, xj = x[i], x[j]
            if y[xi: xi + k] > y[xj: xj + k]:
                x[i], x[j] = x[j], x[i]


# map function
def map(fuc, arr):
    #n = len(arr)
    #return [fuc(arr[elem]) for elem in xrange(n)]
    return [fuc(elem) for elem in arr]


# the izip function
def izip(seqs):
    iseqs = [iter(elem) for elem in seqs]
    while 1:
        try:
            yield [seq.next() for seq in iseqs]
        except:
            break

# the pack
def pack(dtype, val):
    t = dtype.lower()
    if t == 'h':
        n = 2
    elif t == 'i':
        n = 4
    elif t == 'l':
        n = 8
    else:
        n = 1

    string = [''] * n
    for i in xrange(n):
        string[i] = chr(val & 0xFF)
        val = val >> 8

    return ''.join(string)



# open a file
# readline from an open file
#def fopen(name, chk = 128):
def readline(f, chk = 128):

    #qry = args[0]
    #f = open(name, 'r')
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
            yield '%s\n'%line

        c = b[-1]

    #f.close()

# use rmmap to make realine code simple
def readline_mmap(f):
    handle = rmmap.mmap(f.fileno(), 0, access = rmmap.ACCESS_READ)
    for i in xrange(handle.size):
        line = handle.readline()
        if line != '':
            yield line
        else:
            break



class Seq:
    def __init__(self, id = '', seq = '', description = '', qual = ''):
        self.id = id
        self.seq = seq
        self.description = description
        self.qual = qual
    def update(self, id = '', seq = '', description = '', qual = ''):
        self.id = id
        self.seq = seq
        self.description = description
        self.qual = qual


# parse the fastq or fasta
#def parse(f, mode = 'fastq'):
def parse(f, dtype = 'fastq'):

    seq = Seq()
    output = []
    #for i in f:
    for i in readline_mmap(f):
        #print 'seq is', i
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
    n = len(s)
    N = 0
    output = 0
    #for i in s[::-1]:
    for c in xrange(n - 1, -1, -1):
        i = s[c]
        #print 'char', c, i, len(s), code[ord(i)]
        output += code[ord(i)] * int(pow(scale, N))
        #output += code[ord(i)] * pow(scale, N)
        N += 1

    return int(output)

s2n_next = lambda start, ksize, scale, code, char: start % int(pow(scale, (ksize - 1))) * scale + code[ord(char)]

# convert kmer to number from a seq
def seq2n(s, k = 15, scale = 4, code = code):
    n = len(s)
    if n <= k:
        yield k2n(s, scale, code)
    else:
        start = k2n(s[0: k], scale, code)
        yield start
        for i in s[k: n]:
            start = s2n_next(start, k, scale, code, i)
            yield start


# kmer count function
def kmc(seq, bucket, ksize, scale):
    for i in seq2n(seq, ksize, scale):
        val = intmask(bucket[i])
        val = val < 65535 and val + 1 or val
        bucket[i] = r_ushort(val)


# count the kmer from the paired fastq file
def fq2c(qry):
    ksize = max(2, int(qry[0]))
    names = qry[1: ]
    D = len(names)
    scale = 4
    #buck = [array('H', [0]) * pow(scale, k) for elem in xrange(D // 2)]
    #buck = [array('H', [0]) * scale ** k for elem in xrange(D // 2)]
    #buck = np.zeros((D//2, scale ** k), 'int32')
    N = int(pow(scale, ksize))
    buck = [[r_ushort(0) for elem0 in xrange(N)] for elem1 in xrange(D // 2)]

    print 'ksize, names, D, N, buck shape', ksize, names, D, N, len(buck), len(buck[0])

    flag = 0
    itr = 1
    init0 = init1 = time()
    #for name0, name1 in zip(names[0::2], names[1::2]):
    for i in xrange(0, len(names), 2):
        name0, name1 = names[i], names[i + 1]
        print 'qry0', name0, 'qry1', name1
        f0 = open(name0, 'r')
        f1 = open(name1, 'r')
        seqs0 = parse(f0, 'fastq')
        #seqs0 = SeqIO.parse(f0, 'fastq')
        seqs1 = parse(f1, 'fastq')
        #seqs1 = SeqIO.parse(f1, 'fastq')
        for seq0, seq1 in izip([seqs0, seqs1]):
            if seq0.qual.count('N') > 20 or seq1.qual.count('N') > 20:
                continue

            kmc(seq0.seq, buck[flag], ksize, scale)
            kmc(seq1.seq, buck[flag], ksize, scale)

            # timing
            if itr % 1000000 == 0:
                print 'iteration', itr, 'time', time() - init1, 'total time', time() - init0
                init1 = time()
            itr += 1

        f0.close()
        f1.close()

        flag += 1

    #print 'buck', [[intmask(b1) for b1 in b0] for b0 in buck]

    # convert the each reads to a freq hist
    flag = 0
    itr = 1
    init0 = init1 = time()
    #for name0, name1 in zip(names[0::2], names[1::2]):
    for i in xrange(0, len(names), 2):
        name0, name1 = names[i], names[i + 1]
        print 'qry0', name0, 'qry1', name1
        f0 = open(name0, 'r')
        f1 = open(name1, 'r')
        seqs0 = parse(f0, 'fastq')
        #seqs0 = SeqIO.parse(f0, 'fastq')
        seqs1 = parse(f1, 'fastq')
        #seqs1 = SeqIO.parse(f1, 'fastq')
        for seq0, seq1 in izip([seqs0, seqs1]):
            if seq0.qual.count('N') > 20 or seq1.qual.count('N') > 20:
                continue

            output = []
            fwd = [intmask(buck[flag][elem]) for elem in seq2n(seq0.seq, ksize, scale, code)]
            rev = [intmask(buck[flag][elem]) for elem in seq2n(seq1.seq, ksize, scale, code)]
            med = len(fwd) // 2
            step = 16
            start, end = max(0, med - step), max(0, med + step)
            # for test, only keep positive
            output = fwd[start: end] + rev[start: end]
            #output = fwd[start: end]
            #output.sort()
            TimSort(output).sort()

            if len(output) < step * 2:
                output.extend([0] * (step * 2 - len(output)))

            #print max(output), min(output), mean(output), std(output), output
            #A, B, C, D = khist(output)
            #print ' '.join(map(str, A + ['|', B, C, D]))
            #print ' '.join(map(str, output))
            print 'output', ' '.join([str(out) for out in output])


            # timing
            if itr % 1000000 == 0:
                print 'iteration', itr, 'time', time() - init1, 'total time', time() - init0
                init1 = time()
            itr += 1

        f0.close()
        f1.close()

        flag += 1

    return buck





# cdf of normal distribution
ncdf = lambda x : erfc(-x / 1.4142135623730951) / 2


# the sum, mean and std of an array/list
def sum(x):
    flag = 0
    #a = [[1] * 100 for elem in xrange(pow(10, 6))]
    for i in x:
        flag += i
        #a.append(i)
    return flag


mean = lambda x: 1. * sum(x) / len(x)

def std(x):
    if len(x) <= 1:
        return 0
    else:
        m = mean(x)
        var = sum([(elem - m) * (elem - m) for elem in x]) / (len(x) - 1.)
        return sqrt(var)


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


# the euclidean distance
def euclidean(x, y):
    d = 0
    for i in xrange(len(x)):
        z = (x[i] - y[i])
        d += z * z

    return sqrt(d)

# dist between x and y
# the dist can be normalized by the norm of x or y
def dist(x, y, norm = True):
    nx, ny, nz = [sqrt(elem) for elem in map(sum, [[xi * xi for xi in x], [yi * yi for yi in y], [(yi - xi) * (yi - xi) for xi, yi in zip(x, y)]])]
    return norm and nz / nx or nz / ny

# the pearson relationship
def pearson(x, y):
    N, M = len(x), len(y)
    assert N == M
    x_m, y_m = mean(x), mean(y)
    a = b = c = 0.
    for i in xrange(N):
        xi, yi = x[i] - x_m, y[i] - y_m
        a += xi * yi
        b += xi * xi
        c += yi * yi

    d = sqrt(b * c)
    if d > 0:
        return a / d
    else:
        return b == c and 1. or 0.



# the canopy algoithm
# for euc, t1 > t2
# for cor, t1 < t2
#def canopy(data, t1 = 2., t2 = 1.5, dist = pearson):
def canopy(data, t1 = .2, t2 = .6, dist = mannwhitneyu):
#def canopy(data, t1 = 0, t2 = 1e-3, dist = pearson):

    canopies = []
    _o = open('canopy.npy', 'w')
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
        #x = [intmask(val) for val in data[idx]]
        x = map(intmask, data[idx])

        can = [idx]
        #can = array('i', [idx])
        #del can[:]
        #can.append(idx)

        keep = []
        #keep = array('i')
        #del keep[:]

        if dist == mannwhitneyu:
            print 'mann'
            for elem in idxs:
            #while idxs:
                #elem = idxs.pop()
                #y = [intmask(val) for val in data[elem]]
                y = map(intmask, data[elem])
                u, p = dist(x, y)
                if p > t1:
                    can.append(elem)
                if p < t2:
                    keep.append(elem)
                #print 'size', len(idxs)

        elif dist == pearson:
            print 'pearson'
            for elem in idxs:
            #while idxs:
                #elem = idxs.pop()
                #y = [intmask(val) for val in data[elem]]
                y = map(intmask, data[elem])
                cor = dist(x, y)
                if cor > t1:
                    can.append(elem)
                if cor < t2:
                    keep.append(elem)
                #print 'size', len(idxs)


        else:
            for elem in idxs:
            #while idxs:
                #elem = idxs.pop()
                #y = [intmask(val) for val in data[elem]]
                y = map(intmask, data[elem])

                d = dist(x, y)
                if d < t1:
                    can.append(elem)
                if d > t2:
                    keep.append(elem)

        # use -1 as sep and save to disk
        can.append(-1)
        canopies.extend(can)
        if len(canopies) >= 1000000:
            string = ''.join([pack('i', elem) for elem in canopies])
            _o.write(string)
            canopies = []
        #canopies.write(string)
        #print 'can size', len(can), len(keep)
        print 'can size', len(can), len(keep), len(canopies)
        #del x, can, idxs
        idxs = keep
        #del keep
        gc.collect()
        print 'use time', time() - init1_time, 'total', time() - init0_time
        init1_time = time()

    if len(canopies) > 0:
        string = ''.join([pack('i', elem) for elem in canopies])
        _o.write(string)
        canopies = []

    #canopies.close()
    return canopies


#def run(x, y, n):
def run(n, qry):
    ksize = max(2, n)
    rg = rrandom.Random()
    a = []
    for i in xrange(n):
        #b = [r_ushort(int(rg.random() * pow(2, 15) - 1)) for elem0 in xrange(32)]
        b = [int(rg.random() * pow(2, 15) - 1) for elem0 in xrange(32)]
        #b.sort()
        TimSort(b).sort()
        #a.append(b)
        a.append([r_ushort(elem) for elem in b])
        #a.append(d)
    #a = [TimSort([r_ushort(int(rg.random() * pow(2, 15) - 1)) for elem0 in xrange(32)]).sort() for elem1 in xrange(n)]
    #print 'short add', a[0][0] + a[0][0]
    u = p = 0
    Atmp = range(6)
    dna = 'atgcgc'
    qsort(Atmp, dna)
    print 'qsort', Atmp, [dna[elem] for elem in Atmp]
    for i in xrange(1):
        #x = [0] * 32
        #y = [0] * 32
        #for j in xrange(32):
        #    x[j] = rg.random()
        #    y[j] = rg.random()
        #x = [rg.random() for elem in xrange(32)]
        #y = [rg.random() for elem in xrange(32)]
        x = [intmask(elem) for elem in [r_ushort(0)] * 32]
        #x = map(intmask, [0] * 32)
        y = [intmask(elem) for elem in xrange(32)]
        #rgs = range(32)
        #y = map(intmask, rgs)
        #x = y = [0] * 32
        u, p = mannwhitneyu(x, y, True)
        p = pearson(x, y)

    print 'p value is', u, p

    print k2n('ta' * 12), intmask(int('123'))
    test_seq = 'tgatcgctgtagctgatgctcatgctatgctatcgtagtcgtgctagctagcatcgatcgatcgctagaaacagctgcgtatctatctatatatatattaggagaatgtgagaga'
    test_n = seq2n(test_seq)
    #for i in test_n:
    #    print 'seq2n test', i
    canopy(a)

    print [r_uint(elem) for elem in test_n]
    buck = [[r_ushort(0) for elem0 in xrange(pow(scale, ksize))] for elem1 in xrange(len(qry))]
    f0 = open(qry, 'r')
    f1 = open(qry, 'r')
    seqs0 = parse(f0)
    seqs1 = parse(f1)
    for seq0, seq1 in izip([seqs0, seqs1]):
        #print seq0.seq, seq1.seq
        print 'test seq2n'
        #print [nb for nb in seq2n(seq0.seq)]
        kmc(seq0.seq, buck[0], ksize, scale)
        kmc(seq1.seq, buck[0], ksize, scale)


    f0.close()
    f1.close()

    print [intmask(elem) for elem in buck[0][:15]]
    '''
    while 1:
        try:
            seq0 = seqs0.next()
            seq1 = seqs1.next()
            print seq0.seq, seq1.seq
        except:
            break

    return 1
    '''

def entry_point(argv):
    try:
        K = int(argv[1])
    except:
        K = 10

    qry = argv[2]

    fac = lambda x: sum(range(x))
    print fac(K)

    #x = range(32)
    #y = range(32, 64)
    #run(x, y, K)
    #run(K, qry)
    d2 = [[r_ushort(1)] * K, [r_ushort(4), r_ushort(4)] * K]
    print 'D2 shape', len(d2), len(d2[1])
    del d2;
    gc.collect()
    buck = fq2c(argv[1:])
    return 0

def target(*args):
    return entry_point, None

if __name__ == "__main__":
    import sys
    entry_point(sys.argv)

