#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#  Copyright © XYM
# CreateTime: 2017-01-24 20:58:41

import numpy as np
from sklearn.neighbors import KDTree, BallTree
from time import time
import sys

m, n = map(int, map(eval, sys.argv[1: ]))

x = np.random.randn(m // 2, n)
x = np.concatenate([x, x], 0)

t0 = time()
clf = KDTree(x)
print 'construct time', time() - t0

t0 = time()
print 'single point', clf.query(x[: m]), time() - t0

t0 = time()
print 'radius point', clf.query_radius(x[: 100000] + 1e-5, 1)
print 'query time', time() - t0

