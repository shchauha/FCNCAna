import uuid
import os
import json
import pprint
import ROOT as r
import numpy as np

id1 = 11
id2 = 13
pt1 = 33.
pt2 = 21.
ht = 350.
eta1 = 1.1
eta2 = 2.1
nvtx = 22

r.gROOT.ProcessLine(".L common_utils.h")

for year in [2016, 2017, 2018]:
    print year
    print r.leptonScaleFactor(year, id1, pt1, eta1, ht)
    print r.leptonScaleFactorError(year, id1, pt1, eta1, ht)
    print r.flipRate(year, pt1, eta1)
    print r.flipRateError(year, pt1, eta1)
    print r.fakeRate(year, id1, pt1, eta1, ht, 1) 
    print r.fakeRateError(year, id1, pt1, eta1, ht, 1) 
    print r.alternativeFakeRate(year, id1, pt1, eta1, ht, 1) 
    print r.getTruePUw(year, nvtx)
    print r.getTruePUw(year, nvtx, 1)
    print r.getTruePUw(year, nvtx, -1)
    print r.triggerScaleFactor(year, id1, id2, pt1, pt2, eta1, eta2, ht, 1)
    print r.getLumi(year)

