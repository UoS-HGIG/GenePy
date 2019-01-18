# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:35:01 2018

@author: Enrico
"""
import numpy as np
import sys

gdi = ([line.rstrip().split("\t") for line in open('GDI_scores.txt',"r")])
gdi.pop(0)
gdi = np.array(gdi)
_f = dict(gdi[:,[0,2]])

data = np.loadtxt(sys.argv[1], dtype='str', delimiter='\t')
header = data[0,1:]
pzid = data[1:,0].T
genes = data[0,1:]
matrix = data[1:,1:].astype(float)

C=np.copy(matrix)
failed=[]
success=[]
for _i,_n in enumerate(genes):
    try:
        C[:,_i] = np.true_divide(C[:,_i], float(_f[_n]))
        success.append(_i)
    except:
        failed.append(_i)

N=np.vstack((genes[success],C[:,success]))
a = np.hstack((np.reshape(data[:,0],(404,1)),N))

np.savetxt(sys.argv[2],a, fmt="%s")
