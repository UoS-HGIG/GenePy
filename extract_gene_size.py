# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 11:16:52 2018

@author: Enrico
"""

import numpy as np

# PART 1
#%%
def get_genesize():
   bed = np.loadtxt('INTERSECTED_TRGT_nodash.bed', dtype='O', delimiter='\t')
   c = bed[:,2].astype(int) - bed[:,1].astype(int)

   geni = np.unique(bed[:,-1])

   size=[]
   for _g in geni:
       _s = np.sum(c[np.where(bed[:,-1]==_g)])
       size.append([_g,_s])
   size = np.array(size)
   return size
size = get_genesize()
np.savetxt('hg19_sizes.txt', size, fmt='%s', delimiter='\t')


#PART 2
 size = np.loadtxt('hg19_sizes.txt',dtype='str', delimiter='\t')
 def norm_by_size(X, genes):
     C = np.copy(X)
     _f = dict(size)
     for _i,_n in enumerate(genes):
         C[:,_i] = np.true_divide(C[:,_i], int(_f[_n]))*np.mean(size[:,1].astype(int))
     return C

 data = np.loadtxt('REVEL_GenePy_matrix', dtype='str', delimiter='\t')
 header = data[0,1:]
 pzid = data[1:,0].T
 genes = data[0,1:]
 matrix = data[1:,1:].astype(float)

 W = norm_by_size(matrix, genes)

# F = np.vstack((header,W))
# g = data[:,0,np.newaxis]
# F = np.hstack((g, F))

# np.savetxt('REVEL_GenePy_norm_matrix', F, fmt='%s', delimiter='\t')
