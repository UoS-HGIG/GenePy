import numpy as np
ids = np.loadtxt('IDs', dtype='string')
geni = np.loadtxt('geni',dtype='string')

M = np.zeros((ids.shape[0], geni.shape[0]))

dic_geni={}
dic_pz={}

for i,n in enumerate(geni):
    dic_geni[n]=i
for i,n in enumerate(ids):
    dic_pz[n]=i
for n in ids:
	for x in open(n+'_per_gene_mean_coverage.out','r'):
		x = x.rstrip().split(' ')
		if x[0] != 'MIR4253':
			i = dic_pz[n]
			j = dic_geni[x[0]]
			M[i,j] = float(x[1])


ids = ids.reshape((ids.shape[0],1))
geni = geni.reshape((1, geni.shape[0]))

geni = np.insert(geni, 0, 'Samid')

M = np.hstack((ids, M))
M = np.vstack((geni, M))

np.savetxt('full_mean_cov_matrix.txt', M, delimiter='\t', fmt='%s')
