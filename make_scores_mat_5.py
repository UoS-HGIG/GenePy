# -*- coding: utf-8 -*-
"""
Created on Wed Sept 27 12:02:04 2017

@author: Enrico

Change log:
This version is based on ALL frequencies and not EUR.
gnomeAD interrogation has been switched, first genome then exome

14.11.18 Removed the request to ensmbl, using only bi-allelic loci to simplify
script. Much faster. Replaced macap (only for rare SNVs) with eigen. gnomAD exome
is now the pricipal dataset for frequencies. Novel frequencies are set to one 
individual out 125,748 in gnomAD_exome.
"""
import sys, re
import numpy as np



# IMPORT data passed through the first argument
# the gene variable is equal to the gene name input provided
# the header is removed from the data variable and turned into an array to make it easier to parse through
def import_data():
    try:
        X =([line.rstrip() for line in open(sys.argv[1],'r')])
#       example: X = ([line.rstrip() for line in open('NOD2.meta','r')])
    except:
        print('Not imported')
    return X
data = import_data()
gene=sys.argv[2]
# example: gene = 'NOD2'
#%%
header=np.array(data.pop(0).split('\t'))

# since the header has been removed from the data variable, there's only data in that variable
# now the data will be reformated by replacing 0/0 or ./. (which means missing) with 0 (homozygous reference)
# 0/any number is replaced as 1 (heterozygous); 
# any number >0/any number >0 will be 2 (homozygous alternative){should include only 1/1 as per previous filter on biallelic variants. Left for a future development on multiallelic scoring

for _i,_n in enumerate(data):
    data[_i] = re.sub('0/0[\S]+','0', data[_i])
    data[_i] = re.sub('0/[123456789][\S]+','1', data[_i])
    data[_i] = re.sub('[123456789]/[123456789][\S]+','2', data[_i])
    data[_i] = re.sub('\./\.[\S]*','0', data[_i])
    data[_i] = data[_i].split('\t')

# this portion is to clean up the data and header information. First turn the data from a list type to array to make it easier to deal with. 
# then delete portions from the data and header file: positions 20 and 25 (19:CADD_RawScore, 20:CADD13_PHRED, 21: Eigen, 22: REVEL, 
# 23:GWAVA_region_score, 24:GWAVA_tss_score, 25:GWAVA_unmatched_score, 26:dann, 27:Ctrl_1)
data = np.array(data)
data = np.delete(data,[20,25],axis=1)
header = np.delete(header,[20,25])

# after this portion of code the header and data: (19:CADD_RawScore, 20: Eigen, ... 23: GWAVA_tss_score, 24: dann)

scores = data[:,19:25]
scores[scores=='.'] = np.nan
scores = scores.astype('float')
scores_names = header[19:25]
# reformat CADD range
scores[:,0] = (scores[:,0] - (-7.535037))/(35.788538-(-7.535037))

#%% gnomADexome FREQUENCIES (assume biallelic)
# this is the list of all the gnomAD_exome_ALL data which is in position 10 of the data list
# also turns all the non-existant data points into nans and conver to float variables
known_fa_ALL = data[:,10]
known_fa_ALL[known_fa_ALL=='.']=np.nan
known_fa_ALL = known_fa_ALL.astype('float')

# freqs is a 2D array with two positions, which are initialized as [nan nan], for each variants. 
# the second value in that array is replaced by the known_fa_ALL list created before for each data point
# for instance if the known_fa_ALL was [2,3,4,5,6,7,8,9] then the freqs = [[0,2], [0,3], [0,4],..., [0,9]]
# if any of the values in known_fa_ALL were 'nan' or '.' or zero then it is replaced with 7.95e-6 (1 out 125,748 indiv in gnomADexome)
# Any values that are 1 are replaced by 1 - 7.95e-6
# this value is the the frequency of the mutation and the first position is the frequency of the reference allele
# at the first position (the nan) of the 2D array [nan somevalue], the nan is replace with the
# 1 - somevalue (the value at the second position of the array)
freqs = np.zeros((data.shape[0], 2))
freqs[:] = np.nan
freqs[:,1]= known_fa_ALL
freqs[:,1][np.isnan(freqs[:,1])]=7.95e-6  # 1 out 125,748 indiv in gnomADexome
freqs[:,1][freqs[:,1]==0]=7.95e-6
freqs[:,1][freqs[:,1]==1]=1-7.95e-6
freqs[:,0]= 1-freqs[:,1]


#%%   
# genotypes array which is a subset of the data object from the positions from 25 (first sample in the vcf) up til the end. Also takes each value and divides by 2
samples=data[:,25:].astype(float)/2.00
samples_header=header[25:]        


#%% Clalculate GenePy score
# loop though all the score y values (variants) / number of columns e.g. Each loop will consist of CADD_RawScore, Eigen, ..., GWAVA_tss_score, dann
# It will compute a score using the score_db method as long as there is at least one variant
# It will then create a document for each score: and then store them in the folders that was created in the GenePy_Pipline script
def score_db(samples,score,freq):
	# first make copies of the score and samples into S and db1, respectively
    S=np.copy(score)
    db1=np.copy(samples)

    out1=[]
    for i in range(db1.shape[0]):
        if ~np.isnan(S[i]): #if deleteriousness score is available
            deleter=float(S[i])# store the score value into the into the deleter variable 
            db1[i][db1[i]==0.5]=deleter*-np.log10(float(freq[i,0])*float(freq[i,1]))#compute GenePy score for heterozygous variants
            db1[i][db1[i]==1]=deleter*-np.log10(float(freq[i,1])*float(freq[i,1]))
            out1.append(db1[i]) #compute GenePy score for homozygous variants
            
    out1=np.array(out1) # then these values will be stored into the out1 array.

    out1 = np.sum(out1,axis=0) # the out1 is then condensed by suming the columns for each sample.
	
	# formmating the data into columns of the Sample Heading, Value, Gene Name. 
    gg = np.array([gene]*len(samples_header))
    U = np.vstack((samples_header,out1,gg)).T
    
    return U
#%%
for i in range(scores.shape[1]):
    if (np.isnan(scores[:,i]).sum()) < (scores.shape[0]-1): #compute metascores if at least 1 variant
        U = score_db(samples,scores[:,i],freqs)
        np.savetxt('./'+scores_names[i]+'/'+gene+'_'+scores_names[i]+'_matrix',U, fmt='%s', delimiter='\t')
