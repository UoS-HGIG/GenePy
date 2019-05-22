# -*- coding: utf-8 -*-
"""
Created on Wed Sept 27 12:02:04 2017

@author: Enrico
"""
import sys, requests
import numpy as np
import collections as cl
from scipy.stats import ranksums as rks



server = "http://rest.ensembl.org/variation/human/"

#%% INPUT files

def import_data():
    IBD=([line.rstrip().split('\t') for line in open(sys.argv[1],'r')])
    LEV=([line.rstrip().split('\t') for line in open(sys.argv[2],'r')])
    return IBD, LEV

cases, levin = import_data()
gene=sys.argv[3]

#def bubu():
#    cases=([line.rstrip().split('\t') for line in open('test_IBD_combined.txt','r')])
#    levin=([line.rstrip().split('\t') for line in open('test_lev_combined.txt','r')])
#    gene='NOD2'
#    return cases, levin, gene
#cases, levin, gene = bubu()
cases_header=np.array(cases.pop(0))
levin_header=np.array(levin.pop(0))

#%% consolidate databases
### IBD CD CASES  ###
cases=np.array(cases)

scores=cases[:,6:22]
scores_names=cases_header[6:22]

rsids = cases[:,23]

proxies = {'http' : "socks5://localhost:6789"} #on if on iridis
#proxies = {}

freqs = np.zeros((rsids.shape[0], 2))

def get_alts_freq(r, a, rs):
    if rs != 'nan': # is variant associated to rsid?
        ref = np.nan
        alt = np.nan
        ext = rs+"?pops=1"
#        v = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        v = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, proxies=proxies) #if on iridis
        v = v.json()
        fr = 0
        fa = 0
        if '1000Genomes' in v['evidence']:
            for x in v['populations']:
                if x['population']=='1000GENOMES:EUR' and x['allele']==r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population']=='1000GENOMES:EUR' and x['allele']==a and fa == 0:
                    alt = x['frequency']
                    fa = 1
                
        elif '1000Genomes' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == '1000GENOMES:phase_3:EUR' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == '1000GENOMES:phase_3:EUR' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1
                    
        elif 'HapMap' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == 'CSHL-HAPMAP:HapMap-CEU' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == 'CSHL-HAPMAP:HapMap-CEU' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1         
                    
        elif 'ESP' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == 'ESP6500:European_American' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == 'ESP6500:European_American' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1
    
        if fa == 0 and fr == 0:
            ref = 0.9999
            alt = 0.0001
            
        elif fa == 0 and fr == 1:
            alt = 1-fr
            
        elif fa == 1 and fr == 0:
            ref = 1-fa
    else:# add novel freq if no rsid entry
        ref =0.9999
        alt = 0.0001

    return ref, alt

for i, n in enumerate(rsids):
    freqs[i,0], freqs[i,1] = get_alts_freq(cases[i,2], cases[i,3], n)
    
freqs[freqs==0]=0.0001 # set 0 as minimum otherwise log fails


cases=cases[:,24:].astype(float)/2.00
cases_header=cases_header[24:]        


####  LEVIN CONTROLS ###
levin=np.array(levin)
levin=levin[:,24:].astype(float)/2.00
levin_header=levin_header[24:]

final_sample_list=np.hstack((cases_header,levin_header))

#%% Clalculate the distributions and test

def score_db(ibd,lev,score,freq,sname):
    S=np.copy(score)
    db1=np.copy(ibd)
    db2=np.copy(lev)

    out1=[]
    for i in range(db1.shape[0]):
        if S[i]!='nan':
            deleter=float(S[i])
            db1[i][db1[i]==0]=deleter*-np.log10(float(freq[i,0])*float(freq[i,0]))
            db1[i][db1[i]==0.5]=deleter*-np.log10(float(freq[i,0])*float(freq[i,1]))
            db1[i][db1[i]==1]=deleter*-np.log10(float(freq[i,1])*float(freq[i,1]))
            out1.append(db1[i])
    out2=[]
    for i in range(db2.shape[0]):
        if S[i]!='nan':
            deleter=float(S[i])
            db2[i][db2[i]==0]=deleter*-np.log10(float(freq[i,0])*float(freq[i,0]))
            db2[i][db2[i]==0.5]=deleter*-np.log10(float(freq[i,0])*float(freq[i,1]))
            db2[i][db2[i]==1]=deleter*-np.log10(float(freq[i,1])*float(freq[i,1]))
            out2.append(db2[i])
    
    out1=np.array(out1)
    out2=np.array(out2)

    
    out1 = np.sum(out1,axis=0)
    out2 = np.sum(out2,axis=0)
    
    tmp = np.hstack((out1,out2))
    gg  =np.array([gene]*len(final_sample_list))
    U = np.vstack((final_sample_list,tmp,gg)).T
    
    return rks(out1,out2), U


print '#Gene\tScore\tIBD_vs_Levin'

for i in range(scores.shape[1]):
    if cl.Counter(scores[:,i])['nan'] < (scores.shape[0]-1): #compute metascores if at least 1 variant
        pvs, U =score_db(cases,levin,scores[:,i],freqs, scores_names[i])
        print '%s\t%s\t%.5f' %(gene,scores_names[i],pvs[1])
        np.savetxt('./'+scores_names[i]+'/'+gene+'_'+scores_names[i]+'_matrix',U, fmt='%s', delimiter='\t')
    else:
        print '%s\t%s\t%s\n' %(gene,scores_names[i],'-')
