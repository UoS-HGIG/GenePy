# -*- coding: utf-8 -*-
"""
Created on Thu Feb 09 15:18:45 2017

@author: Enrico
"""
## FIRST FILE IS VCF, SECOND IS ANNOVAR, THIRD OUTPUT
import numpy as np
import re
import sys
#print sys.argv
#%% import and format VCF file
vcf=[]
c=0
#for line in open('foo','r'):
for line in open(sys.argv[1],'r'):
    line.strip()
    if c == 0:
        v_header=np.array(line.rstrip().split('\t'))
        c = 1
    else:
        line=re.sub('0/0[\S]+','0', line)
        line=re.sub('0/1[\S]+','1', line)
        line=re.sub('1/1[\S]+','2', line)
        line=re.sub('\./\.[\S]+','.', line)
        vcf.append(line.rstrip().split('\t'))
vcf=np.array(vcf)



#%% import and format ANNOVAR FILE
anno=[]
for line in open(sys.argv[2],'r'):
    line.strip()
    if line.startswith('Chr'):
        a_header=np.array(line.rstrip().split('\t'))
    else:
        anno.append(line.rstrip().split('\t'))
anno=np.array(anno)
keep = [0,1,3,4,6,8,10,12,14,18,19,22,24,26,27,29,31,33,35,36,44,46,48,52,53]
anno = anno[:,keep]
a_header = a_header[keep]
anno[anno == '.'] = np.nan
def fix_scores():
    anno[:,6] = 1-anno[:,6].astype('float') # Correct SIFT to 1-score
    anno[:,9][anno[:,10]=='N'] = 1-anno[:,9][anno[:,10]=='N'].astype('float') #MutationTaster 1-score for N
    anno[:,9][anno[:,10]=='P'] = 1-anno[:,9][anno[:,10]=='P'].astype('float') #MutationTaster 1-score for P
    anno[:,11] = 1-((anno[:,11].astype('float') - (-16.13))/(10.64-(-16.13))) #FATHMM 1-score
    anno[:,12] = 1-((anno[:,12].astype('float') - (-14))/(14-(-14))) #PROVEAN
    anno[:,14] = (anno[:,14].astype('float') - (-2))/(3-(-2)) #MetaSVM
    anno[:,17] = (anno[:,17].astype('float') - (-7.535037))/(35.788538-(-7.535037)) #CADD_raw
    anno[:,20] = (anno[:,20].astype('float') - (-12.3))/(6.17-(-12.3)) #GERP++_RS
    anno[:,21] = (anno[:,21].astype('float') - (-13.282))/(1.199-(-13.282)) #phyloP20way_mammalian
fix_scores()
anno = np.delete(anno,10,1)#remove MutationTaster prediction
a_header = np.delete(a_header,10,0)

 
_fsi_loc = np.where(anno[:,5]=='frameshift insertion')[0]
_fsd_loc = np.where(anno[:,5]=='frameshift deletion')[0]
anno[_fsi_loc, 6:22] = 1.0 # set 1 for FSI
anno[_fsd_loc, 6:22] = 1.0 # set 1 for FSD

_stg_loc = np.where(anno[:,5]=='stopgain')[0] #search and modify stopgain
for _x in _stg_loc:
    _k  = anno[_x, 6:22].astype(float)
    _mask = np.isnan(_k)
    _k[_mask] = 1.0
    anno[_x,6:22]= _k
    
_stl_loc = np.where(anno[:,5]=='stoploss')[0]  #search and modify stoploss
for _x in _stl_loc:
    _k  = anno[_x, 6:22].astype(float)
    _mask = np.isnan(_k)
    _k[_mask] = 1.0
    anno[_x,6:22]= _k

_ukn_var = np.where(anno[:,5]=='nan')[0] #remove uknown variants
anno = np.delete(anno, _ukn_var,0)
vcf = np.delete(vcf, _ukn_var,0)
_syn_var = np.where(anno[:,5]=='synonymous SNV')[0] #remove synonymous variants
anno = np.delete(anno, _syn_var,0)
vcf = np.delete(vcf, _syn_var,0)
K = np.hstack((anno,vcf))
k_header = np.hstack((a_header,v_header))

combi=np.vstack((k_header,K)).astype('string')

#np.savetxt('test_lev_combined.txt', combi, delimiter='\t', fmt='%s')
np.savetxt(sys.argv[3], combi, delimiter='\t', fmt='%s')
print 'DONE!'
