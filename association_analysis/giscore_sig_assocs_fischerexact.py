import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control 

HRD = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/metadata/HRD_cohort_wSigs.xlsx')
sigs = ['CN1', 'CN2', 'CN3', 'CN9', 'CN11', 'CN17', 'CN18', 'CN20']

# binarize GI_score and sigs 
for sig in sigs:
    HRD[f'{sig}_bin'] = '' 
    for i in range(len(HRD)):
        if HRD.loc[i, sig] > 0:
            HRD.loc[i, f'{sig}_bin'] = 1
        else:
            HRD.loc[i, f'{sig}_bin'] = 0

HRD['GI_score_bin'] = ''
for i in range(len(HRD)):
    if HRD.loc[i, 'GI_score'] >= 83:
        HRD.loc[i, 'GI_score_bin'] = 1
    else:
        HRD.loc[i, 'GI_score_bin'] = 0

# define fisher exact 
sigs_bin = ['CN1_bin', 'CN2_bin', 'CN9_bin', 'CN11_bin', 'CN17_bin', 'CN18_bin', 'CN20_bin']

allres = []
for sig in sigs_bin: 
    print(sig)
    res = fisher(sig)
    OR = res.statistic
    p = res.pvalue
    allres.append([sig, OR, p])


resdf = pd.DataFrame(data=allres, columns= ['Sig', 'OR', 'p'])
resdf['q'] = false_discovery_control(resdf['p'])

resdf.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/sig_giscore/HRD_sig_giscore_fisherexact.xlsx')






