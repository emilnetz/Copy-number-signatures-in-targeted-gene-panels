from scipy.stats import mannwhitneyu
from scipy.stats import false_discovery_control
import pandas as pd
import numpy as np


muts = pd.read_csv('/Users/emilnetz/Desktop/cnv_sigs/HRD/cbioportal/HRD_cbiop_muts_wCNA3.csv', sep='\t')
muts.columns = ['ID', 'gene', 'variant', 'type']

muts2 = muts[muts['type'] != 'CNA']

meta = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/metadata/HRD_cohort_wSigs2.xlsx')
meta['Entität'].unique()

sigs = ['CN1', 'CN2', 'CN3', 'CN9', 'CN11', 'CN17', 'CN18', 'CN20']
genes= muts2['gene'].unique()
muts = muts2['type'].unique()
ents = meta['Entität'].unique()

### define FC function
def FoldChange(act_mut, act_no_mut):
    pesudocount = 1e-6
    FC = (np.log2((np.mean(act_mut)+pesudocount)/ (np.mean(act_no_mut)+ pesudocount)))
    return FC 

### loop over all 
allres = []
for sig in sigs:
    for ent in ents:  
        for gene in genes: 
            subdf1 = muts2[(muts2['gene'] == gene)]
            ids1 = subdf1['ID'].unique()
            mutsubdf = meta[(meta['ID'].isin(ids1)) & (meta['Entität'] == ent)]
            notmutsubdf = meta[(~meta['ID'].isin(ids1)) & (meta['Entität'] == ent)]
            act1 = mutsubdf[sig]
            act2 = notmutsubdf[sig]
            if len(mutsubdf) >= 10 and len(notmutsubdf) >= 10:
                U1, p = mannwhitneyu(act1, act2, alternative = 'two-sided')
                fc = FoldChange(act1, act2)
                allres.append([sig, ent, gene, fc, p, len(mutsubdf), len(notmutsubdf)])
            
df = pd.DataFrame(data=allres, columns= ['Signature', 'Ent', 'Gene', 'FC', 'p', 'num mut', 'num not mut'])

df['q'] = false_discovery_control(df['p'])

df.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_mut_sig_corr_mannwhitney/HRD_mutsigcorr_mannwhitney_entsplit.xlsx', index=False)




### comboine pvals where FCs point in the same direction 

resdf = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_mut_sig_corr_mannwhitney/HRD_mutsigcorr_mannwhitney_entsplit.xlsx')

genes = resdf['Gene'].unique()
sigs = resdf['Signature'].unique()

allres = []
for gene in genes: 
    for sig in sigs: 
        subdf = resdf[(resdf['Gene'] == gene) & (resdf['Signature'] == sig)]
        subdf.reset_index(drop=True, inplace=True)
        FCs = list(subdf['FC'])
        lenfc = len(FCs)
        if lenfc > 0: 
            if min(FCs) > 0: 
                res = combine_pvalues(subdf['p'])
                allres.append([gene, sig, res.pvalue, 'pos', lenfc])
            elif max(FCs) < 0: 
                res = combine_pvalues(subdf['p'])
                allres.append([gene, sig, res.pvalue, 'neg', lenfc])


combpdf = pd.DataFrame(data=allres, columns= ['Gene','Sig', 'p', 'direction', 'number of assocs'])
combpdf['q'] = false_discovery_control(combpdf['p'])

combpdf.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_mut_sig_corr_mannwhitney/HRD_mutsigcorr_mannwhitney_combined_pvals.xlsx', index=False)
