from scipy.stats import mannwhitneyu
from scipy.stats import false_discovery_control
import pandas as pd
import numpy as np

muts = pd.read_csv('/Users/emilnetz/Desktop/cnv_sigs/HRD/cbioportal/HRD_cbiop_muts_wCNA3.csv', sep='\t')
meta = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/metadata/HRD_cohort_wSigs.xlsx')

genes= muts['gene'].unique()
ents = meta['Entität'].unique()


### define FC function
def FoldChange(act_mut, act_no_mut):
    pesudocount = 1e-6
    FC = (np.log2((np.mean(act_mut)+pesudocount)/ (np.mean(act_no_mut)+ pesudocount)))
    return FC 


### for homdel
muts2 = muts[muts['variant'] == 'HOMDEL']
genes= muts2['gene'].unique()

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

df.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_cna_sig_corr_mannwhitney/HRD_cnasigcorr_mannwhitney_HOMDEL_entsplit.xlsx', index=False)



### for amplifications
muts2 = muts[muts['variant'] == 'AMP']
genes= muts2['gene'].unique()

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
            if len(mutsubdf) >= 5 and len(notmutsubdf) >= 5:
                U1, p = mannwhitneyu(act1, act2, alternative = 'two-sided')
                fc = FoldChange(act1, act2)
                allres.append([sig, ent, gene, fc, p, len(mutsubdf), len(notmutsubdf)])
            
df = pd.DataFrame(data=allres, columns= ['Signature', 'Ent', 'Gene', 'FC', 'p', 'num mut', 'num not mut'])
df['q'] = false_discovery_control(df['p'])

df.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_cna_sig_corr_mannwhitney/HRD_cnasigcorr_mannwhitney_AMP_entsplit.xlsx', index=False)


### for ampflifications + gain
muts2 = muts[muts['variant'].isin(['AMP', 'GAIN'])]
genes= muts2['gene'].unique()

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
            if len(mutsubdf) >= 5 and len(notmutsubdf) >= 5:
                U1, p = mannwhitneyu(act1, act2, alternative = 'two-sided')
                fc = FoldChange(act1, act2)
                allres.append([sig, ent, gene, fc, p, len(mutsubdf), len(notmutsubdf)])
            
df = pd.DataFrame(data=allres, columns= ['Signature', 'Ent', 'Gene', 'FC', 'p', 'num mut', 'num not mut'])
df['q'] = false_discovery_control(df['p'])

df.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_cna_sig_corr_mannwhitney/HRD_cnasigcorr_mannwhitney_AMPGAIN_entsplit.xlsx', index=False)
