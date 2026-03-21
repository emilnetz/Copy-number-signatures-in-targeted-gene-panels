import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


##### fisherexact for sig per cancer type

meta = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/metadata/HRD_cohort_wSigs2.xlsx')

# create df cancer type as columns

df = pd.DataFrame(0, columns= list(meta['Entität'].unique()), index= meta['ID'].unique())

row_idx = df.index.get_indexer(meta['ID'])
col_idx = df.columns.get_indexer(meta['Entität'])


df.values[row_idx, col_idx] = 1

df['ID'] = df.index
df.reset_index(drop=True, inplace=True)

# merge mutation cbioportal with meta sheet for signatures 

df2 = pd.merge(df, meta, on='ID', how='left')


# binarise signatures
sigs = ['CN1', 'CN2', 'CN3', 'CN9', 'CN11', 'CN17', 'CN18', 'CN20']

for sig in sigs: 
    df2[sig] = df2[sig].apply(lambda x: 1 if x > 0 else 0)


# Fishers exact test for ent / signature  
def fisher(sig, ent):
    table = pd.crosstab(df2[ent], df2[sig])
    res = fisher_exact(table)
    return res

ents = meta['Entität'].unique()

allres = []
for sig in sigs: 
    for ent in ents: 
            df_sub = df2[df2[sig] == 1]
            num = sum(df_sub[ent])
            if num > 0: 
                res = fisher(sig, ent)
                OR = res.statistic
                p = res.pvalue
                allres.append([sig, ent, OR, p, num])


resdf = pd.DataFrame(data=allres, columns= ['Sig', 'Cancer Type', 'OR', 'p', 'cases'])

resdf['q'] = false_discovery_control(resdf['p'])

resdf[resdf['q'] < 0.05]

resdf.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_sigcancertype_distr/HRD_sig_cancertype_fisherexact.xlsx', index=False)


















meta = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/metadata/HRD_cohort_wSigs2.xlsx')

### get number of cancertypes per sig

sigs = ['CN1', 'CN2', 'CN3', 'CN9', 'CN11', 'CN17', 'CN18', 'CN20']

alldf = pd.DataFrame()
for sig in sigs: 
    sub = meta[meta[sig] > 0]
    sub.reset_index(drop=True, inplace=True)

    entcounts = sub['Entität'].value_counts()
    entcounts = entcounts[0:10]

    sub2 = sub.copy()
    for i in range(len(sub2)):
        if sub2.loc[i, 'Entität'] not in entcounts.index:
            sub2.loc[i, 'Entität'] = 'Other'

    counts = sub2['Entität'].value_counts()
    counts.index
    df = pd.DataFrame({'type':list(counts.index), 'count':counts.values})
    df['sig'] = sig

    sums = sum(df['count'])
    df['rel'] = df['count'] / sums * 100
    alldf = pd.concat([alldf, df], ignore_index=True)


alldf.to_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/HRD_sigcancertype_distr/HRD_cancertype_per_sig_distr.xlsx', index=False)
