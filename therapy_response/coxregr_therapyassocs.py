import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import false_discovery_control



comb = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_clin2_wMeds_filtered_wAge4_folfadjust20260212.xlsx')


# everything after the 10 most comon will be OTHER 
entcounts = comb['Entität'].value_counts()
entcounts = entcounts[0:10]

comb2 = comb.copy() 
for i in range(len(comb2)):
    if comb2.loc[i, 'Entität'] not in entcounts.index:
        comb2.loc[i, 'Entität'] = 'Other'


# define ents, meds and sigs
ents = comb2['Entität'].unique() 
ents = list(ents)
ents.remove('Other')
sigs = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']

#define med variable 
pos1 = comb.columns.get_loc('Entität')
pos2 = comb.columns.get_loc('age')
meds = comb.columns[pos1+1:pos2]



# everything after the 10 most comon will be OTHER 
def cox (sig, med):
    cols = [sig, 'event', 'time_to_event', 'stage', 'age']
    clin_sub = cocomb.loc[(cocomb[med] == 1), cols]
    cph = CoxPHFitter()
    cph.fit(clin_sub, duration_col='time_to_event', event_col='event', show_progress=True)
    res = cph.summary
    return res


cox_df = pd.DataFrame()
for ent in ents: 

    cocomb = pd.read_excel(f'/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/cotreatment/cotreatment_adjfolf_20260220/MH_{ent}_clin_wCotreatmentGroups.xlsx')

    #define single meds
    pos1 = cocomb.columns.get_loc('Entität')
    pos2 = cocomb.columns.get_loc('age')
    meds =  cocomb.columns[pos1+1:pos2]

    # define combined meds
    pos3 = cocomb.columns.get_loc('age')
    cos = cocomb.iloc[:,pos3+1:].columns
        
    meds2 = list(meds) + list(cos)

    for sig in sigs: 
        for med in meds2:
            subdf = cocomb[cocomb[med] == 1]
            pos = sum(subdf[sig] > 0) 
            neg = len(subdf) - pos 
            if ((pos >= 10) & (neg >= 10)):  
                try:
                    res = cox(sig, med)
                    res['Signature'] = sig
                    res['Therapy'] = med
                    res['Entity'] = ent
                    subdf = cocomb[(cocomb['Entität'] == ent) & (cocomb[med] == 1)]
                    pos = sum(subdf[sig] > 0) 
                    neg = len(subdf) - pos
                    res['pos_cases'] = pos
                    res['neg_cases'] = neg
                    res['total_cases'] = pos + neg
                    cox_df = pd.concat([cox_df, res], axis=0)
            
                except Exception as e: 
                    print(f'Signature {sig}, {med} and {ent} didnt work, skipping...')
                continue 
                

cox_df2 = cox_df[~cox_df.index.isin(['stage', 'age'])]
cox_df2['q'] = false_discovery_control(cox_df2['p'], method='bh')
cox_df2.to_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_allres_SigMedEnt_wCoGroups20260220_2.xlsx')




### Cox regression with binary CN sigs 

MH = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_clin2_wMeds_filtered_wAge4_folfadjust20260212.xlsx')

# define 'other' category for seldom entities 
comb = MH.copy()

entcounts = comb['Entität'].value_counts()
entcounts = entcounts[0:10]

# everything after the 10 most comon will be OTHER 
comb2 = comb.copy() 
for i in range(len(comb2)):
    if comb2.loc[i, 'Entität'] not in entcounts.index:
        comb2.loc[i, 'Entität'] = 'Other'


# define ents, meds and sigs
ents = comb2['Entität'].unique() 
ents = list(ents)
ents.remove('Other')
sigs = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']

#define med variable 
pos1 = comb.columns.get_loc('Entität')
pos2 = comb.columns.get_loc('age')
meds = comb.columns[pos1+1:pos2]



# everything after the 10 most comon will be OTHER 
def cox (sig, med):
    cols = [sig, 'event', 'time_to_event', 'stage', 'age']
    clin_sub = cocomb.loc[(cocomb[med] == 1), cols]
    cph = CoxPHFitter()
    cph.fit(clin_sub, duration_col='time_to_event', event_col='event', show_progress=True)
    res = cph.summary
    return res


cox_df = pd.DataFrame()
for ent in ents: 

    cocomb = pd.read_excel(f'/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/cotreatment/cotreatment_adjfolf_20260220/MH_{ent}_clin_wCotreatmentGroups.xlsx')
    cocomb[sigs] = cocomb[sigs].applymap(lambda x: 0 if x == 0 else 1)
    #define single meds
    pos1 = cocomb.columns.get_loc('Entität')
    pos2 = cocomb.columns.get_loc('age')
    meds =  cocomb.columns[pos1+1:pos2]

    # define combined meds
    pos3 = cocomb.columns.get_loc('age')
    cos = cocomb.iloc[:,pos3+1:].columns
        
    meds2 = list(meds) + list(cos)

    for sig in sigs: 
        for med in meds2:
            subdf = cocomb[cocomb[med] == 1]
            pos = sum(subdf[sig] > 0) 
            neg = len(subdf) - pos 
            if ((pos >= 10) & (neg >= 10)):  
                try:
                    res = cox(sig, med)
                    res['Signature'] = sig
                    res['Therapy'] = med
                    res['Entity'] = ent
                    subdf = cocomb[(cocomb['Entität'] == ent) & (cocomb[med] == 1)]
                    pos = sum(subdf[sig] > 0) 
                    neg = len(subdf) - pos
                    res['pos_cases'] = pos
                    res['neg_cases'] = neg
                    res['total_cases'] = pos + neg
                    cox_df = pd.concat([cox_df, res], axis=0)
            
                except Exception as e: 
                    print(f'Signature {sig}, {med} and {ent} didnt work, skipping...')
                continue 
                

cox_df2 = cox_df[~cox_df.index.isin(['stage', 'age'])]
cox_df2['q'] = false_discovery_control(cox_df2['p'], method='bh')
cox_df2.to_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_allres_SigMedEnt_wCoGroups20260227_binary.xlsx')