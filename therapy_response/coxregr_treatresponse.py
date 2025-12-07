import numpy as np
import pandas as pd
from collections import Counter
from lifelines import CoxPHFitter
from scipy.stats import false_discovery_control


comb = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_clin2_wMeds_filtered_wAge4.xlsx')

entcounts = comb['Entität'].value_counts()
entcounts = entcounts[0:10]

# cancer types below the 10 most comon are combined to 'other' 
comb2 = comb.copy() 
for i in range(len(comb2)):
    if comb2.loc[i, 'Entität'] not in entcounts.index:
        comb2.loc[i, 'Entität'] = 'Other'

# define ents, meds and sigs
ents = comb2['Entität'].unique() 
meds = comb2.columns[97:138]
sigs = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']

# define cox function 
def cox (sig, med, entity):
    clin_sub = cocomb[(cocomb[med] == 1)]
    colnames = [sig, 'event', 'time_to_event', 'stage', 'age']
    clin_sub2 = clin_sub[colnames]
    cph = CoxPHFitter()
    cph.fit(clin_sub2, duration_col='time_to_event', event_col='event', show_progress=True)
    res = cph.summary
    return res

####

cox_df = pd.DataFrame()
for ent in ents: 
    # combine treatments when coadministered in > 80% of cases
    cocomb = pd.read_excel(f'/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/cotreatment/clin_wCotreatmentGroups2/MH_{ent}_clin_wCotreatmentGroups.xlsx')
    meds2 = [med for med in meds if med in cocomb.columns]
    pos = cocomb.columns.get_loc('age')
    cos = cocomb.iloc[:,pos+1:].columns
    meds3 = meds2.copy()
    for el in cos: 
        comeds = el.split('_')
        for comed in comeds:
            try: 
                meds3.remove(comed)
            except Exception as e: 
                print(f'{comed} not removed')
            continue  
    meds4 = meds3 + list(cos)

    for sig in sigs: 
        for med in meds4:
            subdf = cocomb[cocomb[med] == 1]
            pos = sum(subdf[sig] > 0) 
            neg = len(subdf) - pos 
            if (pos >= 10) & (neg >= 10):     
                try:
                    res = cox(sig, med, ent)
                    res['Signature'] = sig
                    res['Therapy'] = med
                    res['Entitiy'] = ent
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
            

cox_df['q'] = false_discovery_control(cox_df['p'])

cox_df.to_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_allres_SigMedEnt_wCoGroups.xlsx')