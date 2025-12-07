import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
import matplotlib as mpl

### set matplotlib settings 
# text settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
# font style
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']


comb = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_clin2_wMeds_filtered_wAge4.xlsx')
comb['time_event_month'] = comb['time_to_event'].apply(lambda x: x / 30.44)

# define sigs and meds 
sigs = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']
meds = comb.columns[97:138]

# read in results df for HR and qvalue 
res = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_allres_SigMedEnt.xlsx')
res2 = res[res['covariate'].isin(sigs)]


### define KM 
def KM(sig, treat, ent):
    sig_ent = comb[(comb[sig] > 0) & (comb[treat] == 1) & (comb['Entität'] == ent)]
    nosig_ent = comb[(comb[sig] == 0) & (comb[treat] == 1) & (comb['Entität'] == ent)]

    fig, ax = plt.subplots(figsize = (6.5,5.5))

    kmf_sigent = KaplanMeierFitter()
    kmf_sigent.fit(sig_ent['time_event_month'], sig_ent['event'], label=f'{sig} Active')
    kmf_sigent.plot(ci_show=True, ci_alpha= 0.15)

    kmf_nosigent = KaplanMeierFitter()
    kmf_nosigent.fit(nosig_ent['time_event_month'], nosig_ent['event'], label=f'{sig} Inactive')
    kmf_nosigent.plot(ci_show= True, ci_alpha= 0.15, ax=ax)

    plt.title(f'{treat} in {re.sub('Neoplasms', '', ent).strip()} Cancer: Survival Stratified by {sig}')
    ax.set_ylabel('Survival Probability', fontsize = 10)
    ax.set_xlabel('Time (months)', fontsize = 10)

    add_at_risk_counts(kmf_sigent, kmf_nosigent, ax=ax)
    
    risk_ax = plt.gcf().axes[-1]

    for label in risk_ax.get_xticklabels() + risk_ax.get_yticklabels():
        label.set_fontsize(14) 

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='both', labelsize=10)

    logrank = logrank_test(sig_ent['time_event_month'], nosig_ent['time_event_month'], sig_ent['event'], nosig_ent['event'], alpha=.99)
    p = logrank.p_value
    p_text = 'log-rank p < 0.001' if p < 0.001 else f'log-rank p = {p:.3f}'
    ax.text( 0.96, 0.79, p_text, transform=ax.transAxes, ha='right', va='top', fontsize=10)

    ax.legend(fontsize=10, loc='upper right', bbox_to_anchor= (0.97, 0.97), frameon=False)
    plt.tight_layout()

    return fig