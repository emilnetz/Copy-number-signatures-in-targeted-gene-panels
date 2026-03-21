import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
import matplotlib.markers as mmarkers
import matplotlib.patches as patches
from scipy.stats import binomtest
from scipy.stats import spearmanr
import matplotlib as mpl


### set matplotlib settings 

# text settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
# font style
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']




### concordance scatter MH fisherexact steele 
mh = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_mut_sig_corr_fisherexact/MH_mutsigcorr_fisherexact_entsplit.xlsx')
mh['Ent'] = mh['Ent'].str.replace('_Neoplasms', '', regex=True)

stel = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/steele/steele_mut_sig_corr.xlsx', sheet_name='SNV_indel driver gene assocs')
tcga = {
    "PRAD": "Prostatic",
    "BRCA": "Breast",
    "UCEC": "Endometrial",
    "UCS":  "Endometrial",
    "OV":   "Ovarian",
    "LUAD": "Lung",
    "LUSC": "Lung",
    "THCA": "Thyroid",
    "COAD": "Colonic",
    "READ": "Colorectal",
    "STAD": "Stomach",
    "PAAD": "Pancreatic",
    "SKCM": "Melanoma",
    "CHOL": "Cholangiocarcinoma",
}

stel2 = stel.copy()
stel2['tumour'] = stel2['tumour'].map(tcga)
stel2 = stel2[stel2['tumour'].notna()]
stel2['OR'] = np.log2(stel2['OR'])

mh2 = mh.copy()
mh2['OR'] = np.log2(mh2['OR'])
mh2 = mh2[['Sig', 'Ent', 'Gene', 'OR', 'p', 'q']].rename(columns={'OR':'mh_OR', 'p':'mh_p', 'q':'mh_q'})

stel3 = stel2[['CNsig', 'tumour', 'gene', 'OR', 'p', 'q']].rename(columns={'CNsig':'Sig', 'tumour':'Ent', 'gene':'Gene', 'OR':'stel_OR', 'p':'stel_p', 'q':'stel_q'})

merged = mh2.merge(stel3, on=['Ent', 'Gene', 'Sig'], how='inner')
merged2 = merged[merged['stel_q'] < 0.05]

merged2.to_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_mut_sig_corr_fisherexact/MH_steele_oncordance.xlsx', index=False)



# statistics  
merged2 = merged[merged['stel_q'] < 0.05]
merged2.reset_index(drop=True, inplace=True)

x = merged2['stel_OR']
y = merged2['mh_OR']
x.reset_index(drop=True, inplace=True)
y.reset_index(drop=True, inplace=True)

agree = np.sign(x) == np.sign(y)
dagree = np.sign(x) != np.sign(y)
n_agree = sum(agree)
n = len(merged2)
p_agree = n_agree / n * 100

inf_df = merged2[np.isinf(merged2['mh_OR']) | np.isinf(merged2['stel_OR'])]
inf_df.reset_index(drop=True, inplace=True)

inf_df.loc[np.isinf(inf_df['mh_OR']), 'mh_OR'] = 0
inf_df.loc[np.isinf(inf_df['stel_OR']), 'stel_OR'] = 0


# plot  
fig, ax = plt.subplots(figsize=(6,6))
R = 6
ax.axhline(0, lw=1, color='0.85')
ax.axvline(0, lw=1, color='0.85')
ax.plot([-R, R], [-R, R], ls='--', lw=1, color='0.5')

# adjust inf values 
inf_df.loc[np.isneginf(inf_df['mh_OR']), 'mh_OR'] = -R + 0.1
inf_df.loc[np.isposinf(inf_df['mh_OR']), 'mh_OR'] = R - 0.1
inf_df.loc[np.isinf(inf_df['stel_OR']), 'stel_OR'] = -R + 0.1


# show all points regardless of significance
ax.scatter(x[agree], y[agree], s=50, alpha=0.9, facecolor='black', edgecolors='none', label='same sign')
ax.scatter(x[~agree], y[~agree], s=50, alpha=0.9, facecolors='none', edgecolors='black', label='opposite sign')

half_circle1 = Path.arc(0, 180)#oben
half_circle2 = Path.arc(90, 270)#links
half_circle3 = Path.arc(180, 360) #unten
half_circle4 = Path.arc(270, 90)#rechts

ax.scatter(inf_df['stel_OR'][1], inf_df['mh_OR'][1], marker=half_circle3, s=50, alpha=0.9, facecolor='black', edgecolors='none')
ax.scatter(inf_df['stel_OR'][0], inf_df['mh_OR'][0], marker=half_circle1, s=50, alpha=0.9, facecolor='black', edgecolors='none')
ax.scatter(inf_df['stel_OR'][2:4], inf_df['mh_OR'][2:4], marker=half_circle1, s=50, alpha=0.9, facecolor='black', edgecolors='none')
ax.scatter(inf_df['stel_OR'][4], inf_df['mh_OR'][4], marker=half_circle4, s=50, alpha=0.9, facecolor='black', edgecolors='none')


ax.set_xlim(-R, R); ax.set_ylim(-R, R)

ax.set_xlabel('Steele - log2(OR)', size = 17, color='black')
ax.set_ylabel('MH-panel - log2(OR)', size = 17, color='black')

# stats summary
p_agree = binomtest(n_agree, n, p=0.5, alternative='greater').pvalue if n>0 else np.nan
rho, p_rho = spearmanr(x, y, nan_policy='omit')

ax.text(0.02, 0.98, f'Spearman ρ={rho:.2f}, p={p_rho:.2g}', transform=ax.transAxes, va='top', color='black')


ax.tick_params(axis='x', colors='black', labelsize=13) 
ax.tick_params(axis='y', colors='black', labelsize=13)

square1 = patches.Rectangle((-R, -R), R, R, facecolor="lightblue", alpha=0.3)
square2 = patches.Rectangle((0,0), R, R, facecolor="lightblue", alpha=0.3)

ax.add_patch(square1)
ax.add_patch(square2)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.3)
ax.spines['bottom'].set_linewidth(1.3)
ax.legend(frameon=False, loc='lower right')
plt.tight_layout()
#plt.show()


plt.savefig('/Users/emilnetz/Desktop/cnv_sigs/MH/figures/MH_steele_mutsigcorr_fisherexact_concordance_plot20260130.pdf', bbox_inches="tight")


