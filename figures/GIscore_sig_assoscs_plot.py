import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import re
from adjustText import adjust_text
import matplotlib as mpl

### set matplotlib settings 

# text settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
# font style
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']



### create scarsigcorr

res = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/HRD/sig_giscore/HRD_sig_giscore_fisherexact3.xlsx', index_col=0)

res['log2_OR'] = np.log2(res['OR'])
res['-log10_qvals'] = -np.log10(res['q'])

# rename Sigs 
res['Sig'] = res['Sig'].apply(lambda x: re.sub('_bin', '', x))

# reorder signatures in res2 df 
order = ['CN1', 'CN2', 'CN9', 'CN11', 'CN17', 'CN18', 'CN20']
res['Sig'] = pd.Categorical(res['Sig'], categories=order, ordered=True)
res = res.sort_values('Sig')


# plot 
fig, ax = plt.subplots(figsize = (7,5))
colors = res['Sig'].apply(lambda sig:  '#D55E00' if ((sig == 'CN17') | (sig == 'CN18')) else '#0072B2' if ((sig == 'CN1') | (sig == 'CN2')) else 'grey')
ax.scatter(x=res['log2_OR'], y=res['-log10_qvals'], s = 110, alpha = 0.65, edgecolors='none', color=colors)
for i, row in res.iterrows():
   if i < 5: 
      ax.text(row['log2_OR']+0.1, row['-log10_qvals'], row['Sig'], fontsize=15, ha='left', va='bottom')
   elif i == 5: 
     ax.text(row['log2_OR'] + 0.1 , row['-log10_qvals'], row['Sig'], fontsize=15, ha='left', va='bottom')
   elif i == 6: 
     ax.text(row['log2_OR'] + 0.1, row['-log10_qvals'], row['Sig'], fontsize=15, ha='left', va='bottom')

ax.set_xlim([-3.2, 3.2])
ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth = 1.2)
ax.axvline(0, color="grey", linestyle="--", linewidth= 1.2)
#ax.axvspan( 0, 4, color='red', alpha= 0.1)
ax.text(0.7, 6.5, 'GI-score enriched', fontweight='bold', size= 14)
ax.text(-2.5, 6.5, 'GI-score depleted', fontweight='bold', size=14)

ax.tick_params(axis='x', colors='black', labelsize=13) 
ax.tick_params(axis='y', colors='black', labelsize=13)

ax.set_xlabel('log2(OR)', size = 17, color='black')
ax.set_ylabel('-log10(q)', size=17, color='black')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.3)
ax.spines['bottom'].set_linewidth(1.3)

#adjust_text(texts, clip_on=False)
plt.tight_layout()
#plt.show()

plt.savefig('/Users/emilnetz/Desktop/cnv_sigs/HRD/figures/HRD_GIscoresigcorr_plot3.pdf', bbox_inches="tight")
