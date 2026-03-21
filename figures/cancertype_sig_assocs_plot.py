import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import re
from matplotlib.colors import LinearSegmentedColormap

### set matplotlib settings 

# text settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
# font style
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']



### mannwhitney typesigcorr bubble plot

res = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_type_sig_corr_mannwhitney/MH_typesigcorr_mannwhitney_nofilt.xlsx')

res['Entity'] = res['Entity'].str.replace('_Neoplasms', '', regex=True)

# filter for significant results 
#res2 = res[res['q'] < 0.05]
res2 = res[res['q'] < 0.1]

# add log2(p) columns to res2 
res2['-log2_qvals'] = -np.log2(res2['q'])
res2['-log10_qvals'] = -np.log10(res2['q'])


# reorder signatures in res2 df 

sord = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']
sord2 = [sig for sig in sord if sig in res2['Signature'].unique()]
sig_pos = {sig: i for i, sig in enumerate(sord2)}
res2['Sig_num'] = res2['Signature'].map(sig_pos)

gene_order = sorted(res2['Entity'].unique())
res2['Gene'] = pd.Categorical(res2['Entity'], categories=gene_order, ordered=True)

res2 = res2.sort_values(['Entity'])

#### create the plot like in steele publication 
cmap_steele = LinearSegmentedColormap.from_list("darkorange_white_deep_purple", ['#FF8C00', "white", "#4B0082"])
const = 100

fig, ax = plt.subplots(figsize = (11,9))

# center the color bar around 0 
R = 3.0
norm = mcolors.TwoSlopeNorm(vmin= -R, vcenter=0, vmax=R)

scatter = ax.scatter(x = res2['Sig_num'], y = res2['Gene'], c = res2['FC'], s = res2['-log10_qvals'] * const, cmap=cmap_steele, norm=norm, alpha = 1, edgecolors='none')

# add legend for color 
cbar = plt.colorbar(scatter, ax=ax, shrink = 0.3, anchor=(0.0,0.15), aspect=5) 
cbar.set_label("log2(FC)")
cbar.ax.yaxis.label.set_size(13)
cbar.set_ticks([-2, 0, 2])

# add legend for -log2(q) size 
sizes = [const, 5 * const, 10 * const]
labels = ['1', '5', '10']

for size, label in zip(sizes, labels):
    ax.scatter([], [], s=size, c= 'black', label=label)

ax.set_xticks(list(sig_pos.values()))
ax.set_xticklabels(sord2)

ax.tick_params(axis='x', colors='black', labelsize=16) 
ax.tick_params(axis='y', colors='black', labelsize=16)
ax.set_xticklabels(sord2, rotation=90)


ax.spines['top'].set_linewidth(1.3)
ax.spines['right'].set_linewidth(1.3)
ax.spines['left'].set_linewidth(1.3)
ax.spines['bottom'].set_linewidth(1.3)

legend = ax.legend(title = '-log10(q)', bbox_to_anchor=(1.34, 0.8), frameon=False, labelspacing=2)
legend.get_title().set_fontsize(13)  

plt.tight_layout()

plt.savefig('/Users/emilnetz/Desktop/cnv_sigs/MH/figures/MH_typesigcorr_mannwhitney_0.1_plot.pdf', bbox_inches="tight")