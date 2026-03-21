
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib as mpl

### set matplotlib settings 

# text settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['text.usetex']  = False
# font style
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']


### bubble plots for mannwhitney entsplit

res = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_mut_sig_corr_mannwhitney/MH_mutsigcorr_mannwhitney_entsplit.xlsx')


# filter for significant results 
res2 = res[res['q'] < 0.1]
res2['Ent'] = res2['Ent'].str.replace('_Neoplasms', '', regex=True)

# add log2(p) columns to res2 
res2['-log2_qvals'] = -np.log2(res2['q'])
res2['-log10_qvals'] = -np.log10(res2['q'])

ents = res2['Ent'].unique()

#### create the plot like in steele publication 
cmap_steele = LinearSegmentedColormap.from_list("darkorange_white_deep_purple", ['#FF8C00', "white", "#4B0082"])
const = 130


for ent in ents: 
    res3 = res2[res2['Ent'] == ent]

    print(res3)
    # reorder signatures in res2 df 
    sord = ['CN1', 'CN2', 'CN9', 'CN12', 'CN13', 'CN14', 'CN16', 'CN17', 'CN18']
    sord2 = [sig for sig in sord if sig in res3['Signature'].unique()]
    sig_pos = {sig: i for i, sig in enumerate(sord2)}
    res3['Sig_num'] = res3['Signature'].map(sig_pos)

    gene_order = sorted(res3['Gene'].unique())
    res3['Gene'] = pd.Categorical(res3['Gene'], categories=gene_order, ordered=True)

    res3 = res3.sort_values(['Gene'])

    fig, ax = plt.subplots(figsize = (5,5))

    # center the color bar around 0 
    #color_vals = res3['FC']  # assuming already log2FC
    #vmax = np.nanmax(np.abs(color_vals))
    R = 3.0
    norm = mcolors.TwoSlopeNorm(vmin= -R, vcenter=0, vmax=R)

    scatter = ax.scatter(x = res3['Sig_num'], y = res3['Gene'], c = res3['FC'], s = res3['-log10_qvals'] * const, cmap=cmap_steele, norm=norm, alpha = 1, edgecolors='none')

    # set title 
    ax.set_title(ent, fontsize=18, pad=20)

    # add legend for -log10(q) size 
    sizes = [const, 3 * const, 5 * const]
    labels = ['1', '3', '5']

    for size, label in zip(sizes, labels):
        ax.scatter([], [], s=size, c= 'black', label=label)

    ax.set_xticks(list(sig_pos.values()))
    ax.set_xticklabels(sord2)

    ax.set_xlim(-0.5, len(sord2) - 0.5)
    ax.set_ylim(-0.5, len(gene_order) - 0.5)

    ax.tick_params(axis='x', colors='black', labelsize=16) 
    ax.tick_params(axis='y', colors='black', labelsize=16)
    ax.set_xticklabels(sord2, rotation=90)

    ax.spines['top'].set_linewidth(1.3)
    ax.spines['right'].set_linewidth(1.3)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)

    legend = ax.legend(title = '-log10(q)', bbox_to_anchor=(1.5, 0.8), frameon=False, labelspacing=2)
    legend.get_title().set_fontsize(13) 
    

    plt.tight_layout()
    plt.savefig(f'/Users/emilnetz/Desktop/cnv_sigs/MH/figures/MH_sigmutcorr_mannwhitney_entsplit_plots/MH_sigmutcorr_mannwhitney_{ent}_plot.pdf', bbox_inches="tight")
