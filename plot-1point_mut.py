# -*- coding: utf-8 -*-
"""
This script plots 1-point mutational maps of relative fluactuation free energy change and global allostery for 6lu7 (Mpro) ENM.

Written by:
            Igors Dubanevics (id583 at york dot ac dot uk)
            Prof. Tom McLeish Group
            University of York
            Apr, 2020

"""
# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# I/O files
data_apo = "6lu7_1D-apo.xlsx"
data_allo = "6lu7_1D-allo.xlsx"

outname_apo = "6lu7_1D-map_apo"
outname_allo = "6lu7_1D-map_allo"

def globalMap(dataFrame, cbar_lbl):
    """ Plots 1-point mutational map for 6lu7 (Mpro) ENM. 
    
    >>> globalMap(df_apo, lbl_apo)
    """
    # Set style
    sns.set_style("ticks")
    # Font scale and line width
    color_line = 'black'
    f = 2.5
    l = 3
    sns.set_style("ticks")
    sns.set(rc={'figure.figsize':(13,8)})
    sns.set_context("notebook", font_scale=f, rc={"lines.linewidth": l})

    # Set axis and title
    xlbl='Amino Acid Number'
    ylbl=r'$k_R/k$'

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = sns.heatmap(dataFrame,center = dataFrame.loc[1.0,10], yticklabels = 2, cmap=plt.cm.RdBu_r, 
                     cbar_kws={'label': cbar_lbl})
        
    start = 1
    end = 306

    # x-axis ticks
    l = list(np.subtract(np.arange(start,end+50,50), start))
    if l[0] < 50:
        del l[0]
    if l[-1] > end:
        del l[-1]  
    xlabels = tuple(map(str, l))
    xpositions = tuple(np.subtract(l, start-1))
    ax.set(xticks=xpositions, xticklabels=xlabels)
    
    ax.set_ylim(0, 22)
    ypositions = (0, 5, 8, 10, 14, 18, 22)
    ylabels = ("0.25", "0.50", "0.75", "1.00", "2.00", "3.00", "4.00")
    ax.set(yticks=ypositions, yticklabels=ylabels)
    
    fig = ax.get_figure()

    # ax.invert_yaxis()
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)
    
    # Twin x-axis
    ax2 = ax.twiny()
    # Active residues
    res_act = np.concatenate((np.arange(41,42), np.arange(49,50), np.arange(143,146), np.arange(163,168), np.arange(187,193), np.arange(214,215), np.arange(284,287)))
    # Set x-axis limit
    ax2.set_xlim(0, aa_num-1)
    # Set tick style
    plt.style.use('seaborn')
    forest = (0.2,0.6,0.2)
    ax2.tick_params(axis="x", direction="out", length=5*f, width=f, color=forest)
    # Hide tick labels
    ax2.set_xticks(res_act)
    ax2.set_xticklabels([])

    return fig

#%% Variables
res_start = 1
res_end = 306

aa_num = res_end - res_start + 1 
aa_list = list(range(res_start,res_end+1))

#%% Make dataframe
#spring_str = ['0.25', '0.30', '0.35', '0.40', '0.45', '0.5', '0.58', '0.67', '0.75', '0.88', '1.00', '1.25', '1.50', '1.75', '2.00', '2.25', '2.50', '2.75', '3.00', '3.25', '3.50', '3.75', '4.00']

df_apo = pd.read_excel(data_apo, header=0, index_col=0)
df_allo = pd.read_excel(data_allo, header=0, index_col=0)

lbl_apo = r"$(G_{mut}-G_{wt})/|G_{wt}|}$ [ $\%$ ]"
lbl_allo = r"$K_2 / K_1$"

fig = globalMap(df_apo, lbl_apo)
fig.savefig(outname_apo+'.png', dpi=600, bbox_inches='tight', format='png')
plt.close('all')

fig = globalMap(df_allo, lbl_allo)
fig.savefig(outname_allo+'.png', dpi=600, bbox_inches='tight', format='png')
plt.close('all')