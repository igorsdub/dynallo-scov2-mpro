# -*- coding: utf-8 -*-
"""
This script plots 2-point mutational maps of relative fluactuation free energy change and global allostery for 6lu7 (SARS-CoV-2 Main Proteasse) ENM.

Written by:
            Igors Dubanevics (id583 at york dot ac dot uk)
            Prof. Tom McLeish Group
            University of York
            Apr, 2020

"""
# Import packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

inname_freeEn_wt = "6lu7_freeEn_wt.xlsx"

innames_apo = ("6lu7_2D-apo_kRk=0_25.xlsx","6lu7_2D-apo_kRk=4_00.xlsx")
innames_allo = ("6lu7_2D-allo_kRk=0_25.xlsx","6lu7_2D-allo_kRk=4_00.xlsx")

outnames_apo = ("6lu7_2D-apo_map_kRk=0_25","6lu7_2D-apo_map_kRk=4_00")
outnames_allo = ("6lu7_2D-allo_map_kRk=0_25","6lu7_2D-allo_map_kRk=4_00")

lbl_apo = r"$(G_{mut}-G_{wt})/|G_{wt}|}$ [ $\%$ ]"
lbl_allo = r"$K_2 / K_1$"

# Total number of modes, amino acids and pairwise mutations
mode_num = 100
aa_num = 612

freeEn_wt = pd.read_excel(inname_freeEn_wt, header=0, index_col=0)  
freeEn_sum_wt = freeEn_wt[:25].sum(axis=0)

#%% Fill 2D mutational array
alloFreeEn_wt = freeEn_sum_wt[2] - 2 * freeEn_sum_wt[1] + freeEn_sum_wt[0] 
k2_k1_wt = np.exp(alloFreeEn_wt) 

def globalMap(dataframe, lbl, ctr):
    xlbl='Amino Acid Number'
    ylbl='Amino Acid Number'

    # Font scale and line width
    color_line = 'black'
    f = 2.5
    l = 0.5
    sns.set_style("ticks")
    sns.set(rc={'figure.figsize':(13,8)})
    sns.set_context("notebook", font_scale=f, rc={"lines.linewidth": l})


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax = sns.heatmap(dataframe.iloc[::-1], center = ctr, square=True, cmap=plt.cm.RdBu_r, cbar_kws={'label': lbl})
    fig = ax.get_figure()

    # Set axis labels
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    plt.xticks(rotation=0)
    ax.invert_yaxis()
    
    # Plot active residues
    chain_len = int(aa_num/2)
    positions = ( 99, 199, 299, 99+chain_len, 199+chain_len, 299+chain_len)
    labels = ( "100", "200", "300", "100", "200", "300")
    ax.set(yticks=positions, yticklabels=labels)
    ax.set(xticks=positions, xticklabels=labels)

    line_pos = [aa_num/2]
    plot_end = aa_num - 1

    # Draw vertical and horizontal lines
    ax.hlines(line_pos, *ax.get_xlim(), color=color_line)
    ax.vlines(line_pos, *ax.get_ylim(), color=color_line)
    # Draw diagonal lines to show 1D scan area
    ax.plot([0.5, 1], [0, 0.5], transform=ax.transAxes, ls='--', color=color_line)
    ax.plot([0, 0.5], [0.5, 1], transform=ax.transAxes, ls='--', color=color_line)

    # Create second axis
    ax2 = fig.add_axes(ax.get_position())
    ax2.set_facecolor("None")
    ax2.set_aspect('equal')
    ax2.tick_params(bottom=0, top=1, left=0, right=1, 
            labelbottom=0, labeltop=1, labelleft=0, labelright=1)
    # Active residues
    res_act = np.concatenate((np.arange(41,42), np.arange(49,50), np.arange(143,146), np.arange(163,168), np.arange(187,193), np.arange(214,215),np.arange(284,287)))
    tick_pos = [*(res_act), *(res_act+chain_len)]

    # Set axis limits
    ax2.set_xlim(0, plot_end)
    ax2.set_ylim(0, plot_end)
    # Set tick style
    plt.style.use('seaborn')
    forest = (0.2,0.6,0.2)
    ax2.tick_params(axis="x", direction="out", length=5*f, width=f, color=forest)
    ax2.tick_params(axis="y", direction="out", length=5*f, width=f, color=forest)
    ax2.grid(False)
    # Hide tick labels
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels([])
    ax2.set_yticks(tick_pos)
    ax2.set_yticklabels([])

    ax2.hlines([0,plot_end], *ax2.get_xlim())
    ax2.vlines([0,plot_end], *ax2.get_ylim())
    
    return fig

for i in range(2):
    array_apo = None
    array_allo = None
    fig = None
    
    array_apo = pd.read_excel(innames_apo[i], header=0, index_col=0)
    fig = globalMap(array_apo, lbl_apo, 0)
    plt.savefig(outnames_apo[i]+'.png',dpi=600, format='png')
    plt.close('all')

    array_allo = pd.read_excel(innames_allo[i], header=0, index_col=0)
    fig = globalMap(array_allo, lbl_allo, k2_k1_wt)
    plt.savefig(outnames_allo[i]+'.png',dpi=600, format='png')
    plt.close('all')

print('Figures saved.')
