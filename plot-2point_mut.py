# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:23:08 2020

@author: lolo
"""
# Import packages
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


outname = '2D-allo'

# Total number of modes, amino acids and pairwise mutations
mode_num = 100
aa_num = 612
mut_num = int(aa_num * (aa_num-1) / 2) + aa_num


   
for i in range(3):
    infile_Gwt = 'results-'+str(i)+'-wt.txt'
    inlines_Gwt = np.loadtxt(infile_Gwt)
    G_allo_wt.append(inlines_Gwt)
#%% Fill 2D mutational array
ddGwt_sum = np.sum(G_allo_wt[2][:mode_sum] - 2 * G_allo_wt[1][:mode_sum] + G_allo_wt[0][:mode_sum]) 

xlbl='Amino Acid Number'
ylbl='Amino Acid Number'

# Font scale and line width
color_line = 'black'
f = 2.5
l = 0.5
sns.set_style("ticks")
sns.set(rc={'figure.figsize':(13,8)})
sns.set_context("notebook", font_scale=f, rc={"lines.linewidth": l})


k2_k1 = np.exp(ddG_2D)
k2_k1wt = np.exp(ddGwt_sum)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax = sns.heatmap(k2_k1, center = k2_k1wt, square=True, cmap=plt.cm.RdBu_r, cbar_kws={'label': r'$K_2 / K_1$'})
fig = ax.get_figure()

# Set axis labels
ax.set_xlabel(xlbl)
ax.set_ylabel(ylbl)

# Invert y-axis
ax.invert_yaxis()
plt.xticks(rotation=0)
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

plt.savefig(outname+'_'+str(mode_sum)+'modes.png',dpi=600, format='png')