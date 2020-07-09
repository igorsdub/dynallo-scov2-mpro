# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 10:43:15 2020

@author: Igors Dubanevics
"""

# Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl


#%% Helpers fucntions

def globalMap(dataFrame, start, end):
    """ Plots global map of allostery for 4hzf (CAP) as in Rodgers et al PLoS paper. 
    Figure is saved as PNG. Title conatins number of modes used to calculate the map.
    
    >>> globalMap(df_100modes, 100)
    """
    
    # Set axis and title
    xlbl='Amino Acid Number'
    ylbl=r'$k_R/k$'
    
    # Set style
    sns.set_style("ticks")
    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = sns.heatmap(dataFrame,center = dataFrame.loc[1.0,10], yticklabels = 2, cmap=plt.cm.RdBu_r, 
                     cbar_kws={'label': r"$K_2 / K_1$"})
    
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
    
    return ax, fig

#%% Variables
res_start = 1
res_end = 306


aa_num = res_end - res_start + 1 
aa_list = list(range(res_start,res_end+1))

#%% Make dataframe
spring_str = np.loadtxt('spring_strength.txt')
spring_num = len(spring_str)

# Font scale and line width
color_line = 'black'
f = 2.5
l = 3
sns.set_style("ticks")
sns.set(rc={'figure.figsize':(13,8)})
sns.set_context("notebook", font_scale=f, rc={"lines.linewidth": l})

# Dataframe
# Calculate dissocitaion constant ratio
k2_k1 = np.exp(ddG_cum[mode_cust])
# Construct data frame
df_cust = pd.DataFrame(k2_k1,columns=aa_list)
df_cust.insert(loc=0, column='k_R/k', value=spring_str)
df_cust = df_cust.set_index('k_R/k')

ax, fig = globalMap(df_cust, res_start, res_end)

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

fig_path = '1point-allo_'+str(mode_cust)+'modes'
fig.savefig(fig_path+'.png', dpi=600, bbox_inches='tight', format='png')
# fig.savefig(fig_path+'.eps', dpi=300, bbox_inches='tight', format='eps')
fig.savefig(fig_path+'.svg', dpi=600, bbox_inches='tight', format='svg')
# fig.savefig(fig_path+'.pdf', dpi=600, bbox_inches='tight', format='pdf')