import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

forms = ("apo","holo1","holo2")

def crosscor_dist(form):
    # I/O files
    infile_dist = '6lu7_dist_'+ form +'.xlsx'
    infile_cros = '6lu7_cross-cor_'+ form +'.xlsx'
    
    aa_num = 612
    het_num = 21
    
    xlbl='Amino Acid Number'
    ylbl='Amino Acid Number'

    df_dist= pd.read_excel(infile_dist, header=0, index_col=0)
    df_cros= pd.read_excel(infile_cros, header=0, index_col=0)

    array_dist = np.flipud(df_dist.to_numpy())
    array_cros = np.flipud(df_cros.to_numpy())
    
    sns.set()
    
    # Create masking matrix
    mask_dist = np.triu(array_dist)
    mask_cros = np.tril(array_cros)
    
    color_line = 'black'
    f = 1
    l = 1
    sns.set_style("ticks")
    sns.set(rc={'figure.figsize':(13,8)})
    sns.set_context("notebook", font_scale=f, rc={"lines.linewidth": l})
    
    nice_plot ={
        'font.size': 30,
        'xtick.labelsize': 22,
        'ytick.labelsize': 22,
        'figure.autolayout': False,
        'axes.titlesize' : 22,
        'axes.labelsize' : 26,
        'lines.linewidth' : 0.5,
        'lines.markersize' : 6,
        'legend.fontsize': 13
        }
    
    plt.rcParams.update(nice_plot)
    
    fig=plt.figure(1)
    ax=fig.add_subplot(111)

    # Plot distance map
    ax1 = sns.heatmap(array_dist, vmin=0,vmax=16, mask=mask_dist, square=True,cmap=plt.cm.gist_yarg_r, 
                    cbar_kws={'label':'Distance [ $\AA{}$ ]','ticks':np.arange(0,18,2)})
    # Plot cross-correlation map 
    ax2 = sns.heatmap(array_cros,  vmin=-1,vmax=1, mask=mask_cros, square=True,cmap=plt.cm.RdBu_r, 
                    cbar_kws={'label': 'Cross-correlation','ticks':np.arange(-1.00,1.25,0.25)})
    
    
    ax.set_xlabel(xlbl)
    
    ax.set_ylabel(ylbl)
    
    plt.gca()
    
    # Invert y-axis
    plt.gca().invert_yaxis()
    plt.xticks(rotation=0)
    
    positions = (99, 199, 299, 407, 507, 607)
    labels = ("100", "200", "300", "100", "200", "300")
    ax.set(yticks=positions, yticklabels=labels)
    ax.set(xticks=positions, xticklabels=labels)
    
    if form == "apo":
        line_pos = [0,array_dist.shape[1], aa_num/2]
        plot_end = aa_num
    elif form == "holo1":
        line_pos = [0,array_dist.shape[1],aa_num/2, aa_num]
        plot_end = aa_num + het_num
    else:
        line_pos = [0,array_dist.shape[1],aa_num/2, aa_num, aa_num + het_num]
        plot_end = aa_num + 2*het_num
        
    ax.hlines(line_pos, *ax.get_xlim())
    ax.vlines(line_pos, *ax.get_ylim())

    
    # Twin x-axis
    ax2 = fig.add_axes(ax.get_position())
    ax2.set_facecolor("None")
    ax2.set_aspect('equal')
    ax2.tick_params(bottom=0, top=1, left=0, right=1, 
            labelbottom=0, labeltop=1, labelleft=0, labelright=1)
    
    # Active residues
    res_act = np.concatenate((np.array([41,49,214]), np.arange(143,146), np.arange(163,168), np.arange(187,193), np.arange(284,287)))
    tick_pos = [*(res_act), *(res_act + aa_num/2)]
    
    # Set x-axis limit
    ax2.set_xlim(0, plot_end)
    ax2.set_ylim(0, plot_end)
    # Set tick style
    plt.style.use('seaborn')
    forest = (0.2,0.6,0.2)
    ax2.tick_params(axis="x", direction="out", length=10, width=2, color=forest)
    ax2.tick_params(axis="y", direction="out", length=10, width=2, color=forest)
    ax2.grid(False)

    # Hide tick labels
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels([])
    ax2.set_yticks(tick_pos)
    ax2.set_yticklabels([])
    ax.axis("on")
    ax2.axis("on")
    sns.set_style("ticks")    

    # Draw box around heatmap
    ax2.hlines([0,array_dist.shape[1]], *ax2.get_xlim())
    ax2.vlines([0,array_dist.shape[1]], *ax2.get_ylim())
    
    print('Plot generated.')
    return plt, fig

#%% Run crosscor_dist fucntion
for form in forms:
    plt.close('all')
    plt, fig = crosscor_dist(form)
    fig_path = '6lu7_crosscor_dist_'+form
    fig.savefig(fig_path+'.png', dpi=600, bbox_inches='tight', format='png')
print('Figures saved.')        
