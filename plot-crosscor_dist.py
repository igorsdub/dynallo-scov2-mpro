import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *
from numpy import ma
import seaborn as sns

 
def crosscor_dist(form):
    # I/O files
    infile1 = 'dist-'+str(form)+'.dat'
    infile2 = 'crosscor-'+str(form)+'.dat'
    
    aa_num = 612
    het_num = 21
    
    xlbl='Amino Acid Number'
    ylbl='Amino Acid Number'
    mi=[]
    mi2=[]
    mj=[]
    mj2=[]
    ol=[]
    ol2=[]
    i=-1
    j=-1
    		
    inlines=open(infile1,'r').readlines()
    
    if inlines[-1]=='\n':
    	inlines[-1:]=[]
    
    i=i+1
    mi.append([])
    mj.append([])
    ol.append([])
    
    for line in inlines:
    	if line=='\n':
    		i=i+1
    		mi.append([])
    		mj.append([])
    		ol.append([])
    		
    	else:
    		mi[i].append(int(line.split()[0]))
    		mj[i].append(int(line.split()[1]))
    		ol[i].append(float(line.split()[2]))
    				
    mi=np.array(mi)
    mj=np.array(mj)
    ol=np.array(ol)
    
    #------------------------------------------------------------------
    
    inlines=open(infile2,'r').readlines()
    
    if inlines[-1]=='\n':
    	inlines[-1:]=[]
    
    j=j+1
    mi2.append([])
    mj2.append([])
    ol2.append([])
    
    for line in inlines:
    	if line=='\n':
    		j=j+1
    		mi2.append([])
    		mj2.append([])
    		ol2.append([])
    		
    	else:
    		mi2[j].append(int(line.split()[0]))
    		mj2[j].append(int(line.split()[1]))
    		ol2[j].append(float(line.split()[2]))
    		
    mi2=np.array(mi2)
    mj2=np.array(mj2)
    ol2=np.array(ol2)
    
    
    sns.set()
    
    # Create masking matrix
    mask_dist = np.triu(ol)
    mask_cross = np.tril(ol2)
    
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
    ax1 = sns.heatmap(ol, vmin=0,vmax=16, mask=mask_dist, square=True,cmap=plt.cm.gist_yarg_r, 
                    cbar_kws={'label':'Distance [ $\AA{}$ ]','ticks':np.arange(0,18,2)})
    # Plot cross-correlation map 
    ax2 = sns.heatmap(ol2,  vmin=-1,vmax=1, mask=mask_cross, square=True,cmap=plt.cm.RdBu_r, 
                    cbar_kws={'label': 'Cross-correlation','ticks':np.arange(-1.00,1.25,0.25)})
    
    # Alternative brutal method to draw spines
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    
    ax.set_xlabel(xlbl)
    
    ax.set_ylabel(ylbl)
    
    gca()
    
    for item in range(mi.min(), mi.max()):
        if item not in mi:
           gca().add_patch(Rectangle((item,mj.min()),1,mj.max()-mj.min(),color='black'))
           gca().add_patch(Rectangle((mi.min(),item),mi.max()-mi.min(),1,color='black'))
    
    # Invert y-axis
    plt.gca().invert_yaxis()
    plt.xticks(rotation=0)
    
    positions = (99, 199, 299, 407, 507, 607)
    labels = ("100", "200", "300", "100", "200", "300")
    ax.set(yticks=positions, yticklabels=labels)
    ax.set(xticks=positions, xticklabels=labels)
    
    if form == 0:
        line_pos = [0,ol.shape[1], aa_num/2]
        plot_end = aa_num
    elif form == 1:
        line_pos = [0,ol.shape[1],aa_num/2, aa_num]
        plot_end = aa_num + het_num
    else:
        line_pos = [0,ol.shape[1],aa_num/2, aa_num, aa_num + het_num]
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
    ax2.hlines([0,ol.shape[1]], *ax2.get_xlim())
    ax2.vlines([0,ol.shape[1]], *ax2.get_ylim())
    
    print('Plot generated.')
    return plt, fig


#%% Run crosscor_dist fucntion
for f in [0,1,2]:
    plt.close('all')
    plt, fig = crosscor_dist(m,f)
    fig_path = 'graph-'+str(m)+'-'+str(f)
    fig.savefig(fig_path+'.png', dpi=600, bbox_inches='tight', format='png')
print('Figures saved.')        