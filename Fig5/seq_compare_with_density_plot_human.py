

#######################################################################
# generate density plots for each comparison with the full background
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import pickle


from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

# generate background from the .fa file


gc_counter = seekrBasicCounter('gencodevM25_unique_divide_namedchunks_500_withfullairn.fa',
                               mean='gencodevM25_unique_mean_4mers.npy',
                               std='gencodevM25_unique_std_4mers.npy',
                               k=4,silent=True)
gc_counter.make_count_file()
gc_sim = seekrPearson(gc_counter.counts,gc_counter.counts)

# get the values in upper triangle of gc_sim and saved as 1D np array
gc_sim_triu = gc_sim[np.triu_indices(gc_sim.shape[0], k=1)]
# 895808628 values in gc_sim_triu
print(len(gc_sim_triu))

# check how many values in gc_sim_triu are 0
print(len(gc_sim_triu[gc_sim_triu==0]))
# 0

# save gc_sim_triu as npy file
np.save('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn.npy', gc_sim_triu)


def plot_density(bkgpath, comppath, compname):
    
    bkgdist=np.load(bkgpath)

    qlist=pickle.load(open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'rb'))

    q975=qlist['q975']

    # load comppath from csv file
    compdist=pd.read_csv(comppath)
    # get the seekr_rscore column as np array
    compdist=compdist['seekr_rscore'].to_numpy()

    # perform ks test to compare bkgdist and compdist
    # only test for whether compdist is greater than bkgdist
    _, pval = ks_2samp(bkgdist, compdist, alternative='greater')

    # use Airal for fonts
    # Specify the path to your .ttf font file
    font_path = 'arial.ttf'
    # Register the font with matplotlib
    font_manager.fontManager.addfont(font_path)
    # Set the font properties
    prop = font_manager.FontProperties(fname=font_path)

    # plot bkgdist and compdist as density plots in one figure
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['pdf.fonttype'] = 42
    sns.kdeplot(bkgdist, fill=True, 
                label='Background', color='blue', alpha=0.3,
                common_norm=False)
    sns.kdeplot(compdist, fill=True, 
                label=compname, color='red', alpha=0.3,
                common_norm=False)
    plt.xlim(-0.6,0.6)
    plt.ylim(0,3.5)
    # add a solide vertical line at q975
    plt.axvline(x=q975, color='green', linestyle='solid', linewidth=2)
    # add text '97.5%: q975' to the vertical line
    # between 97.5% and the q975 value add a return line
    # align to left upper corner of the text
    plt.text(q975+0.01, 3.4, f'97.5%:\n{q975:.4f}', fontsize=24, ha='left', va='top')
    # add a 35% transparent white rectangle to the left of the vertical line
    # to cover the density plot
    plt.axvspan(-0.6, q975, alpha=0.65, color='white')
    plt.xlabel('Seekr r score', fontsize=36)
    plt.xticks(fontsize=34, rotation=45)
    plt.yticks(fontsize=34)
    plt.ylabel('Density', fontsize=36)
    # plt.title(f'Density Plot of Seekr rscore for {compname} vs background distribution',
    #           fontweight='bold', fontsize=28)
    # set legend fontsize to 28 and fontweight to bold
    # plt.legend(prop={'size':20},loc='upper right')
    # print the pval outside the plot, below the legend
    plt.text(-0.55, 3.4, f'pval: {pval:.2e}', fontsize=24) 
    # hide the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    plt.tight_layout()
    # save figure as pdf file
    # get comppath without extension
    comppath=comppath.split('.')[0]
    plt.savefig(f'{comppath}_bkg_density.pdf', dpi=300, format='pdf')


bkgpath='gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn.npy'


plot_density(bkgpath, comppath='Xist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv', 
             compname='XIST vs KCNQ1OT1')

