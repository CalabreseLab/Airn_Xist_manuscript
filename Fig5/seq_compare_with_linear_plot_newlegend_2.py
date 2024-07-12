


####### plot transcripts as linear parallel tracks ####################
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np
import pandas as pd

from matplotlib.patches import Polygon

from seekr.fasta_reader import Reader as seekrReader

# divide sequence into chunks around a target length without overlap
# length of seq divided by target length, the quotient is the number of chunks
# if remainder is greater than flex_length, add one more chunk
# if remainder is less than flex_length, the total num of chunk is the quotient
# divide the sequence equally into the total num of chunks, the remiander goes to the last chunk
def divide_seq(seq, target_length, flex_length): 
    seq_length = len(seq)
    quotient = seq_length // target_length
    remainder = seq_length % target_length
    if remainder > flex_length:
        num_chunk = quotient + 1
    else:
        num_chunk = quotient
    chunk_length = seq_length // num_chunk
    chunks = [seq[i:i+chunk_length] for i in range(0, seq_length, chunk_length)]
    # save the index range for each chunk
    # if last index is greater than seq_length, save seq_length instead
    chunk_index = [(i, np.min([i+chunk_length, seq_length])) for i in range(0, seq_length, chunk_length)]
    
    # append the last chunk to the second to the last chunk
    # if the last chunk length is less than chunk_length
    if len(chunks[-1]) < chunk_length:
        chunks[-2] = chunks[-2] + chunks[-1]
        # update chunk_index accordingly
        chunk_index[-2] = (chunk_index[-2][0], chunk_index[-1][1])
        del chunks[-1]
        del chunk_index[-1]
    return chunks, chunk_index

# generate color for each link of mXist based on feature coordinates
# generate color for mXist manual chunks
def Xist_assign_color(header):
    # Define the colors for each feature
    xist_color_dict = {
        'mXist_interval1': '#e377c2',
        'mXist_repeatA': '#d62728',
        'mXist_ss234': '#ff7f0e',
        'mXist_repeatFdwn': '#bcbd22',
        'mXist_repeatB': '#2ca02c',
        'mXist_repeatC': '#17becf',
        'mXist_interval2_0': '#1f77b4',
        'mXist_interval2_1': '#1f77b4',
        'mXist_interval2_2': '#1f77b4',
        'mXist_interval2_3': '#1f77b4',
        'mXist_interval2_4': '#1f77b4',
        'mXist_interval2_5': '#1f77b4',
        'mXist_interval2_6': '#1f77b4',
        'mXist_interval2_7': '#1f77b4',
        'mXist_interval2_8': '#1f77b4',
        'mXist_interval2_9': '#1f77b4',
        'mXist_interval2_10': '#1f77b4',
        'mXist_repeatEbroad': '#9467bd',
        'mXist_interval3_0':'#8c564b',
        'mXist_interval3_1':'#8c564b',
        'mXist_interval3_2':'#8c564b',
        'mXist_interval3_3':'#8c564b',
        'mXist_interval3_4':'#8c564b',
        'mXist_interval3_5':'#8c564b',
        'mXist_interval3_6':'#8c564b',
        'mXist_interval3_7':'#8c564b',
        'mXist_interval3_8':'#8c564b',
        'mXist_interval3_9':'#8c564b',
        'mXist_interval3_10':'#8c564b',
        'mXist_interval3_11':'#8c564b',
        'mXist_interval3_12':'#8c564b'
    }
    
    return xist_color_dict[header]


# generate color for A and K
# use color from matplotlib colormaps tab10
# A, 18 chunks (498) per color block, the last two blocks will have 17 chunks (178 in total)
# K, 19 chunks (500) per color block, the last four blocks will have 18 chunks (186 in total)
# airnseq=seekrReader('airn.fa').get_seqs()[0]
# _,ad=divide_seq(airnseq, 500, 200)
# kseq=seekrReader('kcnq1ot1.fa').get_seqs()[0]
# _,kd=divide_seq(kseq, 500, 200)

airn_color_dict = {
    'a1': '#e377c2',
    'a2': '#d62728',
    'a3': '#ff7f0e',
    'a4': '#bcbd22',
    'a5': '#2ca02c',
    'a6': '#17becf',
    'a7': '#1f77b4',
    'a8': '#9467bd',
    'a9': '#8c564b',
    'a10': '#7f7f7f',
}

airn_coord_dict = {
    'a1':[0,8964],
    'a2':[8964,17928],
    'a3':[17928,26892],
    'a4':[26892,35856],
    'a5':[35856,44820],
    'a6':[44820,53784],
    'a7':[53784,62748],
    'a8':[62748,71712],
    'a9':[71712,80178],
    'a10':[80178,88754],
}

k_color_dict = {
    'k1': '#e377c2',
    'k2': '#d62728',
    'k3': '#ff7f0e',
    'k4': '#bcbd22',
    'k5': '#2ca02c',
    'k6': '#17becf',
    'k7': '#1f77b4',
    'k8': '#9467bd',
    'k9': '#8c564b',
    'k10': '#7f7f7f',
}

k_coord_dict = {
    'k1':[0,9500],
    'k2':[9500,19000],
    'k3':[19000,28500],
    'k4':[28500,38000],
    'k5':[38000,47500],
    'k6':[47500,57000],
    'k7':[57000,66000],
    'k8':[66000,75000],
    'k9':[75000,84000],
    'k10':[84000,93092],
}



def general_assign_color(c_start, c_end, coord_dict, color_dict):

    # check which block the c_start and c_end are in from coord_dict
    for key, value in coord_dict.items():
        if c_start >= value[0] and c_end <= value[1]:
            color = color_dict[key]
            break
 
    return color




def plot_AKtranscripts(seq1path, seq2path, seq1name, seq2name, 
                       seq1lk,seq2lk,
                       target_length, flex_length,
                       k, sim_df_path, q_cutoff):
    
    
    seq1 = seekrReader(seq1path).get_seqs()[0]
    seq2 = seekrReader(seq2path).get_seqs()[0]


    _, seq1_chunk_index = divide_seq(seq1, target_length, flex_length)
    _, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)


    # get the first index for each chunk
    seq1_chunk_start = [index[0] for index in seq1_chunk_index]
    seq2_chunk_start = [index[0] for index in seq2_chunk_index]

    # get the second index for each chunk
    seq1_chunk_end = [index[1] for index in seq1_chunk_index]
    seq2_chunk_end = [index[1] for index in seq2_chunk_index]

    # generate chunk number
    seq1_chunk_num = len(seq1_chunk_start)
    seq2_chunk_num = len(seq2_chunk_start)

    # Calculate the total length for each transcript
    total_length1 = len(seq1)
    total_length2 = len(seq2)
    
    # Calculate the x-coordinates for the lines
    # x1 = np.linspace(-total_length1/2, total_length1/2, total_length1)
    # x2 = np.linspace(-total_length2/2, total_length2/2, total_length2)
    if seq1name == 'Airn':
        coord_dict = airn_coord_dict
        color_dict = airn_color_dict
    elif seq1name == 'Kcnq1ot1':
        coord_dict = k_coord_dict
        color_dict = k_color_dict

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(30, 5))

    # Plot the lines representing the lengths of transcript1 and transcript2
    # ax.plot(x1, np.zeros_like(x1), color='red', linewidth=8, label=seq1name)
    # ax.plot(x2, np.ones_like(x2)*5, color='blue', linewidth=8, label=seq2name)


    # Add transcript names on the top or bottom
    # ax.text(0, 5.2, seq1name, ha='center', va='center', color='black', fontsize=10)
    # ax.text(0, -0.25, seq2name, ha='center', va='center', color='black', fontsize=10)


    # generate color for seq1 chunk with name 'mXist' using Xist_assign_color function
    # if not named 'mXist', use general_assign_color function, color block set to 10000
    seq1_chunk_color = [general_assign_color(c_start, c_end, coord_dict, color_dict) for c_start, c_end in zip(seq1_chunk_start, seq1_chunk_end)]
    

    # generate color for seq2 chunk, get dark grey and light grey and cycle through them
    gmap = ['#d3d3d3','#8c8c8c']
    seq2_chunk_color = [gmap[i % 2] for i in range(seq2_chunk_num)]
    
    # get the longer seq length
    maxlen = max(total_length1, total_length2)
    
    # directly plot chunks without the lines
    # Plot transcript1 chunks as rectangular boxes beneath the line
    for i in range(seq1_chunk_num):
        # scale the chunks to the same length of seq1 and seq2
        start = seq1_chunk_start[i]*maxlen/total_length1
        end = seq1_chunk_end[i]*maxlen/total_length1-1
        # start = seq1_chunk_start[i]
        # end = seq1_chunk_end[i]-1
        color = seq1_chunk_color[i]
        # label = str(i+1)
        # plot the chunk box
        ax.add_patch(plt.Rectangle((start - maxlen/2, 4.8), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # ax.add_patch(plt.Rectangle((start - total_length1/2, 4.8), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        
        # add chunk name in the middle of the box
        # ax.text(start - total_length1/2 + (end-start)/2, 0.1, label, ha='center', va='center', color='white', fontsize=4)
        # add chunk start coordinate as a tick
        # ax.plot([(start - total_length1/2), (start - total_length1/2)], [-0.1, -0.05], color='black', linewidth=0.8)
        # add chunk start coordinate as a label
        # ax.text(start - total_length1/2, -0.15, str(start), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq1_chunk-sum) add the end coordinate as a tick
        # if i == seq1_chunk_num - 1:
        #     ax.plot([end - total_length1/2, end - total_length1/2], [-0.1, -0.05], color='black', linewidth=0.8)
        #     ax.text(end - total_length1/2, -0.15, str(end), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        
    # Plot transcript2 chunks as rectangular boxes above the line
    for i in range(seq2_chunk_num):
        # scale the chunks to the same length of seq1 and seq2
        start = seq2_chunk_start[i]*maxlen/total_length2
        end = seq2_chunk_end[i]*maxlen/total_length2-1
        # start = seq2_chunk_start[i]
        # end = seq2_chunk_end[i]-1
        color = seq2_chunk_color[i]
        # label = str(i+1)
        ax.add_patch(plt.Rectangle((start - maxlen/2, 0), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # ax.add_patch(plt.Rectangle((start - total_length2/2, 0), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))

        # ax.text(start - total_length2/2 + (end-start)/2, 4.9, label, ha='center', va='center', color='white', fontsize=4)
        # ax.plot([start - total_length2/2, start - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        # ax.text(start - total_length2/2, 5.15, str(start), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq2_chunk-sum) add the end coordinate as a tick
        # if i == seq2_chunk_num - 1:
        #     ax.plot([end - total_length2/2, end - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        #     ax.text(end - total_length2/2, 5.15, str(end), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)

    # load in saved csv file
    sim_df_lg = pd.read_csv(sim_df_path)
    # Firstly, filter the DataFrame to keep rows where pval < pval_cutoff
    sig_chunks = sim_df_lg[sim_df_lg[q_cutoff] == True]
    #sig_chunks = sim_df_lg[sim_df_lg['seekr_rscore'] > 0.22]

    # check if sig chunks is not empty
    if sig_chunks.empty:
        plt.show()     
    else:
        # separate seq_chunk in sig_chunks by _
        seq1_link = sig_chunks[seq1lk].str.split('_')
        seq2_link = sig_chunks[seq2lk].str.split('_')
        # save the second and third element as start and end coordinates
        seq1_link_start = [int(chunk[1]) for chunk in seq1_link]
        seq1_link_end = [int(chunk[2]) for chunk in seq1_link]
        seq2_link_start = [int(chunk[1]) for chunk in seq2_link]
        seq2_link_end = [int(chunk[2]) for chunk in seq2_link]
        
        # get seekr_rscore as alpha for links
        # link_alpha = sig_chunks['seekr_rscore'].tolist()

        # normaized seekr_rscore with max seekr_rscore as link alpha
        # link_alpha = np.array(sig_chunks['seekr_rscore']/sig_chunks['seekr_rscore'].max())

        # generate color for each link
        # if seq1name == 'mXist', use Xist_assign_color function
        # else use general_assign_color function, color block set to 10000
        link_color = [general_assign_color(c_start, c_end, coord_dict, color_dict) for c_start, c_end in zip(seq1_link_start, seq1_link_end)]

        # plot links
        for i in range (len(seq1_link)):
            # Plot the link as a ploygon
            # generate points a Nx2 array with seq1_link_start, seq1_link_end as x and 0.3 as y
            # scale the start and end to the same scale of chunks
            points = np.array([[((seq1_link_start[i]*maxlen/total_length1) - maxlen/2), 4.76], [((seq1_link_end[i]*maxlen/total_length1) - maxlen/2 -1), 4.76],
                              [((seq2_link_start[i]*maxlen/total_length2) - maxlen/2), 0.24], [((seq2_link_end[i]*maxlen/total_length2) - maxlen/2 -1), 0.24]])
            # points = np.array([[seq1_link_start[i] - total_length1/2, 4.76], [(seq1_link_end[i] - total_length1/2 -1), 4.76],
            #                   [seq2_link_start[i] - total_length2/2, 0.24], [(seq2_link_end[i] - total_length2/2 -1), 0.24]])
                               
            # polygon = Polygon(points, closed=True, color=link_color[i], alpha=link_alpha[i])
            # set alpha to be unique for all links
            polygon = Polygon(points, closed=True, color=link_color[i], alpha=0.15)
            ax.add_patch(polygon)

    # Set axis limits
    ax.set_xlim(-maxlen/2, maxlen/2)
    
    ax.set_ylim(-2,7)
    
    # Remove ticks and labels from x and y axes
    ax.axis('off')
    
    # Remove the spines (outline)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add legend
    # ax.legend(loc='upper left')

    # save figure as svg file
    plt.savefig(f'{seq1name}_{seq2name}_{target_length}nt_k{k}_{q_cutoff}_seekr_linear_manxist_2.svg', dpi=300, format='svg')

    
    #plt.show()
    

plot_AKtranscripts('airn.fa', 'kcnq1ot1.fa', 'Airn', 'Kcnq1ot1', 
                   seq1lk='seq1_chunk',seq2lk='seq2_chunk',
                   target_length=500, flex_length=200, k=4, 
                   sim_df_path='airn_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                   q_cutoff='q975')

plot_AKtranscripts('kcnq1ot1.fa','airn.fa', 'Kcnq1ot1', 'Airn', 
                   seq1lk='seq2_chunk',seq2lk='seq1_chunk',
                   target_length=500, flex_length=200, k=4, 
                   sim_df_path='airn_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                   q_cutoff='q975')


def plot_Xtranscripts(seq2path, seq2name, target_length, flex_length,
                      k, sim_df_path, q_cutoff):
    
    
    seq2 = seekrReader(seq2path).get_seqs()[0]

    _, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)


    # get the first index for each chunk
    seq1_chunk_start = [0,354,745,1497,2859,3080,4692,
                        5193,5694,6195,6696,7197,7698,8199,
                        8700,9201,9702,10213,11630,12115,
                        12600,13085,13570,14055,14540,15025,
                        15510,15995,16480,16965,17450]
    seq2_chunk_start = [index[0] for index in seq2_chunk_index]

    # get the second index for each chunk
    seq1_chunk_end = [354,745,1497,2859,3080,4692,
                      5193,5694,6195,6696,7197,7698,8199,
                      8700,9201,9702,10213,11630,12115,
                      12600,13085,13570,14055,14540,15025,
                      15510,15995,16480,16965,17450,17946]
    seq2_chunk_end = [index[1] for index in seq2_chunk_index]

    # generate chunk number
    seq1_chunk_num = 31
    seq2_chunk_num = len(seq2_chunk_start)

    # Calculate the total length for each transcript
    total_length1 = 17946
    total_length2 = len(seq2)
    
    # Calculate the x-coordinates for the lines
    # x1 = np.linspace(-total_length1/2, total_length1/2, total_length1)
    # x2 = np.linspace(-total_length2/2, total_length2/2, total_length2)


    # Create figure and axis
    fig, ax = plt.subplots(figsize=(30, 5))

    # Plot the lines representing the lengths of transcript1 and transcript2
    # ax.plot(x1, np.zeros_like(x1), color='red', linewidth=8, label=seq1name)
    # ax.plot(x2, np.ones_like(x2)*5, color='blue', linewidth=8, label=seq2name)

    maxlen = max(total_length1, total_length2)

    # Set axis limits
    ax.set_xlim(-maxlen/2, maxlen/2)
    
    ax.set_ylim(-2,7)
    

    # Add transcript names on the top or bottom
    # ax.text(0, 5.2, 'Xist', ha='center', va='center', color='black', fontsize=10)
    # ax.text(0, -0.25, seq2name, ha='center', va='center', color='black', fontsize=10)

    xist_color_dict = {
        'mXist_interval1': '#e377c2',
        'mXist_repeatA': '#d62728',
        'mXist_ss234': '#ff7f0e',
        'mXist_repeatFdwn': '#bcbd22',
        'mXist_repeatB': '#2ca02c',
        'mXist_repeatC': '#17becf',
        'mXist_interval2_0': '#1f77b4',
        'mXist_interval2_1': '#1f77b4',
        'mXist_interval2_2': '#1f77b4',
        'mXist_interval2_3': '#1f77b4',
        'mXist_interval2_4': '#1f77b4',
        'mXist_interval2_5': '#1f77b4',
        'mXist_interval2_6': '#1f77b4',
        'mXist_interval2_7': '#1f77b4',
        'mXist_interval2_8': '#1f77b4',
        'mXist_interval2_9': '#1f77b4',
        'mXist_interval2_10': '#1f77b4',
        'mXist_repeatEbroad': '#9467bd',
        'mXist_interval3_0':'#8c564b',
        'mXist_interval3_1':'#8c564b',
        'mXist_interval3_2':'#8c564b',
        'mXist_interval3_3':'#8c564b',
        'mXist_interval3_4':'#8c564b',
        'mXist_interval3_5':'#8c564b',
        'mXist_interval3_6':'#8c564b',
        'mXist_interval3_7':'#8c564b',
        'mXist_interval3_8':'#8c564b',
        'mXist_interval3_9':'#8c564b',
        'mXist_interval3_10':'#8c564b',
        'mXist_interval3_11':'#8c564b',
        'mXist_interval3_12':'#8c564b'
    }


    # generate color for seq1 chunk with name 'mXist' using Xist_assign_color function
    # if not named 'mXist', use general_assign_color function, color block set to 10000
    seq1_chunk_color = list(xist_color_dict.values())
    # xistchunkname = ['interval1','repeatA','ss234','repeatFdwn',
    #                  'repeatB','repeatC','interval2_0','interval2_1',
    #                  'interval2_2', 'interval2_3','interval2_4','interval2_5',
    #                  'interval2_6','interval2_7','interval2_8','interval2_9',
    #                  'interval2_10','repeatEbroad','interval3_0','interval3_1',
    #                  'interval3_2','interval3_3','interval3_4','interval3_5',
    #                  'interval3_6','interval3_7','interval3_8','interval3_9',
    #                  'interval3_10','interval3_11','interval3_12']


    # generate color for seq2 chunk, get dark grey and light grey and cycle through them
    gmap = ['#d3d3d3','#8c8c8c']
    seq2_chunk_color = [gmap[i % 2] for i in range(seq2_chunk_num)]
    
    
    # directly plot chunks without the lines
    # Plot transcript1 chunks as rectangular boxes beneath the line
    for i in range(seq1_chunk_num):
        start = seq1_chunk_start[i]*maxlen/total_length1
        end = (seq1_chunk_end[i]*maxlen/total_length1)-1
        color = seq1_chunk_color[i]
        # label = str(i+1)
        # plot the chunk box
        ax.add_patch(plt.Rectangle((start - maxlen/2, 4.8), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # add chunk name below the box
        # ax.text(start - total_length1/2 + (end-start)/2, 0.1, label,
        #         va='center', ha='center',
        #         color='white', fontsize=4)
        # add chunk start coordinate as a tick
        # ax.plot([(start - total_length1/2), (start - total_length1/2)], [-0.1, -0.05], color='black', linewidth=0.8)
        # add chunk start coordinate as a label
        # ax.text(start - total_length1/2, -0.15, str(start), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq1_chunk-sum) add the end coordinate as a tick
        # if i == seq1_chunk_num - 1:
        #     ax.plot([end - total_length1/2, end - total_length1/2], [-0.1, -0.05], color='black', linewidth=0.8)
        #     ax.text(end - total_length1/2, -0.15, str(end), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        
    # Plot transcript2 chunks as rectangular boxes above the line
    for i in range(seq2_chunk_num):
        start = seq2_chunk_start[i]*maxlen/total_length2
        end = (seq2_chunk_end[i]*maxlen/total_length2)-1
        color = seq2_chunk_color[i]
        # label = str(i+1)
        ax.add_patch(plt.Rectangle((start - maxlen/2, 0), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # ax.text(start - total_length2/2 + (end-start)/2, 4.9, label, ha='center', va='center', color='white', fontsize=4)
        # ax.plot([start - total_length2/2, start - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        # ax.text(start - total_length2/2, 5.15, str(start), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq2_chunk-sum) add the end coordinate as a tick
        # if i == seq2_chunk_num - 1:
        #     ax.plot([end - total_length2/2, end - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        #     ax.text(end - total_length2/2, 5.15, str(end), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)

    # load in saved csv file
    sim_df_lg = pd.read_csv(sim_df_path)
    sig_chunks = sim_df_lg[sim_df_lg[q_cutoff] == True]
    #sig_chunks = sim_df_lg[sim_df_lg['seekr_rscore'] > 0.22]


    xistchunks = {
        'mXist_interval1': [0,354],
        'mXist_repeatA': [354,745],
        'mXist_ss234': [745,1497],
        'mXist_repeatFdwn': [1497,2859],
        'mXist_repeatB': [2859,3080],
        'mXist_repeatC': [3080,4692],
        'mXist_interval2_0': [4692,5193],
        'mXist_interval2_1': [5193,5694],
        'mXist_interval2_2': [5694,6195],
        'mXist_interval2_3': [6195,6696],
        'mXist_interval2_4': [6696,7197],
        'mXist_interval2_5': [7197,7698],
        'mXist_interval2_6': [7698,8199],
        'mXist_interval2_7': [8199,8700],
        'mXist_interval2_8': [8700,9201],
        'mXist_interval2_9': [9201,9702],
        'mXist_interval2_10': [9702,10213],
        'mXist_repeatEbroad': [10213,11630],
        'mXist_interval3_0':[11630,12115],
        'mXist_interval3_1':[12115,12600],
        'mXist_interval3_2':[12600,13085],
        'mXist_interval3_3':[13085,13570],
        'mXist_interval3_4':[13570,14055],
        'mXist_interval3_5':[14055,14540],
        'mXist_interval3_6':[14540,15025],
        'mXist_interval3_7':[15025,15510],
        'mXist_interval3_8':[15510,15995],
        'mXist_interval3_9':[15995,16480],
        'mXist_interval3_10':[16480,16965],
        'mXist_interval3_11':[16965,17450],
        'mXist_interval3_12':[17450,17946]
    }

    

    # check if sig chunks is not empty
    if sig_chunks.empty:
        plt.show()     
    else:
        # separate seq_chunk in sig_chunks by _
        seq1_link = sig_chunks['seq1_chunk']
        seq2_link = sig_chunks['seq2_chunk'].str.split('_')
        # save the second and third element as start and end coordinates
        seq1_link_start = [xistchunks[chunk][0] for chunk in seq1_link]
        seq1_link_end = [xistchunks[chunk][1] for chunk in seq1_link]
        seq2_link_start = [int(chunk[1]) for chunk in seq2_link]
        seq2_link_end = [int(chunk[2]) for chunk in seq2_link]
        
        # get seekr_rscore as alpha for links
        # link_alpha = sig_chunks['seekr_rscore'].tolist()

        # normaized seekr_rscore with max seekr_rscore as link alpha
        # link_alpha = np.array(sig_chunks['seekr_rscore']/sig_chunks['seekr_rscore'].max())

        # generate color for each link
        # if seq1name == 'mXist', use Xist_assign_color function
        # else use general_assign_color function, color block set to 10000
        link_color = [xist_color_dict[header] for header in seq1_link]
        
        # plot links
        for i in range (len(seq1_link)):
            # Plot the link as a ploygon
            # generate points a Nx2 array with seq1_link_start, seq1_link_end as x and 0.3 as y
            points = np.array([[(seq1_link_start[i]*maxlen/total_length1 - maxlen/2), 4.76], [(seq1_link_end[i]*maxlen/total_length1 - maxlen/2 -1), 4.76],
                              [(seq2_link_start[i]*maxlen/total_length2 - maxlen/2), 0.24], [(seq2_link_end[i]*maxlen/total_length2 - maxlen/2 -1), 0.24]])
                               
            #polygon = Polygon(points, closed=True, color=link_color[i], alpha=link_alpha[i])
            polygon = Polygon(points, closed=True, color=link_color[i], alpha=0.2)
            ax.add_patch(polygon)

    
    # Remove ticks and labels from x and y axes
    ax.axis('off')
    
    # Remove the spines (outline)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add legend
    # ax.legend(loc='upper left')

    # save figure as pdf file
    plt.savefig(f'Xist_{seq2name}_{target_length}nt_k{k}_{q_cutoff}_seekr_linear_manxist_2.svg', dpi=300, format='svg')

    
    #plt.show()



plot_Xtranscripts('airn.fa', 'Airn', 
                  target_length=500, flex_length=200,  q_cutoff='q975',
                  k=4, sim_df_path='mXist_airn_500nt_k4_seekr_rscore_qlist_manxist.csv')


plot_Xtranscripts('kcnq1ot1.fa', 'Kcnq1ot1',
                  target_length=500, flex_length=200, q_cutoff='q975',
                  k=4, sim_df_path='mXist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv')



def plot_againstXtranscripts(seq2path, seq2name, target_length, flex_length,
                             k, sim_df_path, q_cutoff):
    
    
    seq2 = seekrReader(seq2path).get_seqs()[0]

    _, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)


    # get the first index for each chunk
    seq1_chunk_start = [0,354,745,1497,2859,3080,4692,
                        5193,5694,6195,6696,7197,7698,8199,
                        8700,9201,9702,10213,11630,12115,
                        12600,13085,13570,14055,14540,15025,
                        15510,15995,16480,16965,17450]
    seq2_chunk_start = [index[0] for index in seq2_chunk_index]

    # get the second index for each chunk
    seq1_chunk_end = [354,745,1497,2859,3080,4692,
                      5193,5694,6195,6696,7197,7698,8199,
                      8700,9201,9702,10213,11630,12115,
                      12600,13085,13570,14055,14540,15025,
                      15510,15995,16480,16965,17450,17946]
    seq2_chunk_end = [index[1] for index in seq2_chunk_index]

    # generate chunk number
    seq1_chunk_num = 31
    seq2_chunk_num = len(seq2_chunk_start)

    # Calculate the total length for each transcript
    total_length1 = 17946
    total_length2 = len(seq2)

    if seq2name == 'Airn':
        coord_dict = airn_coord_dict
        color_dict = airn_color_dict
    elif seq2name == 'Kcnq1ot1':
        coord_dict = k_coord_dict
        color_dict = k_color_dict
    
    # Calculate the x-coordinates for the lines
    # x1 = np.linspace(-total_length1/2, total_length1/2, total_length1)
    # x2 = np.linspace(-total_length2/2, total_length2/2, total_length2)


    # Create figure and axis
    fig, ax = plt.subplots(figsize=(30, 5))

    # Plot the lines representing the lengths of transcript1 and transcript2
    # ax.plot(x1, np.zeros_like(x1), color='red', linewidth=8, label=seq1name)
    # ax.plot(x2, np.ones_like(x2)*5, color='blue', linewidth=8, label=seq2name)
    maxlen = max(total_length1, total_length2)

    # Set axis limits
    ax.set_xlim(-maxlen/2, maxlen/2)
    
    ax.set_ylim(-2,7)
    

    # Add transcript names on the top or bottom
    # ax.text(0, 5.2, seq2name, ha='center', va='center', color='black', fontsize=10)
    # ax.text(0, -0.25, 'Xist', ha='center', va='center', color='black', fontsize=10)



    # generate color for seq2 chunk use general_assign_color function, color block set to 10000
    seq2_chunk_color = [general_assign_color(c_start, c_end, coord_dict, color_dict) for c_start, c_end in zip(seq2_chunk_start, seq2_chunk_end)]
    

    # generate color for mXist as dark grey and light grey and cycle through them
    gmap = ['#d3d3d3','#8c8c8c']
    seq1_chunk_color = [gmap[i % 2] for i in range(seq1_chunk_num)]
    
    
    # directly plot chunks without the lines
    # Plot transcript1 chunks as rectangular boxes beneath the line
    for i in range(seq1_chunk_num):
        start = seq1_chunk_start[i]*maxlen/total_length1
        end = seq1_chunk_end[i]*maxlen/total_length1-1
        color = seq1_chunk_color[i]
        # label = str(i+1)
        # plot the chunk box
        ax.add_patch(plt.Rectangle((start - maxlen/2, 0), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # add chunk name below the box
        # ax.text(start - total_length1/2 + (end-start)/2, 4.9, label,
        #         va='center', ha='center',
        #         color='white', fontsize=4)
        # add chunk start coordinate as a tick
        # ax.plot([(start - total_length1/2), (start - total_length1/2)], [-0.1, -0.05], color='black', linewidth=0.8)
        # add chunk start coordinate as a label
        # ax.text(start - total_length1/2, -0.15, str(start), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq1_chunk-sum) add the end coordinate as a tick
        # if i == seq1_chunk_num - 1:
        #     ax.plot([end - total_length1/2, end - total_length1/2], [-0.1, -0.05], color='black', linewidth=0.8)
        #     ax.text(end - total_length1/2, -0.15, str(end), ha='center', va='top', color='black', fontsize=5, weight='bold', rotation=90)
        
    # Plot transcript2 chunks as rectangular boxes above the line
    for i in range(seq2_chunk_num):
        start = seq2_chunk_start[i]*maxlen/total_length2
        end = seq2_chunk_end[i]*maxlen/total_length2-1
        color = seq2_chunk_color[i]
        # label = str(i+1)
        ax.add_patch(plt.Rectangle((start - maxlen/2, 4.8), (end-start), 0.2, facecolor=color, edgecolor='black', linewidth=0.5))
        # ax.text(start - total_length2/2 + (end-start)/2, 0.1, label, ha='center', va='center', color='white', fontsize=4)
        # ax.plot([start - total_length2/2, start - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        # ax.text(start - total_length2/2, 5.15, str(start), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)
        # for the last i in range(seq2_chunk-sum) add the end coordinate as a tick
        # if i == seq2_chunk_num - 1:
        #     ax.plot([end - total_length2/2, end - total_length2/2], [5.05, 5.1], color='black', linewidth=0.8)
        #     ax.text(end - total_length2/2, 5.15, str(end), ha='center', va='bottom', color='black', fontsize=5, weight='bold', rotation=90)

    # load in saved csv file
    sim_df_lg = pd.read_csv(sim_df_path)
    sig_chunks = sim_df_lg[sim_df_lg[q_cutoff] == True]
    #sig_chunks = sim_df_lg[sim_df_lg['seekr_rscore'] > 0.22]


    xistchunks = {
        'mXist_interval1': [0,354],
        'mXist_repeatA': [354,745],
        'mXist_ss234': [745,1497],
        'mXist_repeatFdwn': [1497,2859],
        'mXist_repeatB': [2859,3080],
        'mXist_repeatC': [3080,4692],
        'mXist_interval2_0': [4692,5193],
        'mXist_interval2_1': [5193,5694],
        'mXist_interval2_2': [5694,6195],
        'mXist_interval2_3': [6195,6696],
        'mXist_interval2_4': [6696,7197],
        'mXist_interval2_5': [7197,7698],
        'mXist_interval2_6': [7698,8199],
        'mXist_interval2_7': [8199,8700],
        'mXist_interval2_8': [8700,9201],
        'mXist_interval2_9': [9201,9702],
        'mXist_interval2_10': [9702,10213],
        'mXist_repeatEbroad': [10213,11630],
        'mXist_interval3_0':[11630,12115],
        'mXist_interval3_1':[12115,12600],
        'mXist_interval3_2':[12600,13085],
        'mXist_interval3_3':[13085,13570],
        'mXist_interval3_4':[13570,14055],
        'mXist_interval3_5':[14055,14540],
        'mXist_interval3_6':[14540,15025],
        'mXist_interval3_7':[15025,15510],
        'mXist_interval3_8':[15510,15995],
        'mXist_interval3_9':[15995,16480],
        'mXist_interval3_10':[16480,16965],
        'mXist_interval3_11':[16965,17450],
        'mXist_interval3_12':[17450,17946]
    }

    

    # check if sig chunks is not empty
    if sig_chunks.empty:
        plt.show()     
    else:
        # separate seq_chunk in sig_chunks by _
        seq1_link = sig_chunks['seq1_chunk']
        seq2_link = sig_chunks['seq2_chunk'].str.split('_')
        # save the second and third element as start and end coordinates
        seq1_link_start = [xistchunks[chunk][0] for chunk in seq1_link]
        seq1_link_end = [xistchunks[chunk][1] for chunk in seq1_link]
        seq2_link_start = [int(chunk[1]) for chunk in seq2_link]
        seq2_link_end = [int(chunk[2]) for chunk in seq2_link]
        
        # get seekr_rscore as alpha for links
        # link_alpha = sig_chunks['seekr_rscore'].tolist()

        # normaized seekr_rscore with max seekr_rscore as link alpha
        # link_alpha = np.array(sig_chunks['seekr_rscore']/sig_chunks['seekr_rscore'].max())

        # generate color for each link
        # if seq1name == 'mXist', use Xist_assign_color function
        # else use general_assign_color function, color block set to 10000
        link_color = [general_assign_color(c_start, c_end, coord_dict, color_dict) for c_start, c_end in zip(seq2_link_start, seq2_link_end)]

        
        # plot links
        for i in range (len(seq2_link)):
            # Plot the link as a ploygon
            # generate points a Nx2 array with seq1_link_start, seq1_link_end as x and 0.3 as y
            points = np.array([[(seq1_link_start[i]*maxlen/total_length1 - maxlen/2), 0.24], [(seq1_link_end[i]*maxlen/total_length1 - maxlen/2 -1), 0.24],
                              [(seq2_link_start[i]*maxlen/total_length2 - maxlen/2), 4.76], [(seq2_link_end[i]*maxlen/total_length2 - maxlen/2 -1), 4.76]])
                               
            polygon = Polygon(points, closed=True, color=link_color[i], alpha=0.2)
            # polygon = Polygon(points, closed=True, color=link_color[i], alpha=link_alpha[i])
            ax.add_patch(polygon)

    
    # Remove ticks and labels from x and y axes
    ax.axis('off')
    
    # Remove the spines (outline)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add legend
    # ax.legend(loc='upper left')

    # save figure as pdf file
    plt.savefig(f'{seq2name}_Xist_{target_length}nt_k{k}_{q_cutoff}_seekr_linear_manxist_2.svg', dpi=300, format='svg')

    
    #plt.show()

plot_againstXtranscripts('airn.fa', 'Airn', 
                         target_length=500, flex_length=200,  q_cutoff='q975',
                         k=4, sim_df_path='mXist_airn_500nt_k4_seekr_rscore_qlist_manxist.csv')


plot_againstXtranscripts('kcnq1ot1.fa', 'Kcnq1ot1',
                         target_length=500, flex_length=200, q_cutoff='q975',
                         k=4, sim_df_path='mXist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv')
