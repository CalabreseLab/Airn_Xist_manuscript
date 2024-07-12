
# plot one transcript chunks as x axis
# plot sig chunk pairs as color coded histogram, stacked
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import matplotlib.font_manager as font_manager
import matplotlib as mpl

# use Airal for fonts
# Specify the path to your .ttf font file
font_path = 'arial.ttf'
# Register the font with matplotlib
font_manager.fontManager.addfont(font_path)
# Set the font properties
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = prop.get_name()
#plt.rcParams["axes.labelsize"] = 30

# make sure the plot font in pdf can be editted in using illustrator
mpl.rcParams['pdf.fonttype'] = 42

# # divide sequence into chunks around a target length without overlap
# # length of seq divided by target length, the quotient is the number of chunks
# # if remainder is greater than flex_length, add one more chunk
# # if remainder is less than flex_length, the total num of chunk is the quotient
# # divide the sequence equally into the total num of chunks, the remiander goes to the last chunk
# def divide_seq(seq, target_length, flex_length): 
#     seq_length = len(seq)
#     quotient = seq_length // target_length
#     remainder = seq_length % target_length
#     if remainder > flex_length:
#         num_chunk = quotient + 1
#     else:
#         num_chunk = quotient
#     chunk_length = seq_length // num_chunk
#     chunks = [seq[i:i+chunk_length] for i in range(0, seq_length, chunk_length)]
#     # save the index range for each chunk
#     # if last index is greater than seq_length, save seq_length instead
#     chunk_index = [(i, np.min([i+chunk_length, seq_length])) for i in range(0, seq_length, chunk_length)]
    
#     # append the last chunk to the second to the last chunk
#     # if the last chunk length is less than chunk_length
#     if len(chunks[-1]) < chunk_length:
#         chunks[-2] = chunks[-2] + chunks[-1]
#         # update chunk_index accordingly
#         chunk_index[-2] = (chunk_index[-2][0], chunk_index[-1][1])
#         del chunks[-1]
#         del chunk_index[-1]
#     return chunks, chunk_index

######## generate xak_sig_list csv for sig chunk pairs that are shared ########
# check if the chunk has sig pair for both A and K
# if shared both, both column is 1 otherwise 0
def parse_string(s):
    # check if the string is NaN
    if pd.isna(s):
        return 0

    # split the string by ';' and then by '_' 
    parts = [part.split('_')[0] for part in s.split(';')]

    # get the unique parts
    unique_parts = set(parts)

    # if there are 2 unique parts, return 1, otherwise return 0
    return 1 if len(unique_parts) == 2 else 0

def sig_list_generator(xakchunk_path, header):
    xakdf=pd.read_csv(xakchunk_path)
    # drop last two columns
    xakdf=xakdf.iloc[:, :-2]
    # rename columns
    xakdf.columns=['chunk', 'sig_chunk']
    # replace header with '' in chunk column
    xakdf.loc[:, 'chunk'] = xakdf['chunk'].str.replace(header, '', regex=True)
    # add xpos
    xakdf['xpos'] = range(len(xakdf))
    # create the 'both' column by applying the parse_string function to the 'sig_chunk' column
    xakdf['both'] = xakdf['sig_chunk'].apply(parse_string)
    # drop sig_chunk column
    xakdf=xakdf.drop(columns=['sig_chunk'])
    # save to csv
    xakdf.to_csv(f'{header}sig_list.csv', index=False)

sig_list_generator('mXist_XAK_q975chunk.csv',header='mXist_')
sig_list_generator('airn_XAK_q975chunk.csv',header='airn_')
sig_list_generator('kcnq1ot1_XAK_q975chunk.csv',header='kcnq1ot1_')


# mark each chunk of X that are shared in both A and K
# and those only with A or K
# and those are not shared at all

def count_string(s):
    # check if the string is NaN
    if pd.isna(s):
        return 0

    # split the string by ';' and then by '_' 
    parts = [part.split('_')[0] for part in s.split(';')]

    # get the unique parts
    unique_parts = set(parts)

    # convert unique_parts to string
    unique_parts = ';'.join(unique_parts)

    # return the unique rna names
    return unique_parts

def sig_count_generator(xakchunk_path, header):
    xakdf=pd.read_csv(xakchunk_path)
    # drop last two columns
    xakdf=xakdf.iloc[:, :-2]
    # rename columns
    xakdf.columns=['chunk', 'sig_chunk']
    # replace header with '' in chunk column
    # xakdf.loc[:, 'chunk'] = xakdf['chunk'].str.replace(header, '', regex=True)
    # create the 'sig_pair' column by applying the parse_string function to the 'sig_chunk' column
    xakdf['sig_pair'] = xakdf['sig_chunk'].apply(count_string)
    # drop sig_chunk column
    xakdf=xakdf.drop(columns=['sig_chunk'])
    # save to csv
    xakdf.to_csv(f'{header}sig_counts.csv', index=False)
    
sig_count_generator('mXist_XAK_q975chunk.csv',header='mXist_')
sig_count_generator('airn_XAK_q975chunk.csv',header='airn_')
sig_count_generator('kcnq1ot1_XAK_q975chunk.csv',header='kcnq1ot1_')






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

color_order = ['#e377c2', '#d62728', '#ff7f0e', '#bcbd22', '#2ca02c', 
               '#17becf', '#1f77b4', '#9467bd', '#8c564b', '#7f7f7f']

def general_assign_color(c_start, c_end, coord_dict, color_dict):

    # check which block the c_start and c_end are in from coord_dict
    for key, value in coord_dict.items():
        if c_start >= value[0] and c_end <= value[1]:
            color = color_dict[key]
            break
 
    return color


def plot_transcripts(seq1name, chunkcol, refcol,
                     sim_df_path,sig_list_path, q_cutoff,ph,addht):
    
    
    # seq1 = seekrReader(seq1path).get_seqs()[0]
    # total_length1 = len(seq1)

    # load in saved csv file
    sim_df_lg = pd.read_csv(sim_df_path)
    # sig_chunks = sim_df_lg[sim_df_lg['pval'] < pval_cutoff]

    # add a 'color' column to sim_df_lg
    sim_df_lg['chunk_color'] = ''

    # parse seq2_chunk column to get seq2name
    seq2name = sim_df_lg[refcol][0].split('_')[0]

    # if seq2name is mXist, change it to Xist
    # else if seq2name is airn, change it to Airn
    if seq2name == 'mXist':
        seq2name = 'Xist'
    elif seq2name == 'airn':
        seq2name = 'Airn'
    elif seq2name == 'kcnq1ot1':
        seq2name = 'Kcnq1ot1'

    if seq1name == 'Airn':
        color_dict = airn_color_dict
        coord_dict = airn_coord_dict
    elif seq1name == 'Kcnq1ot1':
        color_dict = k_color_dict
        coord_dict = k_coord_dict

    

    # get rid of the first word after split with _ from seq2_chunk column
    sim_df_lg[refcol] = sim_df_lg[refcol].apply(lambda x: '_'.join(x.split('_')[1:]))
    
    # generate color for seq1 chunk with name 'mXist' using Xist_assign_color function
    # if not named 'mXist', use general_assign_color function
    if seq1name == 'Xist':

        # assign color to the row based on Xist_assign_color function
        for index, row in sim_df_lg.iterrows():
            sim_df_lg.loc[index, 'chunk_color'] = Xist_assign_color(row[chunkcol])

    else:
        # color_feature=general_color_feature
        # for each row of sim_df_lg
        # parse c_start and c_end from seq1_chunk column
        # assign color to the row based on general_assign_color function
        for index, row in sim_df_lg.iterrows():
            c_start = int(row[chunkcol].split('_')[1])
            c_end = int(row[chunkcol].split('_')[2])
            sim_df_lg.loc[index, 'chunk_color'] = general_assign_color(c_start, c_end, coord_dict, color_dict)


    # Assuming your DataFrame is named df.
    # Firstly, filter the DataFrame to keep rows where pval < pval_cutoff
    filtered_df = sim_df_lg[sim_df_lg[q_cutoff] == True]

    # Group by 'seq2_chunk' and 'color', then count the number of rows in each group
    grouped_df = filtered_df.groupby([refcol, 'chunk_color']).size().reset_index(name='counts')

    # Pivot the DataFrame so that 'seq2_chunk' is the index, 'color' columns, and 'counts' are the values
    pivoted_df = grouped_df.pivot(index=refcol, columns='chunk_color', values='counts')

    # Reindex the pivoted DataFrame with the unique 'seq2_chunk' values from the original DataFrame
    unique_seq2_chunks = sim_df_lg[refcol].unique()
    pivoted_df = pivoted_df.reindex(unique_seq2_chunks)

    # Sort the DataFrame by 'seq2_chunk'
    # pivoted_df = pivoted_df.sort_index()

    # Fill any NaN values with 0
    pivoted_df = pivoted_df.fillna(0)

    # Filter columns from color_order that exist in pivoted_df
    filtered_colors = [color for color in color_order if color in pivoted_df.columns]

    # Reindex pivoted_df based on the filtered_colors
    pivoted_df = pivoted_df[filtered_colors]


    sig_list = pd.read_csv(sig_list_path)
    # convert to intege
    sig_list['both'] = sig_list['both'].astype(int)
    sig_list['barht'] = sig_list['chunk'].apply(lambda x: pivoted_df.loc[x].sum() if x in pivoted_df.index else 0)

    # Create an identifier for each group of consecutive 'T's
    sig_list['group'] = sig_list['both'].diff().ne(0).cumsum()

    # Calculate group sizes
    group_sizes = sig_list.groupby('group')['both'].transform('size')

    # Assign True to the new 'hotspot' column where group size is >=5
    sig_list['hotspot'] = (group_sizes >= 3) & (sig_list['both'] == 1)


    # use Airal for fonts
    # Specify the path to your .ttf font file
    font_path = 'arial.ttf'
    # Register the font with matplotlib
    font_manager.fontManager.addfont(font_path)
    # Set the font properties
    prop = font_manager.FontProperties(fname=font_path)

    # Plot a stacked bar chart
    # set the plot height with ph
    
    fig, ax = plt.subplots(figsize=(25, ph)) 
    plt.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['pdf.fonttype'] = 42

    pivoted_df.plot(kind='bar', stacked=True, color=[chunk_color for chunk_color in pivoted_df.columns], ax=ax)


    # add * to the xpos that the sig_list both column is T
    # * should be right on top of the cooresponding bar

    # Convert xpos to integer for plotting
    sig_list['xpos'] = sig_list['xpos'].astype(int)

    # Plot * for each row where 'both' is 1
    for _, row in sig_list[sig_list['both'] == 1].iterrows():
        plt.text(row['xpos'], row['barht'], '*', fontsize=38, ha='center', va='center')

    # Add a black horizontal line at each xpos where 'hotspot' is True
    hotspot_df = sig_list[sig_list['hotspot']]
    for i in range(len(hotspot_df)-1):
        current_row = hotspot_df.iloc[i]
        next_row = hotspot_df.iloc[i+1]
        if current_row['group'] == next_row['group']:
            plt.hlines(y=max(sig_list['barht'])+addht, 
                       xmin=current_row['xpos'], xmax=next_row['xpos'], 
                       color='black',linewidth=5)


            
    plt.xlabel(f'{seq2name} chunk', fontsize=38)
    #plt.xticks(fontsize=8)
    if seq2name == 'Xist' :
        plt.xticks(fontsize=38, rotation=45, ha='right')
    else: 
        plt.xticks(np.arange(0, len(pivoted_df.index), 10),
               [chunk.split('_')[0] for chunk in pivoted_df.index[::10]],
               fontsize=38, rotation=45)
    plt.ylabel(f'Sig {seq1name}\nChunk Count', fontsize=38)
    
    plt.yticks(fontsize=38)
    # Use MaxNLocator to set y ticks to integers
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    
    # Format the y-axis tick labels to have leading zeros
    def format_func(value, tick_number):
        return "{:02d}".format(int(value))  # This will add leading zeros to numbers less than 10

    formatter = ticker.FuncFormatter(format_func)
    ax.yaxis.set_major_formatter(formatter)
    #plt.title(f'Stacked Histogram of similar (p<{pval_cutoff}) {seq1name} chunk count per {seq2name} chunk',fontweight='bold', fontsize=12)
    
    # Create a custom legend
    # patches = [mpatches.Patch(color=color, label=feature) for color, feature in color_feature.items()]
    #ax.legend(handles=patches, title='Sequence Feature', loc='lower left', bbox_to_anchor=(0.5, 0.5))
    # ax.legend(handles=patches, prop={'size':26, 'weight':'bold'},
    #           loc='upper left', bbox_to_anchor=(1.0, 1.0))
    # hide legend
    ax.legend().set_visible(False)
  
    plt.tight_layout()
    # save figure as pdf file
    #plt.savefig('test.pdf',dpi=300, format='pdf')
    plt.savefig(f'{seq1name}_{seq2name}_500nt_k4_{q_cutoff}_seekr_hist_manxist_2.pdf', dpi=300, format='pdf')

    #plt.show()

plot_transcripts('Xist', chunkcol='seq1_chunk', refcol='seq2_chunk',
                 sim_df_path='mXist_airn_500nt_k4_seekr_rscore_qlist_manxist.csv', 
                 sig_list_path='airn_sig_list.csv',
                 q_cutoff='q975', ph=6, addht=1)

plot_transcripts('Xist', chunkcol='seq1_chunk', refcol='seq2_chunk', 
                 sim_df_path='mXist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                 sig_list_path='kcnq1ot1_sig_list.csv', 
                 q_cutoff='q975',ph=6, addht=1)


plot_transcripts('Airn', chunkcol='seq2_chunk', refcol='seq1_chunk',
                 sim_df_path='mXist_airn_500nt_k4_seekr_rscore_qlist_manxist.csv', 
                 sig_list_path='mXist_sig_list.csv',
                 q_cutoff='q975',ph=7.2, addht=1)

plot_transcripts('Kcnq1ot1', chunkcol='seq2_chunk', refcol='seq1_chunk',
                 sim_df_path='mXist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                 sig_list_path='mXist_sig_list.csv', 
                 q_cutoff='q975',ph=7.2, addht=5)



def plot_AKtranscripts(seq1name, chunkcol, refcol,
                       sim_df_path, sig_list_path, q_cutoff, addht):
    
    
    # seq1 = seekrReader(seq1path).get_seqs()[0]
    # total_length1 = len(seq1)

    # general_color_feature=general_color_feature_dict(total_length1)

    # load in saved csv file
    sim_df_lg = pd.read_csv(sim_df_path)
    # sig_chunks = sim_df_lg[sim_df_lg['pval'] < pval_cutoff]

    # add a 'color' column to sim_df_lg
    sim_df_lg['chunk_color'] = ''

    # parse seq2_chunk column to get seq2name
    seq2name = sim_df_lg[refcol][0].split('_')[0]

    # if seq2name is mXist, change it to Xist
    # else if seq2name is airn, change it to Airn
    if seq2name == 'mXist':
        seq2name = 'Xist'
    elif seq2name == 'airn':
        seq2name = 'Airn'
    elif seq2name == 'kcnq1ot1':
        seq2name = 'Kcnq1ot1'

    
    if seq1name == 'Airn':
        color_dict = airn_color_dict
        coord_dict = airn_coord_dict
    elif seq1name == 'Kcnq1ot1':
        color_dict = k_color_dict
        coord_dict = k_coord_dict


    # get rid of the first word after split with _ from seq2_chunk column
    sim_df_lg[refcol] = sim_df_lg[refcol].apply(lambda x: '_'.join(x.split('_')[1:]))
    
    # generate color for seq1 chunk with name 'mXist' using Xist_assign_color function
    # if not named 'mXist', use general_assign_color function
    
    # color_feature=general_color_feature

    # for each row of sim_df_lg
    # parse c_start and c_end from seq1_chunk column
    # assign color to the row based on general_assign_color function
    for index, row in sim_df_lg.iterrows():
        c_start = int(row[chunkcol].split('_')[1])
        c_end = int(row[chunkcol].split('_')[2])
        sim_df_lg.loc[index, 'chunk_color'] = general_assign_color(c_start, c_end, coord_dict, color_dict)


    # Assuming your DataFrame is named df.
    # Firstly, filter the DataFrame to keep rows where pval < pval_cutoff
    filtered_df = sim_df_lg[sim_df_lg[q_cutoff] == True]

    # Group by 'seq2_chunk' and 'color', then count the number of rows in each group
    grouped_df = filtered_df.groupby([refcol, 'chunk_color']).size().reset_index(name='counts')

    # Pivot the DataFrame so that 'seq2_chunk' is the index, 'color' columns, and 'counts' are the values
    pivoted_df = grouped_df.pivot(index=refcol, columns='chunk_color', values='counts')

    # Reindex the pivoted DataFrame with the unique 'seq2_chunk' values from the original DataFrame
    unique_seq2_chunks = sim_df_lg[refcol].unique()
    pivoted_df = pivoted_df.reindex(unique_seq2_chunks)

    # Sort the DataFrame by 'seq2_chunk'
    # pivoted_df = pivoted_df.sort_index()

    # Fill any NaN values with 0
    pivoted_df = pivoted_df.fillna(0)

    # Filter columns from color_order that exist in pivoted_df
    filtered_colors = [color for color in color_order if color in pivoted_df.columns]

    # Reindex pivoted_df based on the filtered_colors
    pivoted_df = pivoted_df[filtered_colors]

    sig_list = pd.read_csv(sig_list_path)
    # convert to intege
    sig_list['both'] = sig_list['both'].astype(int)
    sig_list['barht'] = sig_list['chunk'].apply(lambda x: pivoted_df.loc[x].sum() if x in pivoted_df.index else 0)

    # Create an identifier for each group of consecutive 'T's
    sig_list['group'] = sig_list['both'].diff().ne(0).cumsum()

    # Calculate group sizes
    group_sizes = sig_list.groupby('group')['both'].transform('size')

    # Assign True to the new 'hotspot' column where group size is >=5
    sig_list['hotspot'] = (group_sizes >= 3) & (sig_list['both'] == 1)

    # use Airal for fonts
    # Specify the path to your .ttf font file
    font_path = 'arial.ttf'
    # Register the font with matplotlib
    font_manager.fontManager.addfont(font_path)
    # Set the font properties
    prop = font_manager.FontProperties(fname=font_path)

    # Plot a stacked bar chart
    fig, ax = plt.subplots(figsize=(25, 6)) 
    plt.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['pdf.fonttype'] = 42
    pivoted_df.plot(kind='bar', stacked=True, color=[chunk_color for chunk_color in pivoted_df.columns], ax=ax)

    # add * to the xpos that the sig_list both column is T
    # * should be right on top of the cooresponding bar

    # Convert xpos to integer for plotting
    sig_list['xpos'] = sig_list['xpos'].astype(int)

    # Plot * for each row where 'both' is 1
    for _, row in sig_list[sig_list['both'] == 1].iterrows():
        plt.text(row['xpos'], row['barht'], '*', fontsize=38, ha='center', va='center')

    # Add a black horizontal line at each xpos where 'hotspot' is True
    hotspot_df = sig_list[sig_list['hotspot']]
    for i in range(len(hotspot_df)-1):
        current_row = hotspot_df.iloc[i]
        next_row = hotspot_df.iloc[i+1]
        if current_row['group'] == next_row['group']:
            plt.hlines(y=max(sig_list['barht'])+addht, 
                       xmin=current_row['xpos'], xmax=next_row['xpos'], 
                       color='black',linewidth=5)


    plt.xlabel(f'{seq2name} chunk', fontsize=38)
    #plt.xticks(fontsize=8)
    # lable x axis ticks with every 10th chunk
    # ticks label is the first element of pivoted_df.index split by _
    plt.xticks(np.arange(0, len(pivoted_df.index), 10),
                [chunk.split('_')[0] for chunk in pivoted_df.index[::10]],
                fontsize=38, rotation=45)

    #plt.xticks(fontsize=5)
    plt.ylabel(f'Sig {seq1name}\nChunk Count', fontsize=38)
    plt.yticks(fontsize=38)
    # Use MaxNLocator to set y ticks to integers
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    # plt.title(f'Stacked Histogram of similar (p<{pval_cutoff}) {seq1name} chunk count per {seq2name} chunk',fontweight='bold', fontsize=12)
    
    # Format the y-axis tick labels to have leading zeros
    def format_func(value, tick_number):
        return "{:02d}".format(int(value))  # This will add leading zeros to numbers less than 10

    formatter = ticker.FuncFormatter(format_func)
    ax.yaxis.set_major_formatter(formatter)
    # Create a custom legend
    # patches = [mpatches.Patch(color=color, label=feature) for color, feature in color_feature.items()]
    #ax.legend(handles=patches, title='Sequence Feature', loc='lower left', bbox_to_anchor=(0.5, 0.5))
    # set legend outside of the plot region
    # ax.legend(handles=patches, prop={'size':26, 'weight':'bold'},
    #           loc='upper left', bbox_to_anchor=(1.0, 1.0))
    # hide legend
    ax.legend().set_visible(False)
  
    plt.tight_layout()
    # save figure as pdf file
    plt.savefig(f'{seq1name}_{seq2name}_500nt_k4_{q_cutoff}_seekr_hist_manxist_2.pdf', dpi=300, format='pdf')

    #plt.show()

plot_AKtranscripts('Airn', chunkcol='seq1_chunk', refcol='seq2_chunk',
                   sim_df_path='airn_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                   sig_list_path='kcnq1ot1_sig_list.csv',
                   q_cutoff='q975', addht=1)


plot_AKtranscripts('Kcnq1ot1', chunkcol='seq2_chunk', refcol='seq1_chunk',
                   sim_df_path='airn_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',
                   sig_list_path='airn_sig_list.csv',
                   q_cutoff='q975',addht=10)


####################################################################

# plot colored rectangles for xist, airn or k according to their length
# color code the line based on the feature
# label coordinates of the feature
# for xist label the feature name in the middle of the line

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import numpy as np
import matplotlib as mpl

# use Airal for fonts
# Specify the path to your .ttf font file
font_path = 'arial.ttf'
# Register the font with matplotlib
font_manager.fontManager.addfont(font_path)
# Set the font properties
prop = font_manager.FontProperties(fname=font_path)


# Define chunks
chunks = {
    'interval1': [0,354],
    'repeatA': [354,745],
    'ss234': [745,1497],
    'repeatFdwn': [1497,2859],
    'repeatB': [2859,3080],
    'repeatC': [3080,4692],
    'interval2': [4692,10213],
    'repeatEbroad': [10213,11630],
    'interval3': [11630,17946]
}

# Define colors for chunks
colors = {
    'interval1': '#1f77b4',
    'repeatA': '#ff7f0e',
    'ss234': '#2ca02c',
    'repeatFdwn': '#d62728',
    'repeatB': '#9467bd',
    'repeatC': '#8c564b',
    'interval2': '#e377c2',
    'repeatEbroad': '#7f7f7f',
    'interval3':'#bcbd22'
}

# Create figure and axis
fig, ax = plt.subplots(figsize=(45, 3))
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = prop.get_name()
ax.set_xlim([0, 17950])

# Iterate over chunks
for chunk, (start, end) in chunks.items():
    
    rect = patches.Rectangle((start, 0), (end-start), 0.3, facecolor=colors[chunk])
    ax.add_patch(rect)
    # add coordinates at the start of each chunk
    # but not add the 2nd and 5th chunk
    if chunk not in ['ss234', 'repeatC']:
        ax.text(start, 0.31, str(start), fontsize=50, ha='left')
    # plot the 2nd and 5th chunk coordinates a little higher to avoid overlap
    else:
        ax.text(start, 0.07, str(start), fontsize=50,  ha='left')

    # add feature name below the chunk
    ax.text((start + end) / 2, -0.1, chunk, 
            ha='right', va='top', fontsize=56, 
            rotation=45)

    # add 'mXist' in the middle below the line
    ax.text(8973, 0.6, 'Xist', ha='center', fontsize=56)

# add the last end coordinate
ax.text(17947, 0.31, str(17947), fontsize=50, ha='left')     
# Remove y-axis
ax.yaxis.set_visible(False)
# remove xaxis and all borders
ax.xaxis.set_visible(False)
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

#plt.savefig('mXist_500nt_k4_seekr_hist_legend.pdf', dpi=300, format='pdf')
# When you save the figure with plt.savefig(), matplotlib by default only includes the area within the axes
# feature names are plotted under chunks with coordinates as -0.1
# it won't be included for saving
# you can adjust the area included in the output file with the bbox_inches parameter
# if the 'tight' argument cuts off element near the edge
# plt.savefig('mXist_500nt_k4_seekr_hist_legend.pdf', dpi=300, format='pdf', bbox_inches=matplotlib.transforms.Bbox([[x0, y0], [x1, y1]]))
# x0, y0 is the lower-left corner of the bounding box, and x1, y1 is the upper-right corner. 
plt.savefig('mXist_500nt_k4_seekr_hist_legend_manxist.pdf', dpi=300, format='pdf', bbox_inches='tight')

#################
# plot airn legend

airn_color_dict = {
    'a1': '#1f77b4',
    'a2': '#ff7f0e',
    'a3': '#2ca02c',
    'a4': '#d62728',
    'a5': '#9467bd',
    'a6': '#8c564b',
    'a7': '#e377c2',
    'a8': '#7f7f7f',
    'a9': '#bcbd22',
    'a10': '#17becf',
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

# Create figure and axis
fig, ax = plt.subplots(figsize=(45, 3))
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = prop.get_name()
ax.set_xlim([0, 88754])

# Iterate over airn_coord_dict to get start and end
# get corresponding color from airn_color_dict
converted_dict = {airn_color_dict[key]: (value[0], value[1]) for key, value in airn_coord_dict.items()}
for color, (start, end) in converted_dict.items():
    
    rect = patches.Rectangle((start, 0), (end-start), 0.3, facecolor=color)
    ax.add_patch(rect)
    # add coordinates at the start of each chunk
    ax.text(start, 0.31, str(start), ha='left', fontsize=50)

    # add 'airn' in the middle below the line
    ax.text(44375, 0.6, 'Airn', ha='center', fontsize=56)

# add the last end coordinate
ax.text(88754, 0.31, str(88754), ha='left', fontsize=50)     
# Remove y-axis
ax.yaxis.set_visible(False)
# remove xaxis and all borders
ax.xaxis.set_visible(False)
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

plt.savefig('airn_500nt_k4_seekr_hist_legend_manxist.pdf', dpi=300, format='pdf')


#############
# plot k legend


k_color_dict = {
    'k1': '#1f77b4',
    'k2': '#ff7f0e',
    'k3': '#2ca02c',
    'k4': '#d62728',
    'k5': '#9467bd',
    'k6': '#8c564b',
    'k7': '#e377c2',
    'k8': '#7f7f7f',
    'k9': '#bcbd22',
    'k10': '#17becf',
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


# Create figure and axis
fig, ax = plt.subplots(figsize=(45, 3))
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = prop.get_name()
ax.set_xlim([0, 93092])

# Iterate over chunks
converted_dict = {k_color_dict[key]: (value[0], value[1]) for key, value in k_coord_dict.items()}
for color, (start, end) in converted_dict.items():
    
    rect = patches.Rectangle((start, 0), (end-start), 0.3, facecolor=color)
    ax.add_patch(rect)
    # add coordinates at the start of each chunk
    ax.text(start, 0.31, str(start), ha='left', fontsize=50)

    # add 'kcnq1ot1' in the middle below the line
    ax.text(46545, 0.6, 'Kcnq1ot1', ha='center', fontsize=56)

# add the last end coordinate
ax.text(93092, 0.31, str(93092), ha='left', fontsize=50)     
# Remove y-axis
ax.yaxis.set_visible(False)
# remove xaxis and all borders
ax.xaxis.set_visible(False)
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

plt.savefig('kcnq1ot1_500nt_k4_seekr_hist_legend_manxist.pdf', dpi=300, format='pdf')


