# cut gencode and XKA into chunks
# use this to calculate r score for all unique 500 chunks to itself

import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
# import time

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
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

# tc, tc_idx=divide_seq('GGGSSSKKDDDHJASSSSIOOOSDUUEEEYYYHH', 10, 2)


# save all element of a list of seqs into one fasta file
# close the files after written to avoid too many open files error

def save_seqs_to_fasta(seqs, seq_names, seqpath):
    seqfile = open(seqpath, 'w')
    for i, seq in enumerate(seqs):
        # save fasta name from corresponding element in seq_names
        seqfile.write(seq_names[i] + '\n' + seq + '\n')
    seqfile.close()

#tc_h = ['>tc_' + str(index[0]) + '_' + str(index[1]) for index in tc_idx]
#save_seqs_to_fasta(tc,tc_h,'tc.fa')


################## cut gencode into chunck for simulation #############

# read in fasta file and divide all seqs using divide_seq function
# write each divided chunks into a list


def divide_fasta_to_chunks(seqpath, target_length, flex_length): 
    seqs = seekrReader(seqpath).get_seqs()
    headers=seekrReader(seqpath).get_headers()
    div_seqs = []
    div_seq_headers = []
    for i, seq in enumerate(seqs):
        # only keep seqs longer than flex_length
        if len(seq) > flex_length:
            # if seq is longer than target_length, divide it into chunks
            if len(seq) > target_length:
                chunks, chunk_index = divide_seq(seq, target_length, flex_length)
                div_seqs.append(chunks)
                # save the header for all chunks
                chunk_headers=[headers[i] + '_' + str(index[0]) + '_' + str(index[1]) for index in chunk_index]
                div_seq_headers.append(chunk_headers)
            # if seq is shorter than target_length, append it to the list
            else:
                div_seqs.append([seq])
                # add seq length to the header
                div_seq_headers.append([headers[i] + '_0_' + str(len(seq))])
    # flatten the list of lists into a list
    div_seqs = [item for sublist in div_seqs for item in sublist]
    div_seq_headers = [item for sublist in div_seq_headers for item in sublist]
    return div_seqs, div_seq_headers


gencode = 'gencode.vM25.lncRNA_transcripts.unique.genesequence_withfullairn.fa'
gencode_div_seqs, gencode_div_seq_headers = divide_fasta_to_chunks(gencode, 500, 200)
# 42150 chunks in total before
# now 42328 chunks in total

# save these gencode chunks into a fasta file
save_seqs_to_fasta(gencode_div_seqs, gencode_div_seq_headers, 'gencodevM25_unique_divide_namedchunks_500_withfullairn.fa')



##################################### generate gencode vs gencode r-scores ############################


gc_counter = seekrBasicCounter('gencodevM25_unique_divide_namedchunks_500_withfullairn.fa',
                               mean='gencodevM25_unique_mean_4mers.npy',
                               std='gencodevM25_unique_std_4mers.npy',
                               k=4,silent=True)
gc_counter.make_count_file()
gc_sim = seekrPearson(gc_counter.counts,gc_counter.counts)

# get the values in upper triangle of gc_sim and saved as 1D np array
gc_sim_triu = gc_sim[np.triu_indices(gc_sim.shape[0], k=1)]
# 895808628 values in gc_sim_triu

# check how many values in gc_sim_triu are 0
len(gc_sim_triu[gc_sim_triu==0])
# 10

# save gc_sim_triu as npy file
np.save('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn.npy', gc_sim_triu)


# convert gc_sim_upper to dataframe
import pandas as pd
gc_sim_df = pd.DataFrame(gc_sim)
# add column names and row names for gc_sim_df
u500_headers=seekrReader('gencodevM25_unique_divide_namedchunks_500_withfullairn.fa').get_headers()
# remove '>' from the header
u500_headers = [header[1:] for header in u500_headers]

gc_sim_df.columns = u500_headers
gc_sim_df.index = u500_headers

# check [1:10, 1:10] of gc_sim_upper_df
# gc_sim_upper_df.iloc[1:10, 1:10]

# melt gc_sim_df into a long format
# with column names: 'seq1', 'seq2', 'rscore'
gc_sim_df_lg = pd.melt(gc_sim_df.reset_index(), id_vars=['index'])
gc_sim_df_lg.columns = ['seq1', 'seq2', 'rscore']

# only keep rows in gc_sim_df_lg which belong to the upper triangle of gc_sim without diagonal
gc_sim_upper_df = gc_sim_df_lg[gc_sim_df_lg.seq1 < gc_sim_df_lg.seq2]

# check [1:10,] of gc_sim_upper_df
gc_sim_upper_df.iloc[0:10,]
gc_sim_upper_df.iloc[237023,]
gc_sim_df.loc['A830052D11Rik_354_708','Gm47518_434_868']


# save gc_sim_upper_df as pickle file
gc_sim_upper_df.to_pickle('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_dataframe.pkl')

# load gc_sim_upper_df from pickle file
gc_sim_upper_df = pd.read_pickle('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_dataframe.pkl')

# gc_sim_triu = gc_sim[np.triu_indices(gc_sim.shape[0], k=1)]
# import matplotlib.pyplot as plt
# plt.hist(gc_sim_triu, bins=300)

gc_np=np.load('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn.npy')
# get the 95% quantile of gc_np
q95=np.quantile(gc_np, 0.95)

################## get sig r-scores from gc_sim_upper_df #############

# get the 95% quantile of gc_sim_upper_df['rscore']
q95=gc_sim_upper_df['rscore'].quantile(0.95)

# get the 97% quantile of gc_sim_upper_df['rscore']
q97=gc_sim_upper_df['rscore'].quantile(0.97)

# get the 97.5% quantile of gc_sim_upper_df['rscore']
q975=gc_sim_upper_df['rscore'].quantile(0.975)

# get the 98% quantile of gc_sim_upper_df['rscore']
q98=gc_sim_upper_df['rscore'].quantile(0.98)

# save the quantiles into a dictionary
q_dict = {'q95':q95, 'q97':q97, 'q975':q975, 'q98':q98}
# save this dictionary as a pickle file
import pickle
with open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'wb') as f:
    pickle.dump(q_dict, f)

# load the dictionary from pickle file
with open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'rb') as f:
    q_dict = pickle.load(f)



