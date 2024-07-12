
import numpy as np
import pandas as pd
import pickle

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader



########################################################################################################

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
        seqfile.write('>'+seq_names[i] + '\n' + seq + '\n')
    seqfile.close()

#tc_h = ['>tc_' + str(index[0]) + '_' + str(index[1]) for index in tc_idx]
#save_seqs_to_fasta(tc,tc_h,'tc.fa')

###################### new function to cut Xist
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
                chunk_headers=[headers[i] + '_' + str(round(index[0]/500)) for index in chunk_index]
                div_seq_headers.append(chunk_headers)
            # if seq is shorter than target_length, append it to the list
            else:
                div_seqs.append([seq])
                # add seq length to the header
                div_seq_headers.append([headers[i] + '_0'])
    # flatten the list of lists into a list
    div_seqs = [item for sublist in div_seqs for item in sublist]
    div_seq_headers = [item for sublist in div_seq_headers for item in sublist]
    return div_seqs, div_seq_headers

# cut mXist interval 2 and 3 into chunks
xist='mxist_interval.fa'
xist_div_seqs, xist_div_seq_headers = divide_fasta_to_chunks(xist, 500, 200)
save_seqs_to_fasta(xist_div_seqs, xist_div_seq_headers, 'mxist_interval_chunks.fa')
# manually assemble xist chunks: xist repeat features with these interval chunks


########################## X vs A #########################
# calculate seekr r score between mXist and airn
# use mXist manually divided chunks
target_length=500
flex_length=200
mean_path='gencodevM25_unique_mean_4mers.npy'
std_path='gencodevM25_unique_std_4mers.npy'
k=4

# read in fasta files
seq2 = seekrReader('airn.fa').get_seqs()[0]

# divide seqs into chunks
seq2_chunks, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)

# convert chunk index to headers with adding seq name
seq1_chunk_headers = seekrReader('mXist_chunk_manual.fa').get_headers()
# remove the '>' from the headers
seq1_chunk_headers = [header[1:] for header in seq1_chunk_headers]
seq2_chunk_headers = [f'airn_{index[0]}_{index[1]}' for index in seq2_chunk_index]

seq2chunkpath='airn_chunks.fa'

# save chunks into fasta files
save_seqs_to_fasta(seq2_chunks, seq2_chunk_headers, seq2chunkpath)
    
# as BasicCounter is not designed to take in strings

# Count *k*-mers
chk1_counter = seekrBasicCounter('mXist_chunk_manual.fa',mean=mean_path,std=std_path,k=k,silent=True)
chk2_counter = seekrBasicCounter(seq2chunkpath,mean=mean_path,std=std_path,k=k,silent=True)
chk1_counter.make_count_file()
chk2_counter.make_count_file()

# Find similarities
sim = seekrPearson(chk1_counter.counts,chk2_counter.counts)

sim_df = pd.DataFrame(sim, index=seq1_chunk_headers, columns=seq2_chunk_headers)

# melt sim_df to long format
# with column name and row name as new columns
sim_df_lg = sim_df.reset_index().melt(id_vars='index')
sim_df_lg.columns = ['seq1_chunk', 'seq2_chunk', 'seekr_rscore']

# load in the qlist 
# load the dictionary from pickle file
with open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'rb') as f:
    q_dict = pickle.load(f)

# check whether each rscore passes the q95 threshold
sim_df_lg['q95'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q95'])

# check whether each rscore passes the q97 threshold
sim_df_lg['q97'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q97'])

# check whether each rscore passes the q975 threshold
sim_df_lg['q975'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q975'])

# check whether each rscore passes the q98 threshold
sim_df_lg['q98'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q98'])

# save sim_df to csv
sim_df_lg.to_csv(f'mXist_airn_{target_length}nt_k{k}_seekr_rscore_qlist_manxist.csv', index=False)



########################## X vs K #########################
# calculate seekr r score between mXist and K
# use mXist manually divided chunks
target_length=500
flex_length=200
mean_path='gencodevM25_unique_mean_4mers.npy'
std_path='gencodevM25_unique_std_4mers.npy'
k=4

# read in fasta files
seq2 = seekrReader('kcnq1ot1.fa').get_seqs()[0]

# divide seqs into chunks
seq2_chunks, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)

# convert chunk index to headers with adding seq name
seq1_chunk_headers = seekrReader('mXist_chunk_manual.fa').get_headers()
# remove the '>' from the headers
seq1_chunk_headers = [header[1:] for header in seq1_chunk_headers]
seq2_chunk_headers = [f'kcnq1ot1_{index[0]}_{index[1]}' for index in seq2_chunk_index]

seq2chunkpath='kcnq1ot1_chunks.fa'

# save chunks into fasta files
save_seqs_to_fasta(seq2_chunks, seq2_chunk_headers, seq2chunkpath)
    
# as BasicCounter is not designed to take in strings

# Count *k*-mers
chk1_counter = seekrBasicCounter('mXist_chunk_manual.fa',mean=mean_path,std=std_path,k=k,silent=True)
chk2_counter = seekrBasicCounter(seq2chunkpath,mean=mean_path,std=std_path,k=k,silent=True)
chk1_counter.make_count_file()
chk2_counter.make_count_file()

# Find similarities
sim = seekrPearson(chk1_counter.counts,chk2_counter.counts)

sim_df = pd.DataFrame(sim, index=seq1_chunk_headers, columns=seq2_chunk_headers)

# melt sim_df to long format
# with column name and row name as new columns
sim_df_lg = sim_df.reset_index().melt(id_vars='index')
sim_df_lg.columns = ['seq1_chunk', 'seq2_chunk', 'seekr_rscore']

# load in the qlist 
# load the dictionary from pickle file
with open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'rb') as f:
    q_dict = pickle.load(f)

# check whether each rscore passes the q95 threshold
sim_df_lg['q95'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q95'])

# check whether each rscore passes the q97 threshold
sim_df_lg['q97'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q97'])

# check whether each rscore passes the q975 threshold
sim_df_lg['q975'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q975'])

# check whether each rscore passes the q98 threshold
sim_df_lg['q98'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q98'])

# save sim_df to csv
sim_df_lg.to_csv(f'mXist_kcnq1ot1_{target_length}nt_k{k}_seekr_rscore_qlist_manxist.csv', index=False)

########################## A vs K #########################
# calculate seekr r score between airn and K
# should be same with previous file

def seq_rscore_csv(seq1path, seq2path, seq1name, seq2name,
                   target_length, flex_length, k, mean_path,
                   std_path):
    
    # read in fasta files
    seq1 = seekrReader(seq1path).get_seqs()[0]
    seq2 = seekrReader(seq2path).get_seqs()[0]

    # divide seqs into chunks
    seq1_chunks, seq1_chunk_index = divide_seq(seq1, target_length, flex_length)
    seq2_chunks, seq2_chunk_index = divide_seq(seq2, target_length, flex_length)

    # convert chunk index to headers with adding seq name
    seq1_chunk_headers = [f'{seq1name}_{index[0]}_{index[1]}' for index in seq1_chunk_index]
    seq2_chunk_headers = [f'{seq2name}_{index[0]}_{index[1]}' for index in seq2_chunk_index]

    seq1chunkpath=f'{seq1name}_chunks.fa'
    seq2chunkpath=f'{seq2name}_chunks.fa'

    # save chunks into fasta files
    save_seqs_to_fasta(seq1_chunks, seq1_chunk_headers, seq1chunkpath)
    save_seqs_to_fasta(seq2_chunks, seq2_chunk_headers, seq2chunkpath)
        
    # as BasicCounter is not designed to take in strings

    # Count *k*-mers
    chk1_counter = seekrBasicCounter(seq1chunkpath,mean=mean_path,std=std_path,k=k,silent=True)
    chk2_counter = seekrBasicCounter(seq2chunkpath,mean=mean_path,std=std_path,k=k,silent=True)
    chk1_counter.make_count_file()
    chk2_counter.make_count_file()

    # Find similarities
    sim = seekrPearson(chk1_counter.counts,chk2_counter.counts)

    sim_df = pd.DataFrame(sim, index=seq1_chunk_headers, columns=seq2_chunk_headers)
    
    # melt sim_df to long format
    # with column name and row name as new columns
    sim_df_lg = sim_df.reset_index().melt(id_vars='index')
    sim_df_lg.columns = ['seq1_chunk', 'seq2_chunk', 'seekr_rscore']

    # load in the qlist 
    # load the dictionary from pickle file
    with open('gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl', 'rb') as f:
        q_dict = pickle.load(f)

    # check whether each rscore passes the q95 threshold
    sim_df_lg['q95'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q95'])

    # check whether each rscore passes the q97 threshold
    sim_df_lg['q97'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q97'])

    # check whether each rscore passes the q975 threshold
    sim_df_lg['q975'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q975'])

    # check whether each rscore passes the q98 threshold
    sim_df_lg['q98'] = sim_df_lg['seekr_rscore'].apply(lambda x: x>q_dict['q98'])

    # subset sim_df_lg to only the first 10 rows
    # sim_df_lg = sim_df_lg.iloc[:10,:]

    # save sim_df to csv
    sim_df_lg.to_csv(f'{seq1name}_{seq2name}_{target_length}nt_k{k}_seekr_rscore_qlist_manxist.csv', index=False)


seq_rscore_csv('airn.fa', 'kcnq1ot1.fa', 'airn', 'kcnq1ot1', 
               target_length = 500, flex_length =200, k=4,
               mean_path = 'gencodevM25_unique_mean_4mers.npy',
               std_path = 'gencodevM25_unique_std_4mers.npy')

