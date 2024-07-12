# simulation with random selected chunks from gencode chunks
# and check how many chunk pairs pass the threshold of rscores
# test whether the num of chunk pairs pass threhold of X vs A is sig

import random
import numpy as np
import pandas as pd
import os

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

import pickle



def count_string(s):
    # check if the string is NaN
    if pd.isna(s):
        return np.array([]), np.array([])

    # split the string by ';' and then by '_' 
    parts = [part.split('_')[0] for part in s.split(';')]

    # count the number of each unique part
    unique_parts = np.unique(parts, return_counts=True)
    # get the cooresponding unique parts
    u_name=unique_parts[0]
    u_counts=unique_parts[1]

    # return the unique rna names
    return u_name, u_counts

# save all element of a list of seqs into one fasta file
# close the files after written to avoid too many open files error

def save_seqs_to_fasta(seqs, seq_names, seqpath):
    seqfile = open(seqpath, 'w')
    for i, seq in enumerate(seqs):
        # save fasta name from corresponding element in seq_names
        seqfile.write('>'+seq_names[i] + '\n' + seq + '\n')
    seqfile.close()



def chunk_r_simulation (seq1path, seq2path, simu_n,
                        qpath, mean_path, std_path, sigfilepath):
    

    seq1counter = seekrBasicCounter(seq1path,mean=mean_path,std=std_path,k=4,silent=True)
    seq2counter = seekrBasicCounter(seq2path,mean=mean_path,std=std_path,k=4,silent=True)
    seq1counter.make_count_file()
    seq2counter.make_count_file()

    # get chunk headers
    seq1headers = seekrReader(seq1path).get_headers()
    seq2headers = seekrReader(seq2path).get_headers()
    # remove the '>' in the headers
    seq1headers = [header[1:] for header in seq1headers]
    seq2headers = [header[1:] for header in seq2headers]

    # Find similarities
    simu_r = seekrPearson(seq1counter.counts,seq2counter.counts)

    # Store the original shape
    simu_shape = simu_r.shape


    # convert simu_r to dataframe
    # simu_r = pd.DataFrame(simu_r, index=seq1headers, columns=seq2headers)

    # read in a q975 chunks files
    adf=pd.read_csv(sigfilepath)
    # drop airn_seq and sig_seq columns
    adf.drop(['airn_seq','sig_seq'], axis=1, inplace=True)

    # add column xcount and kcount to adf
    adf['mXist'] = 0
    adf['kcnq1ot1'] = 0

    # loop through each row in adf
    for n in range(len(adf)):
        u_name, u_counts = count_string(adf.iloc[n]['sig_chunk'])
        if len(u_name) == 1:
            adf.loc[n, u_name[0]] = u_counts[0]
        elif len(u_name) == 2:
            adf.loc[n, u_name[0]] = u_counts[0]
            adf.loc[n, u_name[1]] = u_counts[1]
    
    # set airn_chunk to be the index of adf 
    # so that later the kmask can be applied to other dataframes
    adf.set_index('airn_chunk', inplace=True)

    # generate the mask that has True for places that has Kcnq1ot1 > 0
    # and False for places that has Kcnq1ot1 == 0
    
    kmask = adf['kcnq1ot1'] > 0

    # get adf counts for Xist
    # first apply kmask to xist counts
    adf_counts = adf[kmask]['mXist']
    # then calculate total as ones that has > 0 counts
    # also count 1,2,3,4
    adf_total=sum(adf_counts > 0)
    adf_1=sum(adf_counts == 1)
    adf_2=sum(adf_counts == 2)
    adf_3=sum(adf_counts == 3)
    adf_4=sum(adf_counts == 4)

    qdict = pickle.load(open(qpath, 'rb'))
    rthresh = qdict['q975']

    # load in previous stored data
    # gencode_chunks = seekrReader(bkgpath).get_seqs()

    # initiate a list to store the mXist counts for each simulation
    # total counts, count 1 for each chunk, count 2,3,4 for each chunk
    share_total = np.zeros(simu_n)
    share_1 = np.zeros(simu_n)
    share_2 = np.zeros(simu_n)
    share_3 = np.zeros(simu_n)
    share_4 = np.zeros(simu_n)


    # do random simulation n times for the same length of seqs
    for i in range(simu_n):
        #print('simulation round:' + str(i))
        print('simulation round:' + str(i)) if i%100 == 0 else None
        # start = time.time()

        #randomize all values in simu_r matrix

        # Flatten the array
        flat_simu = simu_r.flatten()
        # Shuffle the flattened array
        np.random.shuffle(flat_simu)

        # Reshape the array back to its original shape
        simu_array = flat_simu.reshape(simu_shape)

        # convert simu_r to dataframe
        simu_array = pd.DataFrame(simu_array, index=seq1headers, columns=seq2headers)

        # for each row get the number of chunks that pass the rthresh
        simu_sig = simu_array.apply(lambda x: sum(x >= rthresh), axis=1)

        # filter with kmask
        simu_counts=simu_sig[kmask]

        # get total counts, count 1,2,3,4 for each chunk
        # and append to the lists
        share_total[i]=sum(simu_counts > 0)
        share_1[i]=sum(simu_counts == 1)
        share_2[i]=sum(simu_counts == 2)
        share_3[i]=sum(simu_counts == 3)
        share_4[i]=sum(simu_counts == 4)

    
    # simu_num_sig=np.array(simu_num_sig)
    # seq_r is a 1D np array
    # simu_r is a 2D np arrays with rows as simulations and columns as r scores
    
    pval_total=np.sum(share_total > adf_total)/simu_n
    pval_1=np.sum(share_1 > adf_1)/simu_n
    pval_2=np.sum(share_2 > adf_2)/simu_n
    pval_3=np.sum(share_3 > adf_3)/simu_n
    pval_4=np.sum(share_4 > adf_4)/simu_n

    return pval_total, pval_1, pval_2, pval_3, pval_4


p_total,p1,p2,p3,p4 = chunk_r_simulation (seq1path='airn_chunks.fa', 
                                          seq2path='mXist_chunk_manual.fa',
                                          simu_n=10000,
                                          qpath='gencodevM25_unique_divide_namedchunks_500_k4_rscores_withfullairn_q_dict.pkl',
                                          sigfilepath='airn_XAK_q975chunk.csv',
                                          mean_path='gencodevM25_unique_mean_4mers.npy',
                                          std_path='gencodevM25_unique_std_4mers.npy')

# p_total: 0.0021; p1: 0.1998; p2:0.1602; p3:0.0002; p4: 0.0