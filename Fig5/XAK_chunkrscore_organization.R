# organize data of chunk seekr comparisons between XAK
# Pin on X and assign all q975 chunks of A and K to each chunk of X, 
# with chunkname(rval); chunkname(rval) 
# and sequences of X and ATGTGTATC;ATGTCTGTA as A and K sequences. 
# Do the same for AK.
library(seqinr)
library(Biostrings)


xa<-read.csv('mXist_airn_500nt_k4_seekr_rscore_qlist_manxist.csv',header=T)
xk<-read.csv('mXist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',header=T)
ak<-read.csv('airn_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',header=T)

# filer to get only q975 True rows
xa.sig<-xa[which(xa$q975=='True'),]
xk.sig<-xk[which(xk$q975=='True'),]
ak.sig<-ak[which(ak$q975=='True'),]

# get unique list of XAK
xchunk<-unique(xa$seq1_chunk)
achunk<-unique(ak$seq1_chunk)
kchunk<-unique(ak$seq2_chunk)

# initialize dataframe
x.df<-as.data.frame(xchunk)
colnames(x.df)<-'xist_chunk'
x.df$sig_chunk<-''
x.df$xist_seq<-''
x.df$sig_seq<-''

a.df<-as.data.frame(achunk)
colnames(a.df)<-'airn_chunk'
a.df$sig_chunk<-''
a.df$airn_seq<-''
a.df$sig_seq<-''

k.df<-as.data.frame(kchunk)
colnames(k.df)<-'kcnq1ot1_chunk'
k.df$sig_chunk<-''
k.df$k_seq<-''
k.df$sig_seq<-''

xist<-read.fasta('mXist.fa',forceDNAtolower = F)
xist<-unlist(getSequence(xist))
airn<-read.fasta('airn.fa',forceDNAtolower = F)
airn<-unlist(getSequence(airn))
kot1<-read.fasta('kcnq1ot1.fa',forceDNAtolower = F)
kot1<-unlist(getSequence(kot1))

# set xist coordinates
xistcoords <- list(
  mXist_interval1 = c(0,354),
  mXist_repeatA = c(354,745),
  mXist_ss234 = c(745,1497),
  mXist_repeatFdwn = c(1497,2859),
  mXist_repeatB = c(2859,3080),
  mXist_repeatC = c(3080,4692),
  mXist_interval2_0 = c(4692,5193),
  mXist_interval2_1 = c(5193,5694),
  mXist_interval2_2 = c(5694,6195),
  mXist_interval2_3 = c(6195,6696),
  mXist_interval2_4 = c(6696,7197),
  mXist_interval2_5 = c(7197,7698),
  mXist_interval2_6 = c(7698,8199),
  mXist_interval2_7 = c(8199,8700),
  mXist_interval2_8 = c(8700,9201),
  mXist_interval2_9 = c(9201,9702),
  mXist_interval2_10 = c(9702,10213),
  mXist_repeatEbroad = c(10213,11630),
  mXist_interval3_0 = c(11630,12115),
  mXist_interval3_1 = c(12115,12600),
  mXist_interval3_2 = c(12600,13085),
  mXist_interval3_3 = c(13085,13570),
  mXist_interval3_4 = c(13570,14055),
  mXist_interval3_5 = c(14055,14540),
  mXist_interval3_6 = c(14540,15025),
  mXist_interval3_7 = c(15025,15510),
  mXist_interval3_8 = c(15510,15995),
  mXist_interval3_9 = c(15995,16480),
  mXist_interval3_10 = c(16480,16965),
  mXist_interval3_11 = c(16965,17450),
  mXist_interval3_12 = c(17450,17946)
)


xak_getseq<-function(chunkname){
  chunkname_split<-unlist(strsplit(chunkname,'_'))
  if (chunkname_split[1]=='mXist') {
    start.coord<-xistcoords[chunkname][[1]][1]+1
    end.coord<-xistcoords[chunkname][[1]][2]
    temp.seq<-paste0(xist[start.coord:end.coord],collapse = "")
  } else if (chunkname_split[1]=='airn') {
    start.coord<-as.numeric(chunkname_split[2])+1
    end.coord<-as.numeric(chunkname_split[3])
    temp.seq<-paste0(airn[start.coord:end.coord],collapse = "")
  } else {
    start.coord<-as.numeric(chunkname_split[2])+1
    end.coord<-as.numeric(chunkname_split[3])
    temp.seq<-paste0(kot1[start.coord:end.coord],collapse = "")
  }
  
  return (temp.seq)
  
}


chunkdata<-function(df.temp,coln) {
  cn<-''
  cs<-''
  
  if (nrow(df.temp)>0) {
    
    for (i in 1:nrow(df.temp)) {
      cn.temp<-paste0(df.temp[[coln]][i],'(',round(df.temp$seekr_rscore[i],4),')')
      cn<-paste(c(cn,cn.temp),collapse=';')
      
      cs.temp<-xak_getseq(df.temp[[coln]][i])
      cs<-paste(c(cs,cs.temp),collapse=';')
    }
    
    cn<-substr(cn, 2, nchar(cn))
    cs<-substr(cs, 2, nchar(cs))
  }
  
  return(c(cn, cs))
  
}

# organize x
for (n in 1:nrow(x.df)) {
  
  xc<-x.df$xist_chunk[n]
  
  x.df$xist_seq[n]<-xak_getseq(xc)
  
  xa.temp<-xa.sig[which(xa.sig$seq1_chunk==xc),]
  xk.temp<-xk.sig[which(xk.sig$seq1_chunk==xc),]
  
  xa.res<-chunkdata(xa.temp,'seq2_chunk')
  xk.res<-chunkdata(xk.temp,'seq2_chunk')
  
  cname<-paste(c(xa.res[1],xk.res[1]),collapse =';')
  cseq<-paste(c(xa.res[2],xk.res[2]),collapse =';')
  
  # remove ; at the beginning or end if xa.res or xk.res is empty
  if (substr(cname,1,1)==';') {
    cname<-substr(cname, 2, nchar(cname))
  }
  
  if (substr(cseq,1,1)==';') {
    cseq<-substr(cseq, 2, nchar(cseq))
  }
  
  if (substr(cname,nchar(cname),nchar(cname))==';') {
    cname<-substr(cname, 1, nchar(cname)-1)
  }
  
  if (substr(cseq,nchar(cseq),nchar(cseq))==';') {
    cseq<-substr(cseq, 1, nchar(cseq)-1)
  }
  
  x.df$sig_chunk[n]<-cname
  x.df$sig_seq[n]<-cseq
  
  # write the target and related seqs into one fasta file for MEME analysis
  header_list <- unlist(strsplit(cname, ";"))
  seq_list <- unlist(strsplit(cseq, ";"))
  
  header_list<-c(xc,header_list)
  seq_list<-c(x.df$xist_seq[n],seq_list)
  
  # Check if lengths match
  if (length(header_list) != length(seq_list)) {
    stop("The number of headers and sequences do not match.")
  }
  
  # Create a DNAStringSet object
  dna_seqs <- DNAStringSet(seq_list)
  names(dna_seqs) <- header_list
  
  # Write to FASTA file
  writeXStringSet(dna_seqs, paste0(xc,".fasta"))
  
  
}

write.csv(x.df,'mXist_XAK_q975chunk.csv',row.names = F)



# organize a
for (n in 1:nrow(a.df)) {
  
  ac<-a.df$airn_chunk[n]
  
  a.df$airn_seq[n]<-xak_getseq(ac)
  
  ax.temp<-xa.sig[which(xa.sig$seq2_chunk==ac),]
  ak.temp<-ak.sig[which(ak.sig$seq1_chunk==ac),]
  
  ax.res<-chunkdata(ax.temp,'seq1_chunk')
  ak.res<-chunkdata(ak.temp,'seq2_chunk')
  
  cname<-paste(c(ax.res[1],ak.res[1]),collapse =';')
  cseq<-paste(c(ax.res[2],ak.res[2]),collapse =';')
  
  # remove ; at the beginning or end if xa.res or xk.res is empty
  if (substr(cname,1,1)==';') {
    cname<-substr(cname, 2, nchar(cname))
  }
  
  if (substr(cseq,1,1)==';') {
    cseq<-substr(cseq, 2, nchar(cseq))
  }
  
  if (substr(cname,nchar(cname),nchar(cname))==';') {
    cname<-substr(cname, 1, nchar(cname)-1)
  }
  
  if (substr(cseq,nchar(cseq),nchar(cseq))==';') {
    cseq<-substr(cseq, 1, nchar(cseq)-1)
  }
  
  a.df$sig_chunk[n]<-cname
  a.df$sig_seq[n]<-cseq
  
  # write the target and related seqs into one fasta file for MEME analysis
  header_list <- unlist(strsplit(cname, ";"))
  seq_list <- unlist(strsplit(cseq, ";"))
  
  header_list<-c(ac,header_list)
  seq_list<-c(a.df$airn_seq[n],seq_list)
  
  # Check if lengths match
  if (length(header_list) != length(seq_list)) {
    stop("The number of headers and sequences do not match.")
  }
  
  # Create a DNAStringSet object
  dna_seqs <- DNAStringSet(seq_list)
  names(dna_seqs) <- header_list
  
  # Write to FASTA file
  writeXStringSet(dna_seqs, paste0(ac,".fasta"))
  
  
}

write.csv(a.df,'airn_XAK_q975chunk.csv',row.names = F)



# organize k
for (n in 1:nrow(k.df)) {
  
  kc<-k.df$kcnq1ot1_chunk[n]
  
  k.df$k_seq[n]<-xak_getseq(kc)
  
  kx.temp<-xk.sig[which(xk.sig$seq2_chunk==kc),]
  ka.temp<-ak.sig[which(ak.sig$seq2_chunk==kc),]
  
  kx.res<-chunkdata(kx.temp,'seq1_chunk')
  ka.res<-chunkdata(ka.temp,'seq1_chunk')
  
  cname<-paste(c(kx.res[1],ka.res[1]),collapse =';')
  cseq<-paste(c(kx.res[2],ka.res[2]),collapse =';')
  
  # remove ; at the beginning or end if xa.res or xk.res is empty
  if (substr(cname,1,1)==';') {
    cname<-substr(cname, 2, nchar(cname))
  }
  
  if (substr(cseq,1,1)==';') {
    cseq<-substr(cseq, 2, nchar(cseq))
  }
  
  if (substr(cname,nchar(cname),nchar(cname))==';') {
    cname<-substr(cname, 1, nchar(cname)-1)
  }
  
  if (substr(cseq,nchar(cseq),nchar(cseq))==';') {
    cseq<-substr(cseq, 1, nchar(cseq)-1)
  }
  
  k.df$sig_chunk[n]<-cname
  k.df$sig_seq[n]<-cseq
  
  # write the target and related seqs into one fasta file for MEME analysis
  header_list <- unlist(strsplit(cname, ";"))
  seq_list <- unlist(strsplit(cseq, ";"))
  
  header_list<-c(kc,header_list)
  seq_list<-c(k.df$k_seq[n],seq_list)
  
  # Check if lengths match
  if (length(header_list) != length(seq_list)) {
    stop("The number of headers and sequences do not match.")
  }
  
  # Create a DNAStringSet object
  dna_seqs <- DNAStringSet(seq_list)
  names(dna_seqs) <- header_list
  
  # Write to FASTA file
  writeXStringSet(dna_seqs, paste0(kc,".fasta"))
  
  
  
}

write.csv(k.df,'kcnq1ot1_XAK_q975chunk.csv',row.names = F)
