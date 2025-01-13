# organize data of chunk seekr comparisons between XAK
# Pin on X and assign all q975 chunks of A and K to each chunk of X, 
# with chunkname(rval); chunkname(rval) 
# and sequences of X and ATGTGTATC;ATGTCTGTA as A and K sequences. 
# Do the same for AK.
library(seqinr)
library(Biostrings)


xk<-read.csv('Xist_kcnq1ot1_500nt_k4_seekr_rscore_qlist_manxist.csv',header=T)

# filer to get only q975 True rows
xk.sig<-xk[which(xk$q975=='True'),]

xk.sig$q95<-NULL
xk.sig$q97<-NULL
xk.sig$q975<-NULL
xk.sig$q98<-NULL

xk.sig$gene1<-'XIST'
xk.sig$gene2<-'KCNQ1OT1'

colnames(xk.sig)[1:2]<-c('chunk1','chunk2')

xist<-read.fasta('XIST_manual_chunks_withF.fa',forceDNAtolower = F)
xist<-as.data.frame(getName(xist))
colnames(xist)<-'chunk'
xist$pos<-c(1:nrow(xist))

ot1<-read.fasta('KCNQ1OT1_chunks.fa',forceDNAtolower = F)
ot1<-as.data.frame(getName(ot1))
colnames(ot1)<-'chunk'
ot1$pos<-c(1:nrow(ot1))

xk.sig$pos1<-NA
xk.sig$pos2<-NA

for (n in 1:nrow(xk.sig)) {
  
  xk.sig$pos1[n]<-xist$pos[which(xist$chunk==xk.sig$chunk1[n])]
  xk.sig$pos2[n]<-ot1$pos[which(ot1$chunk==xk.sig$chunk2[n])]
  
}


xk.sig<-xk.sig[,c(4,1,6,5,2,7,3)]

write.csv(xk.sig,'XIST_OT1_q975chunk.csv',row.names = F)

