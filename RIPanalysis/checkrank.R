df<-read.csv('all34pr_kallisto_igg_rpm_TPMs_updated_04_10_24.csv')
ma<-read.csv('all34pr_kallisto_igg_rpm_TPMs_updated_11_28_23_nodups_gc.csv')

df_sorted <- df[order(df$alyref_rpm_over_igg, decreasing = TRUE), ]
df_sorted$gene_ID[1:10]

# Get the rank of specific genes 
match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
        "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
        "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"), df_sorted$gene_ID)

ma_sorted <- ma[order(ma$alyref_rpm_over_igg, decreasing = TRUE), ]
ma_sorted$gene_ID[1:10]

match(c("Xist_chrX_103460366_103483254(-)_transcript=ENSMUST00000127786.3", 
        "Airn_chr17_12741311_12860136(+)_transcript=unspliced(+)", 
        "Kcnq1ot1_chr7_143203458_143296549(-)_transcript=unspliced(-)"), ma_sorted$gene_ID)

dfcol<-colnames(df)[grepl('_rpm_over_igg',colnames(df))]

rankdf<-as.data.frame(matrix(nrow=3,ncol=68))

pname<-strsplit(dfcol,'_',fixed=T)
pname<-sapply(pname,'[[',1)
pname<-paste0(rep(pname, each=2), c("_new", "_old"))
colnames(rankdf)<-pname
rownames(rankdf)<-c('Xist','Airn','Kcnq1ot1')

for (cn in 1:length(dfcol)) {
  
  tcol<-dfcol[cn]
  
  df_sorted <- df[order(df[,tcol], decreasing = TRUE), ]
  
  # Get the rank of specific genes 
  rankdf[,(2*cn-1)]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                             "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                             "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$gene_ID)
  
  if (tcol %in% colnames(ma)) {
    ma_sorted <- ma[order(ma[,tcol], decreasing = TRUE), ]
    
    # Get the rank of specific genes 
    rankdf[,(2*cn)]<-match(c("Xist_chrX_103460366_103483254(-)_transcript=ENSMUST00000127786.3", 
                             "Airn_chr17_12741311_12860136(+)_transcript=unspliced(+)", 
                             "Kcnq1ot1_chr7_143203458_143296549(-)_transcript=unspliced(-)"),ma_sorted$gene_ID)
  }
  
  
}

write.csv(rankdf, 'XAK_ranks_kallisto_igg_rpm_TPMs_updated_04_10_24.csv')
