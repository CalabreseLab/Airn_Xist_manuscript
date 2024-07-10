# merge data from each protein to the final table

########################################################################

setwd("/work/users/s/h/shuang9/rip/RIP_data_process/kallisto_igg_rpm_all")
# get all file names
filenames<-dir(pattern='*.csv')

comb<-read.csv(filenames[1],header=T)

comb<-comb[,c(1,2,3,6,8)]

for (n in 2:length(filenames)) {
  
  df<-read.csv(filenames[n],header=T)
  
  #print(sum(comb$gene_ID==df$gene_ID)==137033)
  print(sum(comb$gene_ID==df$gene_ID)==114783)
  
  comb<-cbind(comb,df[,c(6,8)])
  
}

write.csv(comb,'all34pr_kallisto_igg_rpm_TPMs_updated_04_10_24.csv',row.names = F)
