# add chr and cyto tpms to the kallisto TPM data

###############################################
setwd("/work/users/s/h/shuang9/rip/RIP_data_process")

# merge df to prop
df<-read.csv('all34pr_kallisto_igg_rpm_TPMs_updated_04_10_24.csv',header=T)

prop<-read.csv('TSC-exp_ESC-exp_Chrom-assoc_4_9_2024.csv',header=T)

sum(df$gene_ID==prop$target_id)

merged<-merge(prop,df,by.x = "target_id", by.y = "gene_ID")

sum(merged$length.x==merged$length)
sum(merged$eff_length.x==merged$eff_length)

merged$length<-NULL
merged$eff_length<-NULL

write.csv(merged,'TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024.csv',row.names = F)

