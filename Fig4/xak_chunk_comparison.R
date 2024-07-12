#Load packages
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyverse)
library(data.table)
library(tidyr)
library(multcomp)


#Load featureCounts data for rip replicates
xak_fc <- read.table(file="xak_chunks_prc_new.txt", header=TRUE)
rownames(xak_fc) <- xak_fc$Geneid

#load total reads for rip replicates
ttlreads <- read.csv("totalreads.csv", header = TRUE)

#subset only rip data and replace data with 0 reads to 0.01
rip_reps <- xak_fc[7:30]
rip_reps <- replace(rip_reps, rip_reps < 1, 0.01)

# make sure the col of rip_reps and the row of ttlreads are the same
ttlreads$sample
colnames(rip_reps)

# divide each row of rip_reps by ttlreads$ttl_reads then x 1000000, then take log2
log2_rpm<-sweep(rip_reps,MARGIN = 2,ttlreads$ttl_reads,'/')
log2_rpm<-log2_rpm*1000000
log2_rpm<-log2(log2_rpm)
colnames(log2_rpm)<-gsub("Aligned_filteredsq30.out.sam", "log2_rpm",colnames(log2_rpm),fixed=T)


#subset data by lncRNA
x <- log2_rpm[c(1:5),]
a <- log2_rpm[c(6:10),]
k <- log2_rpm[c(11:15),]

# Get unique rip proteins from column names
prefixes <- unique(sub("_.*", "", colnames(log2_rpm)))




# Create a function to perform ANOVA and Tukey's HSD test
anova_and_tukey <- function(chunk_id, data) {
  # Subset the data for the current chunk
  chunk_a <- data %>% filter(Chunk == chunk_id)
  
  # Perform one-way ANOVA
  anova_result <- aov(value ~ condition, data = chunk_a)
  
  # Get the summary of the ANOVA
  anova_summary <- summary(anova_result)
  
  # Extract the p-value
  aov_pval <- anova_summary[[1]]$"Pr(>F)"[1]
  
  if (aov_pval<0.05) {
    # Perform Tukey's HSD test
    tukey_result <- TukeyHSD(anova_result, conf.level = 0.95)
    tukey_result <- as.data.frame(tukey_result$condition)
  } else {
    tukey_result <- NA
  }
  
  return(list(aov = aov_pval, adjpvals = tukey_result))
}


# make the plot and stat function -------------------------------------

stat_and_plot<-function(df,rnaname) {
  #create a row identifier column
  df$Chunk <- rownames(df)
  
  # clean up colnames for converting to long format
  colnames(df)<-gsub('_log2_rpm','',colnames(df))
  
  # Reshape the data to long format, separate the column names, and set condition factors
  long_x <- df %>%
    pivot_longer(cols = -c(Chunk), names_to = c("condition", "rep"), names_sep = "_")
  
  long_x$condition <- factor(long_x$condition, 
                             levels = c('igg','hnrnpk','ring1b','rybp','bmi1','ezh2','suz12','epop','mtf2','jarid2'))
  # Get the unique identifiers in the "Chunk" column
  unique_x_chunks <- unique(long_x$Chunk)
  
  rm(adjp_comb)
  rm(aov_comb)
  # Loop through unique chunks and perform the tests
  for(chunk_id in unique_x_chunks) {
    # Perform ANOVA and Tukey's HSD for the current chunk
    X_results <- anova_and_tukey(chunk_id, long_x)
    
    x_aov<-X_results$aov
    names(x_aov)<-chunk_id
    
    if (exists('aov_comb')) {
      aov_comb<-c(aov_comb,x_aov)
    } else {
      aov_comb<-x_aov
    }
    
    
    if (!is.na(X_results$adjpvals)[1]) {
      
      X_adjp <- X_results$adjpvals
      
      X_adjp$chunk<-chunk_id
      
      X_adjp$set<-rownames(X_adjp)
      
      rownames(X_adjp)<-NULL
      
      X_adjp<-X_adjp[,c(5,6,1:4)]
      
      if (exists('adjp_comb')) {
        adjp_comb<-rbind(adjp_comb,X_adjp)
      } else {
        adjp_comb<-X_adjp
      }
      
    }
    
  }
  
  # Write the aov and adj pval results to a csv file
  write.csv(aov_comb, paste0(rnaname,'_chunk_aov.csv'), quote=FALSE)
  write.csv(adjp_comb, paste0(rnaname,'_chunk_tukey.csv'), quote=FALSE, row.names = F)
  
  #find average of reps
  mean_x <- long_x %>%
    group_by(Chunk, condition) %>%
    summarise(avg = mean(value))
  
  #isolate IgG average
  igg_avg_x <- subset(mean_x, mean_x$condition == "igg")
  
  # filter long_x to keep only the sig ones that is compared to igg
  sig_comb<-adjp_comb[grepl('igg',adjp_comb$set),]
  
  sig_comb<-sig_comb[which(sig_comb$`p adj`<0.05),]
  
  sig_comb_temp<-strsplit(sig_comb$set,'-',fixed=T)
  
  sig_comb$rbp1<-sapply(sig_comb_temp,'[[',1)
  
  rm(sig_x)
  for (n in 1:nrow(sig_comb)) {
    
    ltemp<-long_x[which(long_x$Chunk==sig_comb$chunk[n] & long_x$condition==sig_comb$rbp1[n]),]
    
    ltemp<-as.data.frame(ltemp)
    
    if (exists('sig_x')) {
      sig_x<-rbind(sig_x,ltemp)
    } else {
      sig_x<-ltemp
    }
  }
  
  # Create a subset of data for highest points for each rbp in each chunk in sig_x
  highest_points <- sig_x %>%
    group_by(Chunk, condition) %>%
    summarize(highest_value = max(value)) %>%
    ungroup() %>%
    mutate(yvalue = highest_value+0.2)  
  
  # Reshape the data to long format and separate the column names
  plot_x <- long_x %>%
    group_by(Chunk, condition) %>% # separate into lncRNA chunks
    ggplot(aes(x = condition, y = value, color = condition)) +#plot rip xaxis and enrichment y axis, color code by rip
    #scale_y_log10()+
    geom_segment(data=mean_x, aes(xend = condition, x = condition, 
                                  yend = avg +0.075, y = avg -0.075, 
                                  linetype = "dashed", size = 1, color = condition)) + # add average bar to rips per chunk
    geom_point(aes(x = condition, y = value)) + #add individual rip replicate data points
    geom_text(data=highest_points, color='black',aes(x=condition, y=yvalue, label='*'), size=10) +
    geom_hline(data = igg_avg_x, aes(yintercept = avg), color = "black", linetype = "dashed", linewidth = 0.35) + #add line of igg average to compare other rip averages
    theme_bw() + #sets background of plots
    labs(title = paste0(rnaname, " Region"), x = "RIP", y = "Log2 RPM Enrichment") + #plot labels
    theme(plot.title=element_text(size=22),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x=element_text(size=14, angle = 45, hjust = 1),
          axis.text.y=element_text(size=14),
          strip.text = element_text(size = 16)
    ) + #make x axis labels angled
    facet_wrap(~ Chunk, nrow = 1) # Plots Xist from 5' -> 3'
  
  # Adjust the size of the plots
  plot_x = plot_x + theme(legend.position = "none") +
    expand_limits(y = c(min(long_x$value), max(long_x$value) + 1.5))
  
  # Print the updated plots with data points in a single large plot
  print(plot_x)
  
  #save plot
  ggsave(paste0(rnaname,"_chunk_prc_rips.pdf"), plot = plot_x, width = 12, height = 3, units = "in")
  
}



# Plot_XAK ---------------------------------------------------------------
stat_and_plot(x,'Xist')
stat_and_plot(a,'Airn')
stat_and_plot(k,'Kcnq1ot1')

