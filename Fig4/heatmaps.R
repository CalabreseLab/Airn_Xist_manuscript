# plot heatmaps of rankings and top1k pearsons

library(ggplot2)
library(ComplexHeatmap)
library(circlize)


########### ranking heatmap #################
df<-read.csv('TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024.csv',header=T)
# 114783 in total
# filter filter for length >=500, median exp >0.0625, chrom-fraction >0.75, 
# and also remove ERCC-spike in transcripts: ERCC-00171_transcript=ERCC-00171
df<-df[which(df$length.x >= 500 & df$Median_tpm >0.0625 & df$chrom_enrichment >0.75),]
# 19295
#df<-df[!grepl('ERCC-00',df$gene_ID,fixed=T),]


# get the avg exp of chr
# df$chr_mean<-(df$chr_r1_tpm+df$chr_r2_tpm)/2

write.csv(df,'TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024_filtered.csv',row.names = F)

# df<-read.csv('all34pr_kallisto_igg_rpm_TPMs_updated_11_28_23_nodups_gc_filtered.csv',header=T)
# expression total_tpm; expression of chr (avg (chr_r1_tpm	chr_r2_tpm)); 
# hnrnpk, ring1b,rybp,bmi1, ezh2,suz12,epop,jarid2,mtf2
rankdf<-as.data.frame(matrix(nrow=11,ncol=3))
colnames(rankdf)<-c('Xist','Airn','Kcnq1ot1')
rownames(rankdf)<-c('expression_total','expression_chr','hnrnpk','ring1b','rybp',
                    'bmi1','ezh2','suz12','epop','mtf2','jarid2')

######expression total median tpm
df_sorted <- df[order(df[,'Median_tpm'], decreasing = TRUE), ]

# Get the rank of specific genes 
rankdf[1,]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                    "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                    "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)

######expression chr
df_sorted <- df[order(df[,'chrom_average'], decreasing = TRUE), ]

# Get the rank of specific genes 
rankdf[2,]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                    "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                    "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)


#### rbps
rbplist<-c('hnrnpk','ring1b','rybp','bmi1','ezh2','suz12','epop','mtf2','jarid2')
rbplist<-paste0(rbplist,'_rpm_over_igg')

for (n in 1:length(rbplist)) {
  
  tcol<-rbplist[n]
  df_sorted <- df[order(df[,tcol], decreasing = TRUE), ]
  
  rankdf[(2+n),]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                          "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                          "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)
  
  
}

########### plot
log2_rankdf <- log2(rankdf)
# Create a color mapping from white to purple
color_mapping <- colorRamp2(c(min(log2_rankdf), max(log2_rankdf)), c("#34009a", "white"))

cn = colnames(rankdf)

# Create labels for the legend
legend_labels <- c(as.character(nrow(df)),'4096','512','64','8','1')
legend_values <- c(nrow(df),4096,512,64,8,1)
# Transform these values to match the scale of the heatmap
transformed_legend_values <- log2(legend_values)

log2_rankdf<-as.matrix(log2_rankdf)
# Open a PDF device
pdf("ranking_heatmap.pdf", width = 8, height = 8)  # Adjust 'width' and 'height' as needed

# Plot the heatmap
Heatmap(log2_rankdf, 
        col = color_mapping, 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "#e0e0e0", lwd = 1),
        show_row_names = TRUE,
        show_column_names = FALSE,
        #column_names_side = "top",
        row_names_side = "left",
        #column_names_gp = gpar(fontsize = 18,fontface = 3),
        row_names_gp = gpar(fontsize = 18),
        top_annotation = HeatmapAnnotation(
          text = anno_text(cn, rot = 0, location = unit(0.1, "npc"), just = "center",
                           gp = gpar(fontsize = 18, fontface = "italic")),
          annotation_height = max_text_width(cn, gp = gpar(fontsize = 18,fontface = "italic"))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          text_color <- ifelse(rankdf[i, j] <= 15, "white", "black")
          grid.text(rankdf[i, j], x, y, gp = gpar(col = text_color, fontsize = 18))
        },
        heatmap_legend_param = list(
          title = "Ranking among\nchromatin_enriched\ntranscripts",
          at = transformed_legend_values,
          labels = legend_labels,
          color_bar = "continuous",
          legend_direction = "vertical",
          legend_width = unit(6, "cm"),
          legend_height = unit(6, "cm"),
          labels_gp = gpar(fontsize = 18),
          title_gp = gpar(fontsize = 18))
        
)

dev.off()

############### top1k heatmap ######################

# filter df to only contain rbplist
top1k<-df[,rbplist]

# initiate matrix to store r value
topr<-as.data.frame(matrix(nrow=9,ncol=9))
rbpvec<-c('hnrnpk','ring1b','rybp','bmi1',
          'ezh2','suz12','epop','mtf2','jarid2')
rbplist
# check these two matches
colnames(topr)<-rbpvec
rownames(topr)<-rbpvec

for (p in 1:length(rbpvec)) {
  
  rbpcol<-rbplist[p]
  
  top1k_sorted <- top1k[order(top1k[,rbpcol], decreasing = TRUE), ]
  
  top1k_sorted<-top1k_sorted[1:1000,]
  
  sort_data<-top1k_sorted[,rbpcol]
  
  for (n in 1:ncol(topr)) {
    topr[p,n]<-cor(sort_data,top1k_sorted[,n],method='pearson')
  }
  
}


color_fun <- colorRamp2(c(0,0.5,1), c("white", "#ffff00","#a50000"))

# Create labels for the legend
legend_labels_t <- c('0','0.25','0.5','0.75','1')
legend_values_t <- c(0,0.25,0.5,0.75,1)

topr<-as.matrix(topr)

write.csv(topr, 'top1k_heatmap.csv')
# Open a PDF device
pdf("top1k_heatmap.pdf", width = 10, height = 8)  

# Plot the heatmap
Heatmap(topr, 
        col = color_fun, 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "#e0e0e0", lwd = 1),
        show_row_names = TRUE,
        show_column_names = FALSE,
        #column_names_side = "top",
        row_names_side = "left",
        #column_names_gp = gpar(fontsize = 18,fontface = 3),
        row_names_gp = gpar(fontsize = 24),
        top_annotation = HeatmapAnnotation(
          text = anno_text(rbpvec, rot = 90, location = unit(0.4, "npc"), just = "center",
                           gp = gpar(fontsize = 24)),
          annotation_height = max_text_width(cn, gp = gpar(fontsize = 24))
        ),
        heatmap_legend_param = list(
          title = "Pearson's r\n(top 1k transcripts)",
          at = legend_values_t,
          labels = legend_labels_t,
          color_bar = "continuous",
          legend_direction = "vertical",
          legend_width = unit(6, "cm"),
          legend_height = unit(6, "cm"),
          labels_gp = gpar(fontsize = 20),
          title_gp = gpar(fontsize = 20))
)

dev.off()



