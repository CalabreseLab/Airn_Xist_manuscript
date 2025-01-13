# plot to check whether there are structure conservation among XAK
# separately label chunks that are related to a Xist repeat or interval
# do the correlation with the repeats only chunks and all chunks
####################################

setwd("/to/the/folder/contains/XorAorK_XAK_q975chunk.csv")


library(ggplot2)

x<-read.csv('mXist_XAK_q975chunk.csv',header=T)
a<-read.csv('airn_XAK_q975chunk.csv',header=T)
k<-read.csv('kcnq1ot1_XAK_q975chunk.csv',header=T)

# add seq number
x$pos<-c(1:nrow(x))
a$pos<-c(1:nrow(a))
k$pos<-c(1:nrow(k))

# locate chunks that are shared among XAK
x_share<-x[(grepl('airn',x$sig_chunk) & grepl('kcnq1ot1',x$sig_chunk)),]
a_share<-a[(grepl('mXist',a$sig_chunk) & grepl('kcnq1ot1',a$sig_chunk)),]
k_share<-k[(grepl('mXist',k$sig_chunk) & grepl('airn',k$sig_chunk)),]

x_share$xist_seq<-NULL
x_share$sig_seq<-NULL
a_share$airn_seq<-NULL
a_share$sig_seq<-NULL
k_share$k_seq<-NULL
k_share$sig_seq<-NULL

colnames(x)[1]<-'xak_chunk'
colnames(a)[1]<-'xak_chunk'
colnames(k)[1]<-'xak_chunk'

colnames(x_share)[1]<-'xak_chunk'
colnames(a_share)[1]<-'xak_chunk'
colnames(k_share)[1]<-'xak_chunk'
# get the list of shared chunks
allchunks<-rbind(x[,c(1,5)],a[,c(1,5)],k[,c(1,5)])


# add xist feature color to denote repeat or interval
xistfeatures<-c('mXist_repeatA','mXist_ss234','mXist_repeatFdwn','mXist_repeatB',
                'mXist_repeatC','mXist_repeatEbroad')

x_share$xfeature<-''

for (n in 1:nrow(x_share)) {
  if (x_share$xak_chunk[n] %in% xistfeatures) {
    x_share$xfeature[n]<-'repeats'
  } else {
    x_share$xfeature[n]<-'intervals'
  }
}

a_share$xfeature<-''

for (n in 1:nrow(a_share)) {
  
  temp<-unlist(strsplit(a_share$sig_chunk[n],';',fixed=T))
  temp<-gsub("\\(.*$", "", temp, perl = TRUE)
  
  if (any(temp %in% xistfeatures)) {
    a_share$xfeature[n]<-'repeats'
  } else {
    a_share$xfeature[n]<-'intervals'
  }
}

k_share$xfeature<-''

for (n in 1:nrow(k_share)) {
  
  temp<-unlist(strsplit(k_share$sig_chunk[n],';',fixed=T))
  temp<-gsub("\\(.*$", "", temp, perl = TRUE)
  
  if (any(temp %in% xistfeatures)) {
    k_share$xfeature[n]<-'repeats'
  } else {
    k_share$xfeature[n]<-'intervals'
  }
}

# fucntion to loop thru the xak_share df and only keep the shared sig chunk

shared_chunks<-function(df) {
  
  tdf<-data.frame(gene1=character(),
                  chunk1=character(),
                  pos1=integer(),
                  gene2=character(),
                  chunk2=character(),
                  pos2=integer(),
                  xfeature=character())
  
  for (n in 1:nrow(df)) {
    
    temp<-df$sig_chunk[n]
    
    temp<-unlist(strsplit(temp,';',fixed=T))
    
    temp<-gsub("\\(.*$", "", temp, perl = TRUE)
    
    
    for (m in 1:length(temp)){
      
      tempg<-temp[m]
      
      tidx<-which(allchunks$xak_chunk==tempg)
      
      if(length(tidx)==1) {
        tempdf<-data.frame(gene1=gsub("\\_.*$", "", df$xak_chunk[n], perl = TRUE),
                           chunk1=df$xak_chunk[n],
                           pos1=df$pos[n],
                           gene2=gsub("\\_.*$", "", tempg, perl = TRUE),
                           chunk2=tempg,
                           pos2=allchunks$pos[tidx],
                           xfeature=df$xfeature[n])
        
        tdf<-rbind(tdf,tempdf)
        
      } else {
        print(n)
        print(tempg)
      }
      
    }
    
    
  }
  
  return(tdf)
  
}


a_xk<-shared_chunks(a_share)

a_x<-a_xk[which(a_xk$gene2=='mXist'),]
a_k<-a_xk[which(a_xk$gene2=='kcnq1ot1'),]

a_x_repeats<-a_x[which(a_x$xfeature=='repeats'),]
a_k_repeats<-a_k[which(a_k$xfeature=='repeats'),]

cor.test(a_x_repeats$pos1, a_x_repeats$pos2, method = "kendall",alternative = "greater")
cor.test(a_k_repeats$pos1, a_k_repeats$pos2, method = "kendall",alternative = "greater")

cor.test(a_x$pos1, a_x$pos2, method = "kendall",alternative = "greater")
cor.test(a_k$pos1, a_k$pos2, method = "kendall",alternative = "greater")



p<-ggplot(data=a_x,aes(x=pos1,y=pos2,group=xfeature))+
  geom_point(aes(color=xfeature,fill=xfeature), alpha=0.6,size=6) +
  labs(x = "Airn",y = "Xist") +
  #coord_cartesian(xlim=c(0,500000))+
  scale_y_continuous(breaks=seq(from=0, to=31,by=10))+
  scale_x_continuous(breaks=seq(from=0, to=180,by=40))+
  scale_color_manual(values=c('#762a83','#1b7837'))+
  scale_fill_manual(values=c('#af8dc3','#7fbf7b'))+
  theme(panel.background=element_rect(fill='white'),
        plot.margin = margin(10, 15, 10, 10, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x=element_text(size=28),
        axis.text.y=element_text(size=28))

ggsave("XAK_shared_AX_structure_xfeature.svg", plot = p, width = 5.5, height = 5.5, units = "in")


p<-ggplot(data=a_k,aes(x=pos1,y=pos2,group=xfeature))+
  #geom_point(fill='#7fbf7b',color='#1b7837', alpha=0.6,size=4) +
  geom_point(aes(color=xfeature,fill=xfeature), alpha=0.6,size=6) +
  labs(x = "Airn",y = "Kcnq1ot1") +
  #coord_cartesian(xlim=c(0,500000))+
  scale_x_continuous(breaks=seq(from=0, to=180,by=40))+
  scale_y_continuous(breaks=seq(from=0, to=190,by=40))+
  scale_color_manual(values=c('#762a83','#1b7837'))+
  scale_fill_manual(values=c('#af8dc3','#7fbf7b'))+
  theme(panel.background=element_rect(fill='white'),
        plot.margin = margin(10, 15, 10, 10, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x=element_text(size=28, angle=45, hjust=1),
        axis.text.y=element_text(size=28))

ggsave("XAK_shared_AK_structure_xfeature.svg", plot = p, width = 5.5, height = 5.5, units = "in")


############################################################

#############################################
# plot XK chunks, no need to be shared with A for Fig S3 F

x<-read.csv('mXist_XAK_q975chunk.csv',header=T)
k<-read.csv('kcnq1ot1_XAK_q975chunk.csv',header=T)

# add seq number
x$pos<-c(1:nrow(x))
k$pos<-c(1:nrow(k))

colnames(x)[1]<-'xak_chunk'
colnames(k)[1]<-'xak_chunk'


x_xk<-x[(grepl('kcnq1ot1',x$sig_chunk)),]
k_xk<-k[(grepl('mXist',k$sig_chunk)),]

x_xk$xist_seq<-NULL
x_xk$sig_seq<-NULL

k_xk$k_seq<-NULL
k_xk$sig_seq<-NULL

allchunks<-rbind(x[,c(1,5)],k[,c(1,5)])

# add xist feature color to denote repeat or interval
xistfeatures<-c('mXist_repeatA','mXist_ss234','mXist_repeatFdwn','mXist_repeatB',
                'mXist_repeatC','mXist_repeatEbroad')

x_xk$xfeature<-''

for (n in 1:nrow(x_xk)) {
  if (x_xk$xak_chunk[n] %in% xistfeatures) {
    x_xk$xfeature[n]<-'repeats'
  } else {
    x_xk$xfeature[n]<-'intervals'
  }
}


shared_chunks<-function(df) {
  
  tdf<-data.frame(gene1=character(),
                  chunk1=character(),
                  pos1=integer(),
                  gene2=character(),
                  chunk2=character(),
                  pos2=integer(),
                  xfeature=character())
  
  for (n in 1:nrow(df)) {
    
    temp<-df$sig_chunk[n]
    
    temp<-unlist(strsplit(temp,';',fixed=T))
    
    temp<-gsub("\\(.*$", "", temp, perl = TRUE)
    
    
    for (m in 1:length(temp)){
      
      tempg<-temp[m]
      
      tidx<-which(allchunks$xak_chunk==tempg)
      
      if(length(tidx)==1) {
        tempdf<-data.frame(gene1=gsub("\\_.*$", "", df$xak_chunk[n], perl = TRUE),
                           chunk1=df$xak_chunk[n],
                           pos1=df$pos[n],
                           gene2=gsub("\\_.*$", "", tempg, perl = TRUE),
                           chunk2=tempg,
                           pos2=allchunks$pos[tidx],
                           xfeature=df$xfeature[n])
        
        tdf<-rbind(tdf,tempdf)
        
      } else {
        print(n)
        print(tempg)
      }
      
    }
    
    
  }
  
  return(tdf)
  
}


x_xk_plot<-shared_chunks(x_xk)
x_xk_plot<-x_xk_plot[which(x_xk_plot$gene2=='kcnq1ot1'),]

x_xk_repeats<-x_xk_plot[which(x_xk_plot$xfeature=='repeats'),]

cor.test(x_xk_plot$pos1, x_xk_plot$pos2, method = "kendall",alternative = "greater")

cor.test(x_xk_repeats$pos1, x_xk_repeats$pos2, method = "kendall",alternative = "greater")

p<-ggplot(data=x_xk_plot,aes(x=pos1,y=pos2,group=xfeature))+
  #geom_point(fill='#7fbf7b',color='#1b7837', alpha=0.6,size=4) +
  geom_point(aes(fill=xfeature,color=xfeature), alpha=0.6,size=6) +
  labs(x = "Xist",y = "Kcnq1ot1") +
  #coord_cartesian(xlim=c(0,500000))+
  scale_x_continuous(breaks=seq(from=0, to=31,by=10))+
  scale_y_continuous(breaks=seq(from=0, to=190,by=40))+
  scale_color_manual(values=c('#762a83','#1b7837'))+
  scale_fill_manual(values=c('#af8dc3','#7fbf7b'))+
  theme(panel.background=element_rect(fill='white'),
        plot.margin = margin(10, 15, 10, 10, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x=element_text(size=28),
        axis.text.y=element_text(size=28))

ggsave("XK_shared_XK_structure_xfeature.svg", plot = p, width = 5.5, height = 5.5, units = "in")






