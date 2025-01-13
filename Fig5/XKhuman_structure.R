# plot to check whether there are structure conservation among XK
# separately label chunks that are related to a Xist repeat or interval
# do the correlation with the repeats only chunks and all chunks
####################################
setwd("/to/the/folder/contains/XIST_OT1_q975chunk.csv")

library(ggplot2)

xk<-read.csv('XIST_OT1_q975chunk.csv',header=T)

# add xist feature color to denote repeat or interval
xistfeatures<-c('rA_316_743','Fbroad_743_1807','rB1_1951_2041','rB2_2781_2919',
                'rD_6113_8140','rE_11922_12586')

xk$xfeature<-''

for (n in 1:nrow(xk)) {
  if (xk$chunk1[n] %in% xistfeatures) {
    xk$xfeature[n]<-'repeats'
  } else {
    xk$xfeature[n]<-'intervals'
  }
}

xk_repeat<-xk[which(xk$xfeature=='repeats'),]

cor.test(xk_repeat$pos1, xk_repeat$pos2, method = "kendall",alternative = "greater")

cor.test(xk$pos1, xk$pos2, method = "kendall",alternative = "greater")

p<-ggplot(data=xk,aes(x=pos1,y=pos2,group=xfeature))+
  #geom_point(fill='#7fbf7b',color='#1b7837', alpha=0.6,size=4) +
  geom_point(aes(fill=xfeature,color=xfeature), alpha=0.6,size=6) +
  labs(x = "XIST",y = "KCNQ1OT1") +
  #coord_cartesian(xlim=c(0,500000))+
  scale_x_continuous(breaks=seq(from=0, to=36,by=10))+
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

ggsave("human_XK_structure_xfeature.svg", plot = p, width = 5.5, height = 5.5, units = "in")
