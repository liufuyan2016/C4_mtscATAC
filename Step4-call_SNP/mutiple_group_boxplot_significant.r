library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
Exp_plot <- read.table("all.snp.cell.txt.vaf",header=F,sep="\t")
gene <- c("B cells","T cells","Neutrophils","Monocytes","DC","Eosinophils","Erythroblast","Plasma cells","pDC")
Exp_plot[,3]=Exp_plot[,3]*100;
col <-c("#5CB85C","#337AB7","#F0AD4E")

plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-Exp_plot[which(Exp_plot[,2]==gene[i]),c(3,4)]
  
  colnames(bar_tmp)<-c("Number","sample")#统一命名
  print (bar_tmp)
  my_comparisons1 <- list(c("Aged_spleen", "Aged_bone_marrow")) #设置比较组
  my_comparisons2 <- list(c("Aged_spleen", "Young_spleen"))#设置比较组
  pb1<-ggboxplot(bar_tmp,
                 x="sample",
                 y="Number",
                 color="sample",
                 fill=NULL,
                outlier.shape = NA,
                add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                notch = TRUE,
                 font.label = list(size=30),
                 palette = col)+scale_x_discrete(limits=c("Young_spleen","Aged_spleen","Aged_bone_marrow"))+coord_cartesian(ylim =  c(0, 5))+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")
  pb1<-pb1+stat_compare_means(method="wilcox.test",hide.ns = F,comparisons =c(my_comparisons1,my_comparisons2),label="p.signif",label.y = c(3,4))
  plist2[[i]]<-pb1 
}
pdf("all.comapare.VCF.pdf",width=8,height=6)
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],plist2[[4]],plist2[[5]],plist2[[6]],plist2[[7]],ncol=4)
dev.off()
