argv <- commandArgs(T)
library("reshape2")
library("ggplot2")
library("ggpubr")
data <- read.table(argv[1],header=T,check.name=F)    # sample.BGCs.TPM.FDR0.05.group.xls
data1 <- melt(data)
data1$Group <- factor(data1$Group,levels=c("C","E","Y"))
#data1$Group <- factor(data1$Group,levels=c("centenarian","elderly","young"))
c1<-"#00BA38"
c2<-"#F8766D"
c3<-"#619CFF"
color=c(c1,c2,c3) 
pdf(argv[2],width=4,height=4)     # BGCs.TPM.FDR0.05.group.boxplot.pdf
ggplot(data1, aes(x=variable, y =log10(value+0.1))) + 
                geom_boxplot(aes(color = factor(Group)), width=0.4, outlier.size = 0.5, fatten=1,lwd=0.6,    
                position=position_dodge(width=0.8))  + # geom_jitter(aes(col=Group),size=0.5,position = position_dodge(width = 0.6)) + 
                scale_fill_manual(values = color) +
                scale_color_manual(values= color) + 
                ylim(-1,5) +
                stat_compare_means(aes(group=Group),label = "p.signif",size=4,hide.ns=FALSE, label.y = 4) +
                labs(y=expression(TPM(log["10"])),x="BGC",fill="",color="") + 
                theme(axis.text.x = element_text(colour=c(c3,c3,c1,c1,c1,c1),
                                size=8, vjust=0.92, hjust=0.90, angle=45),  
                                axis.text.y = element_text(colour="black",size=10,hjust=1),
                                axis.line = element_line(color="black",size=0.5),
                                axis.line.y = element_line(color="black", size = 0.5),
                                axis.ticks = element_line(color="black",size=0.5),
                                legend.position=c(1,0.85),
                                legend.justification=c(1,0),
                                legend.key = element_blank(),
                                legend.text = element_text(size=10),
                                legend.background = element_blank(),
                                legend.key.width = unit(0.15, "in"),
                                legend.key.height = unit(0.15, "in"),
                                panel.background = element_blank(),
                                plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "in")
                         )
dev.off()
