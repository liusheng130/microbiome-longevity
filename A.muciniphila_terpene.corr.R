library(ggplot2)
library(ggpubr)

dat<-read.table("terpene_Akk.group.xls",check.name=F,head=T)
pdf("terpene_Akk.group.corr.pdf",width=4.2,height=3.6)
c1<-"#00BA38"
c2<-"#F8766D"
c3<-"#619CFF"
ggplot(dat,aes(log10(Akkermansia_muciniphila+0.001),log10(terpene+1))) + 
geom_point(aes(colour=Group)) +
scale_colour_manual(values=c(c1,c2,c3)) +
geom_smooth(method="lm", se=T) +   
labs(y="terpene",x="Akkermansia muciniphila") +
stat_cor(method = "spearman")  
dev.off()
