argv <- commandArgs(T)
library("ggplot2")
library("reshape2")
library("RColorBrewer")

df<-read.table("allsample.BGCs.TPM.group.average.xls",header=T)
df<-melt(df)
pdf("BGCs.TPM.group.average.barplot.pdf",width=6,height=6)
df$BGC <- reorder(df$BGC,df$value)
df$BGC <- factor(df$BGC, levels=rev(levels(df$BGC)))
col<-c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', 'floralwhite', 'gray24',brewer.pal(12,"Set3"),'#FF996A','#00468B','darkorchid1','darkkhaki','indianred2','cornsilk1','bisque2','chartreuse3','coral3','darkgoldenrod3')
ggplot(data =df, mapping=aes(x =factor(variable), y =value, fill=BGC)) + geom_bar(stat = 'identity') +
scale_fill_manual(values = col)+
labs(x="Age group",y="TPM") + 
guides(fill = guide_legend(ncol = 2)) 
dev.off()
