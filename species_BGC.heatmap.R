argv <- commandArgs(T)
library(reshape)
library(gplots)
library(ggplot2)

dat<-read.table(argv[1],head=F)   # species_BGC.correl.sig.txt
colnames(dat)<-c("BGC","Bacteria","corr","p-value","occ","Number")

pdf(argv[2])    # species_BGC.correl.sig.heatmap.pdf
ggplot(dat, aes(BGC, Bacteria)) +
  geom_tile(aes(fill = corr),colour = "white") +
  scale_fill_gradient2(name="Value",limits = c(-0.5, 0.7),breaks = c(-0.5, 0, 0.4, 0.8),low = "royalblue", 
                       mid = "white", 
                       high = "red3", 
                       midpoint = 0) +
  theme(axis.text.x = element_text(vjust = 0.92, hjust = 0.92, angle = 45))+
  theme(axis.text.y = element_text(face ="italic"))+
  coord_fixed(ratio=1)+
  theme(axis.text= element_text(size = 8))+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), unit="mm"))+
  labs(x = "BGC", y = "Bacteria", title = "correlation")+ 
  theme(plot.title = element_text(size = 13,hjust = 0.5))+
  theme(legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'cm'))+
  theme(legend.title = element_text(size = 8)) 
dev.off()
