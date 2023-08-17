argv <- commandArgs(T)
library(vegan)
library(ggplot2)
#library(ggalt)
dat <- read.table(argv[1],header=T,check.names=F)
group <- read.table(argv[2],sep="\t", header=T)
#dat <- subset(dat,select=as.character(group[,1]))
dat <- dat[rowSums(dat)!=0,]
dat <- as.data.frame(t(dat))

dat_dist <- vegdist(dat, method="bray", binary=F)
dat_pcoa <- cmdscale(dat_dist,k=3,eig=T)
dat_pcoa_points <- as.data.frame(dat_pcoa$points)
sum_eig <- sum(dat_pcoa$eig)
eig_percent <- round(dat_pcoa$eig/sum_eig*100,1)
colnames(dat_pcoa_points) <- paste0("PCoA", 1:3)
dat_pcoa_result <- cbind(dat_pcoa_points, group)
dat.div <- adonis2(dat ~ Group, data = group, permutations = 999, method="bray")
dat_adonis <- paste0("adonis R2: ",round(dat.div$R2,2), "; P-value: ", dat.div$`Pr(>F)`)
write.table(dat_pcoa_result,"pcoa_result.xls",sep="\t",quote=F)

pdf(argv[3],width=4.5,height=4)
c1<-"#00BA38" 
c2<-"#F8766D"
#c3<-"#619CFF"
color<-c(c1,c2)
#color<-c(c1,c2,c3)
ggplot(dat_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Group)) +
    scale_colour_manual(values=c(c1,c2)) +
#   scale_colour_manual(values=c(c1,c2,c3)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title=dat_adonis) +
    theme(axis.text.x = element_text(color='black',size=9),
          axis.text.y = element_text(color='black',size=9)) +
    geom_point(size=2) + 
    stat_ellipse(level=0.9) +
 #  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F) +
    theme_classic()
dev.off()
