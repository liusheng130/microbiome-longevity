argv <- commandArgs(T)
if(length(argv) != 4){ stop(
"Rscript XX [species abundance matrix] [BGCs TPM matrix] [method] [outputPrefix]
NOTICE: the column of the profile must be matched with the row of phenotype
method: method of correlation coefficient including spearman, pearson"
) }
# data prof
dat <- read.table(argv[1],check.names=F,head = T)
# phenotype
phe <- read.table(argv[2],check.names=F,head=T)
method <- as.character(argv[3])
prefix <- argv[4]
dat <- as.matrix(dat)
cn <- colnames(dat)
cn <- gsub("","",cn)
rn <- rownames(phe)

### reordering
pid <- pmatch(cn , rn) 
if(length(pid[is.na(pid)]) != 0) stop("please check the data!")
phe <- phe[pid,]
nr <- nrow(phe)

CORT <- function(x, y){
        v <- vector(length = 2)
        x.c <- x[!is.na(x) & !is.na(y)]
        y.c <- y[!is.na(x) & !is.na(y)]
#       n.y <- length(unique(y.c))
        n.x <- length(x.c)
#       n.y <- length(y.c)
        if(n.x < 3) {
                v<- c(NA,NA)
        }else{
                fit <- cor.test(x, y, method = method)
                v[1] <- fit$est
                v[2] <- fit$p.v
        }
        return(v)
}

colN <- colnames(phe)
for(i in 1:ncol(phe)){
        fn <- paste(prefix, "correl", colN[i], "txt", sep = ".")
        Factor <- phe[,i]       
        if(is.factor(Factor)){
                Le <- levels(Factor)
                tab <- table(Factor)
                if(length(levels(Factor)) == 2 & tab[1] * tab[2] !=0){
                        pv <- apply(dat , 1, function(x){wilcox.test(x ~ Factor)$p.v})
                        mean_name <- paste("mean",Le,sep="_")
                        sd_name <- paste("sd",Le,sep="_")
                        occ_name <- paste("occ",Le,sep="_")
                        Mean <- apply(dat , 1 , function(x){tapply(x , Factor , mean)})
                        Mean <- t(Mean)
                        colnames(Mean) <- mean_name
                        Sd <- apply(dat , 1 , function(x){tapply(x , Factor , sd)})
                        Sd <- t(Sd)
                        colnames(Sd) <- sd_name
                        Occ <- apply(dat , 1 , function(x){tapply(x , Factor , function(z){length(z[z!=0])})})
                        Occ <- t(Occ)
                        colnames(Occ) <- occ_name
                        nLe <- cbind(rep(length(Factor[Factor == Le[1]]) , nrow(dat)) , rep(length(Factor[Factor == Le[2]]) , nrow(dat)))
                        colnames(nLe) <- Le
                        # p-value/mean/sd/occ/sample count
                        m <- cbind(pv, Mean, Sd, Occ, nLe) 
                        write.table(m , fn , sep = "\t" , quote = FALSE)
                }
        } else{
                cat("i: ", i, " ", fn, "\n")
                pv <- apply(dat,1,function(x){CORT(x , Factor)})
                pv <- t(pv)
                occ <- apply(dat , 1 , function(z){length(z[z!=0])})
                sampCount <- rep(ncol(dat) , nrow(dat))
                colnames(pv) <- c("estimates","p-value")
                # estimates/p-value/occ/sample count
                m <- cbind(pv, occ, sampCount)
                write.table(m , fn , sep = "\t" , quote = FALSE)
        }
}
