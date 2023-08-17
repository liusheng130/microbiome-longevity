#fivefold_cross-validation_random_forest_selection_V2.1

#! /usr/bin/perl
=head1 Description
Select the optimal set for classifier by fivefold cross-validation random  forest selection(code come from paper DOI: 10.1038/ncomms7528)
=head1 Author
Modified by huangyf@genomics.org.cn
=head1 Usage
perl fivefold_cross-validation_random_forest_selection.pl

Required:
-train_p    [s]    train sample set profile
-train_m    [s]    train sample set group infomation
-name       [s]    group names use to calculate, only allow two groups

Optional:
-test_p    [s]    test sample set profile
-test_m    [s]    test sample set sample group infomation
-prefix    [s]    output prefix [default All]
-outdir    [s]    out put directory [default ./]
-color     [s]    color for group    [default "#4DAF4A,#E41A1C"]
-ci_color  [s]    color for roc ci shape [default "red"]
-fea_num    [i]    feature number to plot [default 10]
-fea_hei    [f]    height of feature importance plot [default 2]
-fea_wid    [f]    width of feature importance plot [default 6]
-version        print version information

=head1 Example

perl fivefold_cross-validation_random_forest_selection.pl -train_p Train_Genus_profile.xls -train_m Train_sample.info -name Control:Case -test_p Test_Genus_profile.xls -test_m Test_sample.info

=head1 Note

1. The order of -name represents the x axis order
2. Please do not use '-' or space in sample name

=cut

use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my ($train_p,$train_m,$name,$test_p,$test_m,$version);
my ($outdir,$prefix,$color,$cicolor,$fea_number,$fea_height,$fea_width) = ("./","All","#4DAF4A,#E41A1C","red",10,2,6);

GetOptions (
    'train_p:s'    =>    \$train_p,
    'train_m:s'    =>    \$train_m,
    'name:s'    =>    \$name,
    'test_p:s'    =>    \$test_p,
    'test_m:s'    =>    \$test_m,
    'prefix:s'    =>    \$prefix,
    'outdir:s'    =>    \$outdir,
    'color:s'    =>    \$color,
    'ci_color:s'=>    \$cicolor,
    'fea_num:i'    =>    \$fea_number,
    'fea_hei:f'    =>    \$fea_height,
    'fea_wid:f'    =>    \$fea_width,
    'version'    =>    \$version
);

&version && exit if (defined $version);
die (`pod2text $0`) unless (defined $train_p && defined $train_m && defined $name);

my $Rpath = "/data/data1/liusheng/software/R-4.2.0";

chomp (my $pwd = `pwd`);
$outdir = ($outdir eq "./" || $outdir eq ".") ? $pwd."/".$outdir : ($outdir =~ /^\//) ? $outdir : $pwd."/".$outdir;
system ("mkdir -p -m 755 $outdir") unless (-d $outdir);

for ($train_p,$train_m)
{
    $_ = ($_ !~ /^\//) ? $pwd."/".$_ : $_;
}

if (defined $test_p && defined $test_m)
{
    for($test_p,$test_m)
    {
        $_ = ($_ !~ /^\//) ? $pwd."/".$_ : $_;
    }
}

my @group_name = split /:/,$name;
my $group = join "-",@group_name;
my @color = split /,/,$color;

&read_file($train_p,$train_m,"Train");
&read_file($test_p,$test_m,"Test") if (defined $test_p && defined $test_m);

sub read_file {
    my ($profile,$phenotype,$dataset) = @_;
    open INFO, "$phenotype";
    open GROUP, ">$outdir/$prefix\_$group\_$dataset\_Sample_group.info";
    print GROUP "SampleID\tstate\n";
    my %map_info;
    while (<INFO>)
    {
        chomp;
        next if (/^#/);
        my @array = split /\t/;
        my $flag = 0;
        for (my $i=1; $i<=$#array; $i++)
        {
            foreach my $group_name (@group_name)
            {
                if ($group_name eq $array[$i])
                {
                    $map_info{$array[0]} = $array[$i];
                    $flag = 1;
                }
            }
        last if ($flag == 1);
        }
    }
    close INFO;

    open PRO, "$profile";
    open OUT, ">$outdir/$prefix\_$group\_$dataset\_Sample_profile.xls";
    chomp (my $head = <PRO>);
    my @head = split /\t/,$head;

    my ($title,$sample_num);
    foreach my $sample (@head)
    {
        if (exists $map_info{$sample})
        {    
            $title .= "\t$sample";
            $sample_num ++;
            print GROUP "$sample\t$map_info{$sample}\n";
        }
    }
    print OUT "ID$title\n";
    close GROUP;

    while (<PRO>)
    {
        chomp;
        next if (/^#/);
        my @array = split /\t/;
        my $print = "";
        my $zero_num = 0;
        for (my $i=1; $i<=$#array; $i++)
        {
            if (exists $map_info{$head[$i]})
            {
                $print .= "\t$array[$i]";
                $zero_num ++ if ($array[$i] == 0);
            }
        }
        next if ($zero_num == $sample_num);
        $array[0] = "Label$array[0]" if ($array[0] =~ /^\d+/);
        print OUT "$array[0]$print\n";
    }
    close PRO;
    close OUT;
}

chomp (my $species_num = (split /\s+/,`wc -l $outdir/$prefix\_$group\_Train_Sample_profile.xls`)[0]);
$species_num = $species_num - 1;

my $Rscript = <<R;
library("randomForest",lib.loc="$Rpath/library/")
library("ggplot2",,lib.loc="$Rpath/library/")
library("grid",lib.loc="$Rpath/library/")
library("gridBase",lib.loc="$Rpath/library/")
library("pROC",lib.loc="$Rpath/library/")
library("gridExtra",lib.loc="$Rpath/library/")
library("dplyr",lib.loc="$Rpath/library/")

dat1 <- read.table("$outdir/$prefix\_$group\_Train_Sample_profile.xls", head=T,sep="\t",row.names=1)
conf1 <- read.table("$outdir/$prefix\_$group\_Train_Sample_group.info",head=T,row.names=1,sep="\t")

outcome = conf1\$state
outcome <- sub("$group_name[0]","0",outcome)
outcome <- sub("$group_name[1]","1",outcome)
outcome <- as.factor(outcome)
X <- as.data.frame(t(dat1))
X\$outcome = outcome

#----- 5*10_crossvalidation -----
pdf("$outdir/$prefix\_$group\_ROC.pdf",width=9,height=2)
#par(mfrow = c(1,3),mai=c(0.54,0.5,0.29,0.32),tcl=-0.3,lwd=1.4,mgp=c(2.3,0.5,0),las=1)

set.seed(999)
source("$Bin/randomforest.crossvalidation.r")
result <- replicate(5, rfcv1(X[,-ncol(X)], X\$outcome, cv.fold=10,step=0.9), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
error.cv.cbm <- cbind(rowMeans(error.cv), error.cv)
cutoff <- min(error.cv.cbm[,1])+sd(error.cv.cbm[,1])
cutoff.num <- nrow(error.cv.cbm[error.cv.cbm[,1]<cutoff,])
optimal.set.feature.num <- as.numeric(rownames(error.cv.cbm[error.cv.cbm[,1]<cutoff,])[cutoff.num])

#---- layout the picture -----
plot.new()
plotlayout <- grid.layout(nrow=1,ncol=5)
vp1 <- viewport(layout.pos.col=1,layout.pos.row=1)
vp2 <- viewport(layout.pos.col=2,layout.pos.row=1)
vp3 <- viewport(layout.pos.col=3,layout.pos.row=1)
vp4 <- viewport(layout.pos.col=4,layout.pos.row=1)
vp5 <- viewport(layout.pos.col=5,layout.pos.row=1)
pushViewport(viewport(layout=plotlayout))

##### plot figure 1 #####
pushViewport(vp1)
par(new=TRUE,fig=gridFIG(),mar=c(2.15,2.1,1.75,0.3),cex=0.8)
matplot(result[[1]]\$n.var, cbind(error.cv, rowMeans(error.cv)), type="l",lwd=c(rep(1.4, ncol(error.cv)),1.4), col=c("gray","gray","gray","gray","gray","black"), lty=1, log="x",xlab="Number of variables", ylab="CV Error",mgp=c(1,0.04,0),xaxt="n",yaxt="n",cex.lab=0.7,cex.axis=0.7)
axis(side=1,tcl=-0.2,lwd.ticks=1.4,mgp=c(1,0.04,0),cex=0.7)
axis(side=2,tcl=-0.2,lwd.ticks=1.4,mgp=c(1,0.07,0),cex=0.7)
abline(v=optimal.set.feature.num,col="$cicolor",lwd=1.4)
box(lwd=1.4)
popViewport()

#----- pick marker by corossvalidation -----
k = 1
b <- matrix(0,ncol=$species_num,nrow=50)    ## ncol is gene, genus or mlg number
for(i in 1:5)
{
    for(j in 1:10)
    {
        b[k,] <- result[[i]]\$res[[j]]
        k = k+1
    }
}

mlg.list <- b[,1:10]
list <- c()
k = 1
for(i in 1:10)
{
    for(j in 1:50)
    {
        list[k] <- mlg.list[j,i]
        k = k+1
    }
}

mlg.sort <- as.matrix(table(list))
mlg.sort <- mlg.sort[rev(order(mlg.sort[,1])),]
pick <- as.numeric(names(head(mlg.sort,optimal.set.feature.num)))
tmp = X[,-ncol(X)]
mlg.pick <- colnames(tmp)[pick]
write.table(mlg.pick,"$outdir/$prefix\_$group.cross_validation_pick.txt",sep="\t",quote=F)

#----- train.set -----
train1 <- X[,c(pick,$species_num+1)] ##
set.seed(999)
train1 <- data.frame(train1)
train1.rf <- randomForest(outcome~., data = train1,importance = TRUE)
train1.pre <- predict(train1.rf,type="prob")
p.train <- train1.pre[,2]
combine <- as.data.frame(cbind(predict.value=as.matrix(p.train)[match(rownames(as.matrix(p.train)),rownames(as.matrix(conf1)))],as.matrix(conf1)))
write.table(combine,"$outdir/$prefix\_$group.cross_validation.marker.predict.in.train.txt",sep="\t",quote=F)

###### plot figure 2 #####
temp.data <- read.table("$outdir/$prefix\_$group.cross_validation.marker.predict.in.train.txt",head=T)
temp.data\$state <- factor(temp.data\$state,levels=c("$group_name[0]","$group_name[1]"))
colors <- c("$color[0]","$color[1]")

pushViewport(vp2)
par(new=TRUE,fig=gridFIG())
P2 <- ggplot(data=temp.data,aes(x=state,y=predict.value)) +
    stat_boxplot(geom = "errorbar",width=0.3,size=0.5,color=colors) +
    geom_boxplot(aes(color=state),lwd=0.5,outlier.shape=1) +
    labs(x='',y='Probability of Disease',fill= '') +
    scale_color_manual(values= colors) +
    scale_y_continuous(limits = c(0, 1),breaks=seq(0, 1, 0.5)) +
    theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color="black",size=8),
        axis.line = element_line(color='black',size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )
print(P2,newpage=FALSE)
popViewport()

#----- ROC in train -----
##### plot figure 3 #####
pushViewport(vp3)
par(new=TRUE,fig=gridFIG())
train.roc <- roc(outcome,p.train,percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(train.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c("sp","se.low","se.median","se.high"))
write.table(sens.ci,"$outdir/$prefix\_$group\_train_ci_result.txt",sep="\t",quote=F,col.names=T,row.names=F)
train.ci <- read.table("$outdir/$prefix\_$group\_train_ci_result.txt",head=T)
P3 <- ggroc(train.roc,color="$cicolor",size=0.5,legacy.axes=TRUE)+
    annotate(geom="text",x=0.2,y=0.15,label=paste("AUC =",round(train.roc\$ci[2],2)),size=2.5,hjust=0) +
    annotate(geom="text",x=0.2,y=0.05,label=paste("95% CI:",round(train.roc\$ci[1],2),"-",round(train.roc\$ci[3],2)),size=2.5,hjust=0) +
    geom_abline(intercept = 0, slope = 1, color = "gray") +
    geom_ribbon(data=train.ci,aes(x=1-sp,ymin=se.low,ymax=se.high,fill="$cicolor"),fill="$cicolor",alpha=0.4) +
    scale_x_continuous(breaks=seq(0, 1, 0.5)) +
    scale_y_continuous(breaks=seq(0, 1, 0.5)) +
    theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color="black",size=8),
        axis.line = element_line(color='black',size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        legend.position = "none",
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )
print (P3, newpage=FALSE)
popViewport()
R

if (defined $test_p && defined $test_m)
{
    $Rscript .= <<R;

#-----test.set -----
dat2 <- read.table("$outdir/$prefix\_$group\_Test_Sample_profile.xls",header=T,sep="\t",row.names=1) ##
conf2 <- read.table("$outdir/$prefix\_$group\_Test_Sample_group.info",header=T,sep="\t",row.names=1)
dat2 <- data.frame(t(dat2))

set.seed(999)
p.test<-predict(train1.rf, dat2,type="prob")
pre.test <- as.data.frame(cbind(pre.value=as.matrix(p.test[,2])[match(rownames(as.matrix(p.test)),rownames(as.matrix(conf2)))],as.matrix(conf2)))
write.table(pre.test,"$outdir/$prefix\_$group.cross_validation.marker.predict.in.test.txt",sep="\t",quote=F)

##### plot figure 4 #####
pre.data <- read.table("$outdir/$prefix\_$group.cross_validation.marker.predict.in.test.txt",head=T)
rownames(pre.data) <- NULL
pre.data\$state <- factor(pre.data\$state,levels=c("$group_name[0]","$group_name[1]"))
colors <- c("$color[0]","$color[1]")
pre.data <- pre.data[order(pre.data\$pre.value),]

pushViewport(vp4)
number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number
par(new=TRUE,fig=gridFIG())
P4 <- ggplot(pre.data) +
    geom_point(aes(as.numeric(reorder(as.numeric(row.names(pre.data)),pre.value)),pre.value,color=state),size=0.4) +
    geom_hline(yintercept=0.5,linetype="dashed") +
    labs(x="Samples",y="Probability of Disease",color="") +
    scale_color_manual(values= colors) +
    scale_y_continuous(limits = c(0, 1),breaks=seq(0,1,0.5)) +
    scale_x_continuous(breaks=number_ticks(5)) +
    theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',angle=90,hjust=0.5,size=8),
        axis.title = element_text(color="black",size=8),
        axis.line = element_line(color='black',size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.text = element_text(size=8),
        legend.key = element_blank(),
        legend.background = element_blank(),
#        legend.key.width = unit(0.2, "in"),
#        legend.key.height = unit(0.2, "in"),
        legend.key.size = unit(0.1,"in"),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )
print (P4,newpage=FALSE)
popViewport()

#----- test.ROC -----
outcome.test = conf2\$state
outcome.test <- sub("$group_name[0]","0",outcome.test)
outcome.test <- sub("$group_name[1]","1",outcome.test)
##### plot figure 5 #####
pushViewport(vp5)
par(new=TRUE,fig=gridFIG())
test.roc <- roc(outcome.test,p.test[,2],percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(test.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c("sp","se.low","se.median","se.high"))
write.table(sens.ci,"$outdir/$prefix\_$group\_test_ci_result.txt",sep="\t",quote=F,col.names=T,row.names=F)
test.ci <- read.table("$outdir/$prefix\_$group\_test_ci_result.txt",head=T)
P5 <- ggroc(test.roc,color="$cicolor",size=0.5,legacy.axes=TRUE)+
    annotate(geom="text",x=0.2,y=0.15,label=paste("AUC =",round(test.roc\$ci[2],2)),size=2.5,hjust=0) +
    annotate(geom="text",x=0.2,y=0.05,label=paste("95% CI:",round(test.roc\$ci[1],2),"-",round(test.roc\$ci[3],2)),size=2.5,hjust=0) +
    geom_abline(intercept = 0, slope = 1, color = "gray") +
    geom_ribbon(data=test.ci,aes(x=1-sp,ymin=se.low,ymax=se.high),fill="$cicolor",alpha=0.4) +
    scale_x_continuous(breaks=seq(0, 1, 0.5)) +
    scale_y_continuous(breaks=seq(0, 1, 0.5)) +
    theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color="black",size=8),
        axis.line = element_line(color='black',size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        legend.position = "none",
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )
print (P5,newpage=FALSE)
popViewport()
R
} else {
    print "Warning: No test set result without test data!\n\n"
}

$Rscript .= <<R;
#----- plot importance of each MGS -----
all.importance <- randomForest(outcome~., data = X,importance = TRUE)
write.table(importance(all.importance),"$outdir/$prefix\_$group.feature.importance.xls",sep="\t",quote=F)
imp.data <- read.table("$outdir/$prefix\_$group.feature.importance.xls",header=T)
imp.data <-setNames(cbind(rownames(imp.data),imp.data,row.names=NULL),c("Variable","DecreaseAccuracy","DecreaseGini","MeanDecreaseAccuracy","MeanDecreaseGini"))
Accuracy.data <- arrange(imp.data,desc(MeanDecreaseAccuracy))
Accuracy.data <- head(Accuracy.data,n=$fea_number)
number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number

pdf ("$outdir/$prefix\_$group.feature.importance.pdf",width=$fea_width,height=$fea_height)
P6 <- ggplot(Accuracy.data,aes(x=reorder(Variable,MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity",width=0.7) +
    labs(y="Mean Decrease Accuracy",x="",fill="",color="") +
    scale_y_continuous(breaks=number_ticks(3)) +
    coord_flip() +
    theme(axis.text.x=element_text(color="black",size=8),
        axis.text.y=element_text(color="black",size=7),
        axis.title.x=element_text(color="black",size=8),
        axis.line = element_line(color="black",size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )

Gini.data <- arrange(imp.data,desc(MeanDecreaseGini))
Gini.data <- head(Gini.data,n=$fea_number)
P7 <- ggplot(Gini.data,aes(x=reorder(Variable,MeanDecreaseGini),y=MeanDecreaseGini)) +
    geom_bar(stat = "identity",width=0.7) +
    labs(y="Mean Decrease Gini",x="",fill="",color="") +
    scale_y_continuous(breaks=number_ticks(3)) +
    coord_flip() +
    theme(axis.text.x=element_text(color="black",size=8),
        axis.text.y=element_text(color="black",size=7),
        axis.title.x=element_text(color="black",size=8),
        axis.line = element_line(color="black",size=0.5),
        axis.ticks = element_line(color="black",size=0.5),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), "in"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_blank(),
        panel.background = element_blank()
    )

grid.arrange(P6,P7,nrow=1)
R

open R, ">$outdir/$prefix\_$group\_ROC.R";
print R "$Rscript";
close R;
`$Rpath/bin/R -f $outdir/$prefix\_$group\_ROC.R`;
`rm $outdir/$prefix\_$group\_train_ci_result.txt` if (-e "$outdir/$prefix\_$group\_train_ci_result.txt");
`rm $outdir/$prefix\_$group\_test_ci_result.txt` if (-e "$outdir/$prefix\_$group\_test_ci_result.txt");

