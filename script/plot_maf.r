#!/usr/bin/Rscript

library('getopt')
spec <- matrix(c(
  'verbose','v',2,'integer',
  'help','h',0,'logical',
  'infile','i',1,'character',
  'outdir','o',1,'character'
  ),byrow=TRUE,ncol=4)

opt <- getopt(spec,debug=TRUE)
#uasge
print_usage <- function(spec=NULL){
  cat(getopt(spec,usage=TRUE))
  cat("\n")
  cat("
    Usage example:
    #输入文件格式：两列，MAF,Number.\\t分隔
    #例：MAF     Number
    #    0-0.1    1960
    #    0.1-0.2  7606
    Rscript plot_maf.r --infile /lustre/project/pinggu/maf/stat_maf.txt --outdir /lustre/project/pinggu/maf/\n")
  q(status=1 )
}

#======================================================================
# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { cat("please give the input file\n");print_usage(spec) }
if ( is.null(opt$outdir) ) { cat("please give the path to save analysis result\n");print_usage(spec) }
#set some reasonable defaults for the options that are needed
#but were not specified.
#====================================================================================================
opt$infile <- sub("/$","",opt$infile)
opt$outdir <- sub("/$","",opt$outdir)


#评估maf
library(ggplot2)
library(ggthemes)
#library(gcookbook)
#pdf(file="MAF200k.pdf")
data = read.table(file=opt$infile,header = T,sep = "\t")
#windowsFonts(myFont = windowsFont("Times New Roman"))
colnames(data)<-c("MAF","Number")
options(scipen = 200)
p <- ggplot(data, aes(x = MAF, y = Number,fill = MAF))+
  geom_text(aes(label = Number), vjust=-0.2)+
  geom_bar(stat = "identity",width = 0.2)+
  #ggtitle("The distribution of MAF score")+
  scale_fill_wsj("colors6", "")+
  guides(fill = guide_legend(title = NULL,keywidth=0.7,keyheight=0.7))
p + theme(
  axis.ticks.length = unit(0.2,'cm'),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),axis.line = element_line(colour = "black"),
  plot.title = element_text(color = "navy", size = 20, face = "bold",hjust = 0.5),
  axis.title.x = element_text(color = "black", size = 10, face = "bold"),
  axis.title.y = element_text(color = "black",  size = 10, face = "bold"),
  axis.text.x = element_text(size = 12,  color = "black", vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 12,  color = "black", vjust = 0.5, hjust = 0.5),
  legend.text = element_text(size = 10,  colour = "black"))
#dev.off()
ggsave(file=paste(opt$outdir,"maf_bar.png",sep="/"),width=8,height=5)
ggsave(file=paste(opt$outdir,"maf_bar.pdf",sep="/"),width=8,height=5)
