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
    # 输入文件格式：每条染色体每个优先级的个数，以\\t 分隔
    # 输入文件第一列列名必须是Chrom
    #例: Chrom Tiling_order1  Tiling_order2  Tiling_order3  Tiling_order4
    #    chr1  2370  224 2117  201
    #    chr2  1957  251 1946  186
    #    chr3  2387  255 1548  184
    Rscript plot_priority.r --infile /lustre/project/pinggu/priority/stat_priority.txt --outdir /lustre/project/pinggu/priority/\n")
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


## SNPs loci on Chromosomes
library(reshape2)
library(ggplot2)
library(ggthemes)
## library(showtext)
## CairoPDF("SNPs_test.pdf")
options(scipen = 200)
## showtext.auto(enable = TRUE)
## png(file="SNPs_test.png",width=800,height=700,res=100)
data = read.table(file=opt$infile,header = T)
mydata <- melt(data,id.vars = "Chrom",variable.name = "Tiling_order",value.name = "Number")
mydata$Chrom <- as.factor(mydata$Chrom)
# windowsFonts(myFont = windowsFont("Times New Roman"))
p <- ggplot(mydata,aes(Chrom,Number,fill = Tiling_order))+
  geom_bar(stat = "identity",position = "stack", width = 0.3)+
  scale_fill_wsj("colors6", "")+
  guides(fill = guide_legend(title = NULL,keywidth=0.7,keyheight=0.7))
p + theme(
  axis.ticks.length = unit(0.1,'cm'),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),axis.line = element_line(colour = "black"),
  # plot.title = element_text(color = "navy", family = "myFont", size = 20, face = "bold.italic",hjust = 0.5),
  # axis.title.x = element_text(color = "black", family = "myFont", size = 14, face = "bold"),
  # axis.title.y = element_text(color = "black", family = "myFont", size = 14, face = "bold"),
  # axis.text.x = element_text(angle = 60,size = 12, family = "myFont", color = "black", vjust = 0.5, hjust = 0.5),
  # axis.text.y = element_text(size = 12, family = "myFont", color = "black", vjust = 0.5, hjust = 0.5),
  # legend.text = element_text(size = 10, family = "myFont", colour = "black"))

  plot.title = element_text(color = "navy", size = 20, face = "bold.italic",hjust = 0.5),
  axis.title.x = element_text(color = "black",  size = 10, face = "bold"),
  axis.title.y = element_text(color = "black",  size = 10, face = "bold"),
  axis.text.x = element_text(size = 10,  color = "black", vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 10,  color = "black", vjust = 0.5, hjust = 0.5),
  legend.text = element_text(size = 10,  colour = "black"))

ggsave(file=paste(opt$outdir,"priority_bar.png",sep="/"),width=8,height=5)
ggsave(file=paste(opt$outdir,"priority_bar.pdf",sep="/"),width=8,height=5)
