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
    Rscript plot_priority.r --infile /lustre/project/pinggu/snp_type/stat_snptype.txt --outdir /lustre/project/pinggu/snp_type/\n")
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

library(ggplot2)
library(ggthemes)
# library(gcookbook)

data = read.table(file=opt$infile,header = T,sep = "\t")
data$SNP_types <- factor(data$SNP_types, levels = data$SNP_types)
#windowsFonts(myFont = windowsFont("Times New Roman"))
colnames(data)<-c("SNP_types","Number")
options(scipen = 200)
p <- ggplot(data, aes(y = Number,x = SNP_types,fill = SNP_types))+
  geom_text(aes(label = Number), vjust = -0.2)+
  geom_bar(stat = "identity",width = 0.25)+
  #ggtitle("The distribution of SNPs type")+
  scale_fill_manual(values = c("#c72e29" ,"#016392" ,"#be9c2e" ,"#098154", "#fb832d", "#000000","#479DA0"))+
  guides(fill = guide_legend(title = NULL,keywidth=1,keyheight=1))
p<- p + theme(
  axis.ticks.length = unit(0.2,'cm'),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),axis.line = element_line(colour = "black"),
  plot.title = element_text(color = "navy",  size = 20, face = "bold",hjust = 0.5),
  axis.title.x = element_text(color = "black",  size = 10, face = "bold"),
  axis.title.y = element_text(color = "black", size = 10, face = "bold"),
  axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 10,  color = "black", vjust = 0.5, hjust = 0.5),
  legend.text = element_text(size = 10,  colour = "black"))

ggsave(p,file=paste(opt$outdir,"SNP_types.png",sep="/"),width=8,height=5)
ggsave(p,file=paste(opt$outdir,"SNP_types.pdf",sep="/"),width=8,height=5)
