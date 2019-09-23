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
    Rscript plot_anno_type.r --infile /lustre/project/pinggu/anno/stat_annotype.txt --outdir /lustre/project/pinggu/anno/\n")
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



library(grid)
library(ggplot2)
library(ggthemes)
#pdf(file="plot.Functional_sitesCattle600K.pdf")
# library(gcookbook)
data = read.table(opt$infile,header = T,sep = "\t")
data$Type <- factor(data$Type, levels = data$Type)
#windowsFonts(myFont = windowsFont("Times New Roman"))
colnames(data)<-c("Type","Number")
options(scipen = 200)
p <- ggplot(data, aes(x = Type, y = Number,fill = Type))+
  geom_text(aes(label = Number), vjust=-0.2,size=2)+   # 使条形图在图上显示数字
  geom_bar(stat = "identity",position=position_dodge(1),width = 0.7)+
  #ggtitle("The distribution of functional sites")+ ## 后续的评估报告作图就不加标题了
  scale_fill_manual(values = c("#c72e29" ,"#016392" ,"#be9c2e" ,"#098154", "#fb832d", "#000000","#479DA0","#ECAF5E","#7F6CB9","#177778","#2473B8","#FF9F31","#003366"))+
  guides(fill = guide_legend(title = NULL,keywidth=0.7,keyheight=0.7))
p <- p + theme(
  # legend.key.size=unit(1,'cm'),
  axis.ticks.length = unit(0.1,'cm'),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),axis.line = element_line(colour = "black"),
  plot.title = element_text(color = "navy",  size = 20, face = "bold",hjust = 0.5),
  axis.title.x = element_text(color = "black",  size = 10, face = "bold"),
  axis.title.y = element_text(color = "black",  size = 10, face = "bold"),
  axis.text.x = element_text(angle = 70, size = 8,  color = "black", vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(size = 8,  color = "black", vjust = 0.5, hjust = 0.5),
  legend.text = element_text(size = 8,  colour = "black"))

ggsave(p,file=paste(opt$outdir,"anno_type.png",sep="/"),width=8,height=5)
ggsave(p,file=paste(opt$outdir,"anno_type.pdf",sep="/"),width=8,height=5)