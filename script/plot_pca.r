#!/usr/bin/Rscript

library('getopt')
spec <- matrix(c(
  'verbose','v',2,'integer',
  'help','h',0,'logical',
  'infile','i',1,'character',
  'group','g',1,'character',
  'prefix','p',1,'character',
  'outdir','o',1,'character'
  ),byrow=TRUE,ncol=4)

opt <- getopt(spec,debug=TRUE)
#uasge
print_usage <- function(spec=NULL){
  cat(getopt(spec,usage=TRUE))
  cat("\n")
  cat("
    Usage example:
    # --file: 输入文件,为plink --pca 生成的 *.eigenvec 文件
    # --group: 样本分组文件samplegroup.txt，两列，第一列为样本名称，第二列为品种或分组，带表头，以制表符分隔。
    # 例：Individuals breed
    #     altay.A4J102  Altay
    #     altay.A4J103  Altay
    # --prefix: 输出文件名称前缀
    Rscript plot_pca.r --infile /lustre/project/pinggu/pca/pca.eigenvec --group samplegroup.txt --prefix Ovine_50K --outdir /lustre/project/pinggu/pca/\n")
  q(status=1 )
}

#======================================================================
# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { cat("please give the input file\n");print_usage(spec) }
if ( is.null(opt$group) ) { cat("please give the samplegroup file\n");print_usage(spec)}
if ( is.null(opt$prefix) ) { cat("please give the prefix of output file\n");print_usage(spec)}
if ( is.null(opt$outdir) ) { cat("please give the path to save analysis result\n");print_usage(spec) }
#set some reasonable defaults for the options that are needed
#but were not specified.
#====================================================================================================
opt$infile <- sub("/$","",opt$infile)
opt$group <- sub("/$","",opt$group)
opt$prefix <- sub("/$","",opt$prefix)
opt$outdir <- sub("/$","",opt$outdir)

library('grid')
library("RColorBrewer")
library(ggplot2)
infile <- read.table(opt$infile, header=F,sep=' ')
samplegroup <- read.table(opt$group, header=T,sep="\t")
colnames(samplegroup) <- c("Sample","group")
data <- infile[,-1]  #data:pca.eigenvec去掉第一列
colnames(data) <- c("Sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
data <- merge(samplegroup,data,by="Sample") # 将分组和样本的PC*合并 Sample group PC1 PC2....
		  
### 画PCA 图
# PCA 图theme
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(text=element_text(color='black',size=8),
          rect=element_rect(fill=bg),
          line=element_line(color = 'black',linetype = 1,size = 0.2),
          plot.title=element_text(hjust=0.5),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black', vjust=0.1,size=8),
          axis.line.x = element_line(color = 'black',linetype = 1,size = 0.2),
          axis.line.y = element_line(color = 'black',linetype = 1,size = 0.2),
          axis.text.x = element_text(color='black',size=5),
          axis.text.y = element_text(color='black',size=5),
          axis.ticks.length = unit(0.2,"lines"),    #调节坐标轴刻度
          axis.ticks.x = element_line(color = 'black',linetype = 1,size = 0.2),
          axis.ticks.y = element_line(color = 'black',linetype = 1,size = 0.2),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))
}		  

data$group <- as.factor(data$group)
data$groupnum <- as.numeric(data$group)
shape <- data$groupnum[!duplicated(data$groupnum)]

# PC1-PC2
p <- ggplot(data,aes(x=data$PC1,y=data$PC2)) + 
	 geom_point(aes(shape=data$group,color=data$group),size=2) +
	 scale_shape_manual(values=shape)+
	 labs(title="PCA",x="PC1",y="PC2") + theme_zg()
filename <- paste(opt$outdir,opt$prefix,sep="/")
ggsave(p,file = paste(filename,"pc1-pc2.png",sep="."),width=8, height=5)
ggsave(p,file = paste(filename,"pc1-pc2.pdf",sep="."),width=8, height=5)

# PC1-PC3
p <- ggplot(data,aes(x=data$PC1,y=data$PC3)) + 
	 geom_point(aes(shape=data$group,color=data$group),size=2) +
	 scale_shape_manual(values=shape)+
	 labs(title="PCA",x="PC1",y="PC3") + theme_zg()
ggsave(p,file = paste(filename,"pc1-pc3.png",sep="."),width=8, height=5)
ggsave(p,file = paste(filename,"pc1-pc3.pdf",sep="."),width=8, height=5)

# PC2-PC3
p <- ggplot(data,aes(x=data$PC2,y=data$PC3)) + 
	 geom_point(aes(shape=data$group,color=data$group),size=2) +
	 scale_shape_manual(values=shape)+
	 labs(title="PCA",x="PC2",y="PC3") + theme_zg()
ggsave(p,file = paste(filename,"pc2-pc3.png",sep="."),width=8, height=5)
ggsave(p,file = paste(filename,"pc2-pc3.pdf",sep="."),width=8, height=5)

# 不知为什么在集群上不能保存图片，做不了
# library(scatterplot3d)
# colour_group <-rainbow(length(unique(data$group)))
# colour <-colour_group[as.numeric((factor(data$group)))]
# # filename <- paste(opt$dir,"pca_3d.pdf",sep="/")
# pdf(filename)
# scatterplot3d(data$PC1,data$PC2,data$PC3,pch=20,color=colour,angle=45,main="PCA_3D",cex.symbols=2,mar=c(5.1,4.1,4.1,8.1))
# legend("right",legend=unique(data$group),col=colour_group[as.numeric((unique(data$group)))],pch=20,bg="white",xpd=TRUE,inset=-0.5)
# dev.off()

