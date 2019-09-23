library(ggplot2)
library(ggthemes)


read_annotation_file = function(
    snp_file,   #csv文件名称
    column_number, #一共多少列
    chr_col, #染色体是第几列
    pos_col, #物理位置是第几列
    skip = 0, #跳过前多少行
    sep = ",", 
    nrow = -1  #从开始读的位置向后读多少行
){
    # 读取SNP文件，删除NA，删除染色体为0的，得到需要的map 
    type.seq = rep('NULL',column_number)
    # print(type.seq)
    type.seq[chr_col] = 'character'
    type.seq[pos_col] = 'character'
    # print(type.seq)
    map = read.table(snp_file ,colClasses = type.seq,skip = skip,sep = sep,nrow = nrow)  # 列的类型为NULL 和character，读取类型为character的列数据
    if  ( chr_col > pos_col)  map = map[,c(2,1)]

    # 去除 NA
    map[,2] = as.numeric(map[,2])
    map <- map[!is.na(map[, 1]), ]
    map <- map[!is.na(map[, 2]), ]
    # 去除染色体为0的
    map <- map[which(map[,1]!=0),]
    map <- map[order(map[, 1], map[,2]), ]
    return(map)
}

read_map_file = function(
    map_file,  #只包含chr和物理位置的文件,去掉了染色体为0的
    skip = 0,
    sep = ",",
    nrow = -1
){
    map <- read.table(map_file,sep=sep)
    map <- map[order(map[,1],map[,2]), ]
    return(map)
}




draw_bar_plot = function(
    df,
    type = 'type',
    count = 'Count',
    angle = 0,
    file_name = 'bar.png'
){
    bar_plot = ggplot(
        df,aes_string(x=type,y=count)) + geom_bar(stat = 'identity',fill = '#00ffcc'
    )

    bar_plot = bar_plot + theme( 
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
    )

    bar_plot = bar_plot +  geom_text(aes_string(label = count))

    if (angle != 0) {
        print(angle)
        bar_plot = bar_plot  +theme(axis.text.x = element_text(angle = angle, hjust = 0.5, vjust = 0.5)) 
    } # angle指文本倾斜的角度

    ggsave(file_name)

}


draw_density_plot = function(
    map,
    chr_name,
    chrom_length,
    bin,
    legend.len = 15,
    max_num_per_bin = NULL
){
    
    snp_density_figure = function(title, band = 3,width = 5) {
        
        # 画图，标题
        plot(
            NULL, main=title, axes=FALSE,
            xlim=c(0, max(chrom_length) + max(chrom_length)/10), 
            ylim=c(0, length(chrom_length) * band + band), 
            xlab="", ylab="", xaxs="i", yaxs="i"
        )
        
        # 画每条染色体的长度，灰色
        for(i in 1 : length(chrom_length)){
            polygon(
                c(0, 0, chrom_length[i],  chrom_length[i]),
                c(  -width/5    - band * (i - length(chrom_length) - 1),
                    width/5     - band * (i - length(chrom_length) - 1),
                    width/5     - band * (i - length(chrom_length) - 1),
                    -width/5    - band * (i - length(chrom_length) - 1)
                ),
                col="grey", border="grey"
            )
        }
        
        # 画染色体上的分布
        for(i in 1 : length(chr_name)){
            segments(pos.per.chr[[chr_name[i]]], -width/5 - band * (i - length(chrom_length) - 1),
                     pos.per.chr[[chr_name[i]]], width/5 - band * (i - length(chrom_length) - 1),
                     col=col[round(color.index.per_pos [[chr_name[i]]] * length(col) / max.SNP_num.per_bin)], lwd=1
            )
        }
        
        # 染色体名称
        mtext(at=seq(band, length(chrom_length) * band, band),text=rev(chr_name), side=2, las=2, font=1, cex=0.6, line=0.2)
        
        # 染色体长度尺
        axis(
            3, at=seq(0, max(chrom_length), length=10),
            labels=paste(round((seq(0, max(chrom_length), length=10)) / 1e6, 0), "Mb", sep=""),
            cex.axis=0.8, tck=0.01, lwd=2, padj=1.2
        )
        
        # 颜色图例
        legend(
            x=( max(chrom_length) +  max(chrom_length)/100),
            y=( -width/2.5 - band * (length(chrom_length) - length(chrom_length) - 1)),
            legend=legend.label, col=legend.col,pch=15, pt.cex = 3,  bty="n",
            xpd=TRUE , yjust=0, xjust=0
        )
    }
    
    get_top_ten_value_from_matrix = function(num.per.bin,decrease = FALSE){
        index = order(num.per.bin,decreasing = decrease)[1:10]
        col.index = index %/% dim(num.per.bin)[1] + 1
        row.index = index %% dim(num.per.bin)[1]
        col.index[row.index == 0 ] = col.index[row.index == 0 ] -1
        row.index[row.index == 0] = dim(num.per.bin)[1]
        df = data.frame(chr       = colnames(num.per.bin)[ col.index ],
                                  interval  = rownames(num.per.bin)[row.index],
                                  number    = num.per.bin[index]  )
        df
    }
    
    # 读取SNP位置，按染色体拆分，得到每个染色体上的SNP位置
    pos.per.chr =split(map[[2]],map[[1]],drop = TRUE)
    
    
    # 统计每个窗口的SNP位点数，窗口最大位点数    
    color.index.per_pos <- list()
    num.per.bin = list()
    max.SNP_num.per_bin <- 0
    for(i in 1 : length(chrom_length)){
        seq.cut = seq(0,chrom_length[i],bin)
        if (! chrom_length[i] %% bin == 0){seq.cut = c(seq.cut,chrom_length[i])}
        cut.r <- cut(pos.per.chr[[chr_name[i]]],seq.cut) #, labels=FALSE)
        eachbin.num <- table(cut.r)
        eachbin.num[eachbin.num>max_num_per_bin] = max_num_per_bin
        # print(str(eachbin.num))
        num.per.bin[[chr_name[i]]] = eachbin.num
        max.SNP_num.per_bin <- max(max.SNP_num.per_bin, max(eachbin.num))
        color_index = rep(eachbin.num, eachbin.num)

        color.index.per_pos[[chr_name[i]]] <- color_index
    }
    
    bin_interval = c()
    for (i in num.per.bin){
        if (length(i) > length(bin_interval)){
            bin_interval = dimnames(i)$cut.r
        }
    }
    # 输出 SNP number stat，最多，最少SNP的区间
    num.per.bin = do.call(cbind, lapply(lapply(num.per.bin, unlist), `length<-`, max(lengths(num.per.bin))))
    rownames(num.per.bin) = bin_interval

    # print(str(num.per.bin))
    # print(as.numeric(num.per.bin))
    density_summary = summary(as.numeric(num.per.bin))
    stat_name = c('Min','1st Qu','Median','Mean','3rd Qu','Max\n')
    cat( paste(stat_name,sep=","), file='summary.csv' )
    cat( paste( as.numeric(density_summary),sep=","), file='summary.csv',append=TRUE )

    top_ten.df = get_top_ten_value_from_matrix(num.per.bin,TRUE)
    least_ten.df = get_top_ten_value_from_matrix(num.per.bin)
    write.table(top_ten.df,'SNP.num.top10.bin.txt',row.names = FALSE,sep = "\t",quote = FALSE)
    write.table(least_ten.df,'SNP.num.least10.bin.txt',row.names = FALSE,sep = "\t",quote = FALSE)

    
    num.per.bin = as.data.frame(num.per.bin)
    interval = rownames(num.per.bin)
    num.per.bin$Interval = interval
    num.per.bin = num.per.bin[,c(length(num.per.bin),1:length(num.per.bin)-1)]
    write.table(num.per.bin,'SNP_num_per_bin.txt',sep = "\t",quote = FALSE,row.names = FALSE)
    
    # 从0 到最大位点数，生成对应颜色
    col=c("darkgreen", "yellow", "red")
    col=colorRampPalette(col)(max.SNP_num.per_bin)
    
    
    # 创建图例标签，从color中选择图例对应颜色
    if(max.SNP_num.per_bin <= legend.len)	legend.len <- max.SNP_num.per_bin	
    
    legend.label <- round(seq(0, max.SNP_num.per_bin, length=legend.len))
    len <- legend.label[2]
    legend.label <- seq(0, max.SNP_num.per_bin, len)
    
    if(!max.SNP_num.per_bin %in% legend.label){
        legend.label <- c(legend.label, paste(">", max(legend.label), sep=""))
        legend.label.col <- c(legend.label[c(-1, -length(legend.label))], max.SNP_num.per_bin)
    }else{    
        legend.label.col <- c(legend.label[-1])
    }
    
    legend.label.col <- as.numeric(legend.label.col)
    legend.col <- c("grey", col[round(legend.label.col * length(col) / max.SNP_num.per_bin)])
    
    
    
    if (! is.null(dev.list())) dev.off()
    png('SNP_density.png',res = 600,width = 4400,height =3500)
    par(mar = c(5,5,4,2),xpd=TRUE)
    snp_density_figure("SNP Density")
    dev.off()
}


plot_density_manhattan <- function(	map,chr_name,
	col=c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
	xlab="Chromosome",
	ylab="SNP Number",
	ylim = NULL
)
{	
	

	value.pos.list     = split(map[[2]],map[[1]],drop = TRUE)

	chr_num <- length(chr_name)

	band = max(map[[2]]) %/% 5   #横轴两个染色体之间的宽度

	label.position <- c()
	point.pos <- c()
	for(i in 0:(chr_num-1)){
		if (i==0){
			point.pos <- value.pos.list[[chr_name[i+1]]] + band	 #第一条染色体上点的位置+宽度
		}else{
			point.pos <- c(point.pos, max(point.pos) + band + value.pos.list[[chr_name[i+1]]]) #其他染色体上点的位置
		}
		label.position <- c(label.position,max(point.pos)-floor(max(value.pos.list[[chr_name[i+1]]])/2)) #指chr画在横轴上的位置
	}

	# 点的纵坐标
	value=map[,3]
	if (is.null(ylim)) {
	   Max=ceiling(max(value[value!=Inf]))
	} else {
	   Max =ylim
	}
	# print(Max)
	Min=floor(min(value[value!=Inf]))


	# 点的颜色
	pos_num.per.chr <- as.numeric(table(map[,1])[chr_name])  #只保留了每条染色体的频数
	point_color = rep(rep_len(col,chr_num),pos_num.per.chr)

	# # 画图
	png('SNP_density_Manhattan.png', width = 4200,height=1600,res=300)
	par(mar = c(5,6,4,3),xaxs="i",yaxs="r",xpd=TRUE)

	# 散点图
	plot(
		point.pos,value,
		pch=19,cex=0.5,col=point_color,
		xlim=c(0,max(point.pos)),
		ylim=c(Min,Max),
		xlab=xlab,ylab=ylab,
		cex.axis=1,cex.lab=2,font=2,axes=FALSE,main=""
	)


	# 横轴，Chr
	axis(1, at=c(0,label.position),
		cex.axis=1,font=2,labels=c("Chr",chr_name)
	)

	# 纵轴，value
	axis(2,at=seq(Min,(Max),
		ceiling((Max-Min)/5)),cex.axis=1,font=2,
		labels=round(seq(Min,(Max),
		ceiling((Max-Min)/5)), 2)
	)

	dev.off()


}


manhattan_plot  = function(
    map,chr_name,chrom_length,bin,ylim = 30
){
	# 去除 小于0 ， 物理位置大于该染色体最大长度的 行
    map <- map[map[,2] >= 0,]
    for (i in 1:length(chr_name)) {
        bigger_than_chr_length.index = which(map[[1]] == chr_name[i] & map[[2]] > chrom_length[i])
        if (length(bigger_than_chr_length.index) > 0 ) map = map[-bigger_than_chr_length.index,]
    }
    map <- map[order(map[, 1], map[, 2]) ,]
    
    # counting SNP numbers
    pos_per_chr = split(map[[2]],map[1],drop = TRUE)
    
    count_per_chr = list()
    for (i in 1:length(chrom_length))   {
        cut.seq = seq(0,chrom_length[i],bin)
        if (! chrom_length[i]  %% bin == 0) cut.seq = c( cut.seq,chrom_length[i])
		last_bin_length = cut.seq[length(cut.seq)] - cut.seq[length(cut.seq)-1]
        cut.per.chr = cut(pos_per_chr[[chr_name[i]]],cut.seq)
		count = summary(cut.per.chr,maxsum = 1e9)
		# print(count)
		count[length(count)] = count[length(count)] / last_bin_length * bin
		# print(count)
        num_per_chr.df = data.frame(chr =chr_name[i],start = cut.seq[-length(cut.seq)],count = count)
        count_per_chr[[chr_name[i]]] = num_per_chr.df
    }
    
    count.df = do.call(rbind,count_per_chr)
    
    plot_density_manhattan(count.df,chr_name,ylim = ylim)
}


gap = function(
    map,
    chr_name,
    gap_cut_seq= c(0,5e3,10e3,20e3,30e3,40e3,50e3,100e3,500e3,1e8)
    
    # gap_cut_seq = c(0,250,500,750,1000,2000,3000,4e3,5e3,6e3,7e3,8e3,9e3,10e3,11e3,12e3,1e9)
    # gap_cut_seq=  c(0,2e3,4e3,6e3,8e3,10e3,12e3,14e3,16e3,18e3,20e3,25e3,30e3,35e3,40e3,45e3,50e3,1e9)

){

    
    chr_name_list = function(    chr_name,    num){
        chr_name_num = length(chr_name)
        chr_name_name = list()
        for (i in 1:chr_name_num)
        {
            chr_name_name[[chr_name[i]]] = rep(chr_name[i],num[i])
        }
        return(chr_name_name)
    }

        
    gap_merge = function( gap_per_chr,chr_name){
        gap_all = c()
        gap_tmp= gap_all
        for (i in chr_name){
            gap_chr = gap_per_chr[[i]]
            gap_all = c(gap_tmp,gap_chr)
            gap_tmp= gap_all
        }
        return (gap_all)
    }

    num_per_chr = function(    pos_per_chr,    chr_name){
        num = c()
        for (i in 1:length(chr_name)){
            #print(i)
            #print(dim(pos_per_chr[[paste("Chr",i,sep = "")]])[1])
            num[i] = length(pos_per_chr[[chr_name[i]]])    }
        return(num)
    }


    map = unique(map)
    # 统计gap
    print('gap stat')
    pos.per.chr     = split(map[[2]],map[[1]],drop = TRUE) #按染色体分组，得到每个染色体上的SNP物理位置 map[[2]]指map的第二列
    gap.pos.per.chr = lapply(pos.per.chr, function(x){x[-length(x)]}) #每条染色体去掉最后一个snp的位置
    gap.pos.all.chr = gap_merge(gap.pos.per.chr ,chr_name )  #将所有染色体的gap.pos.per.chr合并
    
    gap.per.chr     = lapply(pos.per.chr,function(x){x[-1]-x[-length(x)]}) 
    #列表，每条染色体的gap长度，后一个snp的物理位置减前一个snp的物理位置,读入数据后就已排序
    gap.mean.per.chr = lapply(gap.per.chr, mean)  # 列表，每条染色体的平均gap长度
    gap.mean.per.chr = lapply(gap.mean.per.chr,round,0)  #保留整数
    gap_mean_per_chr = gap_merge(gap.mean.per.chr,chr_name) # 将平均gap按chr_name的顺序排列
    gap_mean_per_chr = data.frame(chr = chr_name,gap_mean_bp = gap_mean_per_chr) 
    gap.all.chr     = gap_merge(gap.per.chr ,chr_name ) # 向量，将所有染色体(chr_name给的)的gap长度放一起
    
    chr_list        = unlist(chr_name_list(chr_name,num_per_chr(pos.per.chr,chr_name)-1),use.names = F)
    gap_df = data.frame(chr = chr_list,pos = gap.pos.all.chr,gap = gap.all.chr)
    
    gap.cut.per.chr = lapply(gap.per.chr,cut,gap_cut_seq) #用cut函数根据分隔点gap_cut_seq把每条染色体的gap划分到各自区间里
    
    gap.cut.all.chr = gap_merge(gap.cut.per.chr ,chr_name ) #将每个gap所属区间进行标记:1,2,3,4,5,6,7,8,9
    stat = table(gap.cut.all.chr)     #统计频次
    gap_length_min = gap_cut_seq[-length(gap_cut_seq)]+1
    gap_length_max = gap_cut_seq[-1]

    interval = paste(gap_length_min,'-',gap_length_max)
    interval = factor(interval,levels = interval)
    gap_stat = data.frame(min = gap_length_min,max = gap_length_max,Interval = interval,Number = as.integer(stat))
    
    print('write table')
    # 输出所有的gap
    write.table(gap_df,"gap_all.txt",sep="\t",row.names =FALSE,quote = FALSE)

    # 输出每条染色体的平均gap
    write.table(gap_mean_per_chr,"gap_mean.txt",sep="\t",row.names = FALSE,quote = FALSE)

    # 输出 前10 的gap
    gap.top10 = head(gap_df[order(gap_df$gap,decreasing = TRUE),],n = 10)
    write.table(gap.top10,'gap_top10.txt',sep = "\t",row.names = FALSE,quote = FALSE)
    
    # 输出gap统计
    write.table(gap_stat,'gap_stat.txt',sep = "\t",row.names = FALSE,quote = FALSE)
    
    #gap summary 统计
    print('gap summary')
    gap_summary = summary(gap_df$gap)
    stat_name = c('Min','1st Qu','Median','Mean','3rd Qu','Max\n')
    cat( paste(stat_name,sep="\t"), file='summary.csv' )
    cat( paste( as.numeric(gap_summary),sep="\t"), file='summary.csv',append=TRUE )



    #做gap条形图
    # gap_stat$group <- as.factor(c(rep("#2473B8",7),rep("#FF9F31",2)))
    bar_plot = ggplot(gap_stat,aes(x = Interval,y = Number)) +
               geom_bar(stat= 'identity',width = 0.4,fill = "#2473B8") +
               # geom_bar(stat= 'identity',width = 0.4,fill = gap_stat$group)  +
               theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
               theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line =element_line(colour = "black"))
    # print(bar_plot)
    ggsave('gap_bar.png',width=8,height=5)
    
    # 做gap密度图
    density_plot=   ggplot(gap_df)+geom_density(aes(x=log10(gap),color = chr,fill = chr))
    # print(density_plot)
    ggsave('gap_density.png',width=8,height=5)
    
}






write_vcf_from_map = function(map,vcf_file_name = 'map.vcf'){
    vcf = map
    vcf$end = map[[2]]
    vcf$ref = '.'
    vcf$alt = '.'
    vcf$qual = '.'
    vcf$filter = '.'
    vcf$info = '.'

    write.table(vcf, file = vcf_file_name, quote = FALSE, sep = "\t",col.names = FALSE,row.names = FALSE)
}

annovar_annotation_from_vcf = function(vcf,annovar_lib_path,buildver){
    annovar_cmd = "perl /opt/annovar/table_annovar.pl %s %s --vcfinput --buildver %s --protocol refGene --operation g --out map.vcf"
    run_cmd = sprintf(annovar_cmd,vcf,annovar_lib_path,buildver)
    print(run_cmd)
    system(run_cmd) # 调用Linux命令，运行run_cmd的命令
}

# plot_annovar_variant_bar_plot= function(file,column_num = 18){
#     col_type = rep("NULL",column_num) 
#     col_type[1] = "character"
#     type = read.table(file,colClasses=col_type,sep = "\t") #读取xxx.variant_function的第一列数据

#     count =  as.data.frame(table(type))

#     write.table(count,'annovar.type.txt',quote = FALSE, sep = "\t",
#         row.names = FALSE)
#     draw_bar_plot (count,type  = 'type',count = 'Freq',angle = 45,
#         file_name ='annovar_bar.png' )
# }

# map.vcf.refGene.variant_function文件的列数可能有问题，导致read.table直接读该文件报错
# 改为先将map.vcf.refGene.variant_function的第一列输出为variant_type,然后用variant_type画图
plot_annovar_variant_bar_plot= function(file){
    cmd = "less %s |awk '{print $1}' > variant_type"
    run_cmd = sprintf(cmd,file)
    system(run_cmd)   #调用Linux命令将map.vcf.refGene.variant_function的第一列输出到type.txt 文件中
    type = read.table("variant_type") #读取type.txt文件
    count =  as.data.frame(table(type))
    write.table(count,'annovar.type.txt',quote = FALSE, sep = "\t",
        row.names = FALSE)
    draw_bar_plot (count,type  = 'type',count = 'Freq',angle = 45,
        file_name ='annovar_bar.png' )
}

plot_annovar_stat_from_map = function(map,annovar_lib,buildver){
   write_vcf_from_map(map)
   annovar_annotation_from_vcf('map.vcf',annovar_lib,buildver)
   plot_annovar_variant_bar_plot('map.vcf.refGene.variant_function')
}

plot_annovar_stat_from_vcf = function(vcf,annovar_lib,buildver){
    annovar_annotation_from_vcf(vcf,annovar_lib,buildver)
    plot_annovar_variant_bar_plot('map.vcf.refGene.variant_function')
}


