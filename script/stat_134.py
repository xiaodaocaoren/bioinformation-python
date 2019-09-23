#!/bin/python
#!-*-coding:utf-8-*-
import os,sys
import collections
import pandas as pd
def stat_priority(infile,chr_column=4,priority_column=13,sep = "\t",outdir="priority"):
	'''	
	统计每条染色体各个优先级的数量
	'''
	df = pd.read_table(infile,sep=sep)
	max_pri = max(df.iloc[:,priority_column-1]) # 优先级那列最大的是几
	print(max_pri)

	with open(infile,"r") as in_file,\
		 open(os.path.join(outdir,"stat_all_134.txt"),"w") as stat:
		chr_pri_dict = {}
		chr_pri_dict = collections.OrderedDict() #将字典变有序，按插入顺序输出
		for line in in_file:
			if line.startswith("Locus_Name"):
				continue			
			else:
				info = line.strip().split(sep)
				chr = info[int(chr_column)-1]  # 染色体
				priority = info[int(priority_column)-1]  # 优先级
				if chr not in chr_pri_dict:
					order_dict = {}
					order_dict = collections.OrderedDict()
					for i in range(1,max_pri + 1):
						order_dict[str(i)] = 0
					order_dict[priority] += 1
					chr_pri_dict[chr] =  order_dict
				else:
					order_dict[priority] += 1
					chr_pri_dict[chr] =  order_dict
				# print(chr_pri_dict)

		order = [str(i) for i in range(1,max_pri+1)]	
		stat.write("\t".join(["Chrom"] + ["Tilingorder" + i for i in order]) + "\n")
		for key in chr_pri_dict:
			out_info = []
			for item in order:
				out_info.append(str(chr_pri_dict[key][item]))
			stat.write("\t".join([key,"\t".join(out_info)]) + "\n")

stat_priority("illumina.1x_70k.final.csv.add.Final_Score.add.MAF.add.Tiling_order.novip_fst.order1.3.4",chr_column=4,priority_column=15,sep = ",",outdir="./")
