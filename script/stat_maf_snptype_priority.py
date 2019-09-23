#!/bin/python
#!-*-coding:utf-8-*-
'''
对定制芯片的选点打分文件进行统计。
该文件最多有一行以Locus_Name开头的表头，其余为位点信息行,文件以逗号分隔
统计优先级、maf、snp类型
'''
import sys,os
import collections
import pandas as pd
import argparse

def stat_priority(infile,chr_column=4,priority_column=13,outdir="priority"):
	'''	
	统计每条染色体各个优先级的数量
	'''
	df = pd.read_table(infile,sep=",")
	max_pri = max(df.iloc[:,int(priority_column)-1]) # 优先级那列最大的是几,优先级为1,2,3,4,5等
	# print(max_pri)

	with open(infile,"r") as in_file,\
		 open(os.path.join(outdir,"stat_priority.txt"),"w") as stat:
		chr_pri_dict = {}
		chr_pri_dict = collections.OrderedDict() #将字典变有序，按插入顺序输出
		for line in in_file:
			if line.startswith("Locus_Name"):
				continue			
			else:
				info = line.strip().split(",")
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
		stat.write("\t".join(["Chrom"] + ["Tiling_order" + i for i in order]) + "\n")
		for key in chr_pri_dict:
			out_info = []
			for item in order:
				out_info.append(str(chr_pri_dict[key][item]))
			stat.write("\t".join([key,"\t".join(out_info)]) + "\n")

def stat_maf(infile,maf_column = 16, outdir="maf"):
	'''
	统计所有位点的maf分布情况
	'''
	with open(infile,"r") as in_file,\
		 open(os.path.join(outdir,"stat_maf.txt"),"w") as stat:
		num_1 = num_2 = num_3 = num_4 = num_5 = 0
		for line in in_file:
			if line.startswith("Locus_Name"):
				continue			
			else:
				info = line.strip().split(",")
				maf = float(info[int(maf_column) - 1])  # maf列索引为maf_column - 1
				if maf <= 0.1:
					num_1 += 1
				elif maf > 0.1 and maf <= 0.2:
					num_2 += 1
				elif maf > 0.2 and maf <= 0.3:
					num_3 += 1
				elif maf > 0.3 and maf <= 0.4:
					num_4 += 1
				elif maf > 0.4 and maf <= 0.5:
					num_5 += 1
		stat.write("MAF\tNumber\n")
		stat_info = "0-0.1\t{}\n0.1-0.2\t{}\n0.2-0.3\t{}\n0.3-0.4\t{}\n0.4-0.5\t{}\n".format(str(num_1),str(num_2),str(num_3),str(num_4),str(num_5))
		stat.write(stat_info)

def stat_snptype(infile,seq_column=3, outdir="snp_type"):
	'''
	统计snp为A/T、C/G、A/C.... 等不同类型的数量
	'''
	with open(infile,"r") as in_file,\
		 open(os.path.join(outdir,"stat_snptype.txt"),"w") as stat:
		num_AC = num_AG = num_AT = num_CG = num_CT = num_GT = num_Indel = 0
		for line in in_file:
			if line.startswith("Locus_Name"):
				continue			
			else:
				info = line.strip().split(",")
				snptype = info[int(seq_column)-1].upper()  # 将所有碱基都统一为大写
				if "A/C" in snptype or "C/A" in snptype:
					num_AC += 1
				elif "A/G" in snptype or "G/A" in snptype:
					num_AG += 1
				elif "A/T" in snptype or "T/A" in snptype:
					num_AT += 1
				elif "C/G" in snptype or "G/C" in snptype:
					num_CG += 1
				elif "C/T" in snptype or "T/C" in snptype:
					num_CT += 1
				elif "G/T" in snptype or "T/G" in snptype:
					num_GT += 1
				elif "-/" in snptype or "/-" in snptype:  
					num_Indel += 1
		stat.write("SNP_types\tNumber\n")
		stat_info = "A/C\t{}\nA/G\t{}\nA/T\t{}\nC/G\t{}\nC/T\t{}\nG/T\t{}\nIndel\t{}\n".format(str(num_AC),str(num_AG),str(num_AT),str(num_CG),str(num_CT),str(num_GT),str(num_Indel))
		stat.write(stat_info)

 

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Statistical information on the loci of customized chips")
	parser.add_argument("--infile",required = True, help = "Input file for statistics")
	parser.add_argument("--chr_column",required = True,help = "the column of chr in input file")
	parser.add_argument("--priority_column",required = True,help = "the column of priority in input file")
	parser.add_argument("--maf_column",required = True,help = "the column of maf in input file")
	parser.add_argument("--seq_column",required = True,help = "the column of seq in input file")
	parser.add_argument("--outdir",required = True, help = "Path of output files")

	args = parser.parse_args()
	infile = args.infile
	chr_column = args.chr_column
	priority_column = args.priority_column
	maf_column = args.maf_column
	seq_column = args.seq_column
	outdir = args.outdir


	if not os.path.exists(outdir):
		os.makedir(outdir) 

	[os.makedirs(os.path.join(outdir,i)) for i in ["priority","maf","snp_type"] if not os.path.exists(os.path.join(outdir,i))]
	stat_priority(infile=infile,chr_column=chr_column,priority_column=priority_column,outdir=os.path.join(outdir,"priority"))
	stat_maf(infile=infile,maf_column = maf_column,outdir=os.path.join(outdir,"maf"))
	stat_snptype(infile=infile,seq_column=seq_column,outdir=os.path.join(outdir,"snp_type"))

