#!/opt/bin/python2
# -*- coding: utf-8 -*-
'''
从全部位点annovar的注释结果*multianno.txt 中提取需要位点的注释结果并对注释类型进行统计
输入文件说明：
multianno:全部位点annovar的注释结果*multianno.txt
map_file: 需提取位点的染色体和物理位置的文件。两列，无表头，chr\tpos 
'''
import os
import sys
import configparser
def extract_anno(multianno=None,chr_column = 0,pos_column = 1,map_file=None,outdir=None):
	with open(multianno,"r") as multi,\
		 open(map_file,"r") as map,\
		 open(os.path.join(outdir,"anno.txt"),"w") as out_file:
		extract_snp = []
		for line in map:
			snp_pos = line.strip().split("\t")[0] + "_" + line.strip().split("\t")[1]
			extract_snp.append(snp_pos)
		for anno_line in multi:
			if anno_line.startswith("Chr"):
				out_file.write(anno_line)
			else:
				anno_info = anno_line.strip().split("\t")
				chr_pos = anno_info[chr_column] + "_" + anno_info[pos_column]
				for i in extract_snp:
					if i == chr_pos:
						out_file.write(anno_line)
					else:
						continue

extract_anno(multianno=sys.argv[1],chr_column=4,pos_column = 5,map_file=sys.argv[2],outdir=sys.argv[3])


