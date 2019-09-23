#!/opt/bin/python2
# -*- coding: utf-8 -*-

import os
import sys
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])


output_dir      = config['Default']['output_dir']
bin             = config['Default']['bin']

snp_file        = config['File']['snp_file']
sep             = config['File']['sep']
skip            = config['File']['skip']
column_number   = config['File']['column_number']
chr_col         = config['File']['chr_col']
pos_col         = config['File']['pos_col']
nrow            = config['File']['nrow']

chr_name        = config['Chr']['chr_name']		#染色体名称
chrom_length    = config['Chr']['chrom_length'] #每条染色体的长度

annovar_lib     = config['Annovar']['annovar_lib']  #annovar数据库路径
buildver        = config['Annovar']['buildver'] #数据库前缀


# cmd = '''
#     source('/lustre/project/test/lijia/script/chip_site_plot/Utils.R')


#     snp_file        = "{snp_file}"
#     sep             = '{sep}'
#     column_number   = {column_number}
#     chr_col         = {chr_col}
#     pos_col         = {pos_col}
#     skip            = {skip}
#     nrow            = {nrow}

#     chr_name = c( {chr_name}) 
#     chrom_length = c( {chrom_length}   )
#     bin = {bin}

#     annovar_lib
#     buildver


#     map = read_annotation_file(snp_file,column_number,chr_col,pos_col,skip = skip,nrow = nrow,sep  = sep)

#     draw_density_plot(map,chr_name ,chrom_length,bin)

#     manhattan_plot(map,chr_name ,chrom_length,bin)

#     gap(map,chr_name)

# '''.format(snp_file = snp_file,sep  =sep, column_number = column_number,chr_col = chr_col , pos_col = pos_col, skip = skip,nrow = nrow,chr_name = chr_name,chrom_length = chrom_length, bin = bin)
# print(cmd)
cmd_common = '''
source('/lustre/project/test/lijia/script/chip_site_plot/Utils.R')

snp_file        = "{snp_file}"
sep             = '{sep}'
column_number   = {column_number}
chr_col         = {chr_col}
pos_col         = {pos_col}
skip            = {skip}
nrow            = {nrow}

chr_name = c({chr_name}) 
chrom_length = c({chrom_length})
bin = {bin}

annovar_lib = "{annovar_lib}"
buildver = "{buildver}"

map = read_annotation_file(snp_file,column_number,chr_col,pos_col,skip = skip,nrow = nrow,sep  = sep)
'''.format(snp_file = snp_file,sep  =sep, column_number = column_number,chr_col = chr_col , pos_col = pos_col, skip = skip,nrow = nrow,chr_name = chr_name,chrom_length = chrom_length, bin = bin, annovar_lib = annovar_lib, buildver = buildver)

cmd_denisity=cmd_common+'\n'+ 'draw_density_plot(map,chr_name ,chrom_length,bin)' +'\n' + 'manhattan_plot(map,chr_name ,chrom_length,bin)'
cmd_gap=cmd_common+'\n' +'gap(map,chr_name)'
cmd_annovar=cmd_common+'\n'+'plot_annovar_stat_from_map(map,annovar_lib,buildver)'
# print(cmd_denisity)
# print(cmd_gap)
# print(cmd_annovar)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
# if not os.path.exists(output_dir+'/scripts'):
#     os.makedirs(output_dir+'/scripts')   

# with open(output_dir+'/scripts/run.r','w') as f:
#     f.write(cmd)
if not os.path.exists(output_dir+'/gap'):
    os.makedirs(output_dir+'/gap')
with open(output_dir+'/gap/plot_gap.r','w') as g:
    g.write(cmd_gap)

if not os.path.exists(output_dir+'/density'):
    os.makedirs(output_dir+'/density')
with open(output_dir+'/density/plot_denisity.r','w') as d:
    d.write(cmd_denisity)

if not os.path.exists(output_dir+'/anno'):
    os.makedirs(output_dir+'/anno')
with open(output_dir+'/anno/plot_anno.r','w') as a:
    a.write(cmd_annovar)


# run_cmd  = """
# cd {}

# Rscript scripts/run.r

# cd - 
# """.format(output_dir)
# #print(run_cmd)
# os.system(run_cmd)
#run_cmd_denisity

run_cmd = "cd {0}/density && Rscript plot_denisity.r && cd {0}/gap && Rscript plot_gap.r && cd {0}/anno && Rscript plot_anno.r && cd {0}".format(output_dir)
os.system(run_cmd)
