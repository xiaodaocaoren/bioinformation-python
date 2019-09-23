source('/lustre/pipeline/DNA_Microarray/chip_site_plot/Utils.R')
data = read_annotation_file( 'data.csv',skip = 8,column_number = 19,chr_col= 10,pos_col= 11, nrow = 159072)

name = paste(1:20)
chrom_length = c( 56831624,48577505,45779781,52389146,42234498,51416486,44630646,47837940,50189764,51566898,34766867,40091314,45874162,49042192,51756343,37887014,41641366,58018742,50746916,47904181   )
bin = 100000


### density
draw_density_plot(data,name,chrom_length,bin,legend.len = 5,max_num_per_bin = 48)


### man
manhattan_plot(data,name,chrom_length,bin,ylim = 50)

## gap
gap(data,name)


## annovar
plot_annovar_stat_from_map(data,"/lustre/project/chips/Illumina/Soybean_200K/company_test/plot/plot/200k/anno/annovar_lib",'Gmax_275_Wm82.a2.v1.gene')
