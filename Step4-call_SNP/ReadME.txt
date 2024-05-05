## make GATK. The inut is from step 1 results.Each samples must be running below.

perl make_mGATKV4.pl -input P1_1 -o 02.GATK

##combine all samples,for example in aged spleen samples:P1_1 P1_2 P2_1 P2_2

perl combine_each_sample_mgatk.pl  -index mm10 -mgatk P1_1/final/ P1_2/final/ P2_1/final/ P2_2/final/ -o P_O

###Find informative mtDNA variants. The "P_combined.harmony_GeneActivityLogNormalize.rds" was from step 2. This step produced "mm10.mt_SNP.txt" and "mm10cell_anno.xls" files.

Rscript SignaC_mgatk.r  -R  P_new_combined.harmony_GeneActivityLogNormalize.rds  -m P_O/ -O P_SNP/ -S sample.list -I mm10

####SNP sequence quality and frequency  filter

bamToFreq-master/bin/bamToFreq -q 30 P_O.mm10.final.bam  ##use mtDNA alignment bam file and produced P_O.final_freqs.csv file 

perl calculation_Alt_ratio.pl  -index mtscATATC  -i  P_O.final_freqs.csv  -fa  mm10_add.fasta -o P_O.final_freqs.xls

perl draw_each_SNP_frequency.pl  -i mm10.mt_SNP.txt  -sc  P_O.final_freqs.xls   -index mm10  -o combine_barplot_SNP.txt -sp mouse  -len 16299 ##mm10.mt_SNP.txt from SignaC_mgatk.r and -len is the mtDNA length

perl filter_SNP_from_mgatk_pipline.pl -bam  P_O.mm10.final.bam -SNP mm10.mt_SNP.txt  -barplot_SNP combine_barplot_SNP.txt.all  -od Final_SNP  ###"combine_barplot_SNP.txt.all" from above step and produce "final_select.SNP.txt" file

##draw heatmap
#you can order the SNP according to the frequency. The "final_select.SNP.order" is the ordered file.
perl  HeatMap/draw_HeatMap.pl  -MT -i final_select.SNP.order  -all -tpm mm10.mt_SNP.txt -N 1111 -group cell_group.xls.old  -cluster 1 -complex -o SNP_heatmap


######SNP annotation use ############ 

#Install the annovar software and then running build database use mtDNA annotation "mt.gff3"

perl -e' open(IN,"mm10.mt.gff3");
while(<IN>){
  chomp;
  my @tmp=split/\t+/;
  if($tmp[2]=~/exon/){
	my ($c,$xx,$name,$type)=$tmp[8]=~/Parent=([^;]+).*gbkey=(\S+);gene=(\S+);product=(\S+)/;
	print "111\t$c\tchrM\t$tmp[6]\t$tmp[3]\t$tmp[4]\t$tmp[4]\t$$tmp[4]\t1\t$tmp[3],\t$tmp[4],\t0\t$type\tnone\tnone\t-1\n";
}
}'>mm10_MT_ensGene.txt 

perl annovar/retrieve_seq_from_fasta.pl mm10_MT_ensGene.txt  -seqdir annovar/mousedb/mm10_seq/  -format ensGene -outfile annovar/mousedb/mm10_MT_ensGeneMrna.fa

###make annotation

##input format for "Aged_spleen.SNP.list" is below, we can prepare from the "final_select.SNP.txt" file.
MT      5941    5941    G       A
MT      5902    5902    C       T

###running annovar for annotation and produced "Aged_spleen.SNP.list.exonic_variant_function" file 

perl annovar/annotate_variation.pl -buildver mm10_MT Aged_spleen.SNP.list  annovar/mousedb/ --geneanno -dbtype ensGene

###annotation stats 

perl calculate_SNP.num.pl -i  final_select.SNP.txt  -cell mm10cell_anno.xls  -od 02.Num_VAF  -index Spleen_Aged  -top snp.top -final final_select.SNP.txt -ann Aged_spleen.SNP.list.exonic_variant_function ##"final_select.SNP.txt" from above and "snp.top" is you pay attention  or split groups in the next step.

perl  SNP_anno/stat_anno_and_draw.pl -i 02.Num_VAF/Spleen_Aged.snp.cell.xls  -index spleen -od 03.STAT_ann

Rscript SNP_anno/pheatmap.r --infile spleen.SNP_mutation_num.for_heatmap.txt --outfile spleen.SNP_mutation_num.for_heatmap.png --show_rownames --show_colnames t 1500  --width 4000 --color.type 1

###Compare the VAF and SNP numbers between different samples  for each cell type
#The input file "*.each_cell_snp.num.draw" and "*.snp.cell.txt.vaf.draw" from the output of  "calculate_SNP.num.pl"
Rscript mutiple_group_boxplot_significant.r ##For common cell types between young spleen, aged spleen and aged bone marrow 
Rscript two_group_boxplot_significant.r ##For common cell types between young spleen and aged spleen.

###Split groups Mutation and WT clones by final three SNP. And then find differentially accessible peaks between clones.The input is "mm10.mt.RDS" and the three SNPs selected.

Rscript make_group_according_SNP.r





