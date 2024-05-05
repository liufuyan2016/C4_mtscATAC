###calculate mtDNA ratio and nuc fragments numbers for each cell. -input is the output directory from step1
perl stat_each_cell_Mt_ratio.pl -input 00.Report/ -o 01.Ratio -Mt chrM

###draw top N cells according to nuc fragments 
perl draw_topN_cell.stat.pl  -input  01.Ratio/ -o 02.Pic/ -top 1000

#####make stat for each cell mtDNA contents. -input is from the 01.Ratio results.
perl draw_each_cell.stat_from_signaC.pl  -input SP.MT.NUC.cell.stat.xls  -cell  mm10cell_anno.xls -o 03.stat -Mt chrM

###draw heatmap for aged spleen, young spleen and aged bone marrow
Rscript pheatmap_for_sig.r  --infile sig_heatmap.txt.Aged.SP  --outfile sig_heatmap.txt.Aged.SP.png --show_rownames --show_colnames --height 3000  --width 3000 --color.type 1  --fontsize 12 --fontsize_row 12 --fontsize_col 12  --is.log --legend
Rscript pheatmap_for_sig.r  --infile sig_heatmap.txt.Young.SP  --outfile sig_heatmap.txt.Young.SP.png --show_rownames --show_colnames --height 3000  --width 3000 --color.type 1  --fontsize 12 --fontsize_row 12 --fontsize_col 12  --is.log --legend
Rscript pheatmap_for_sig.r  --infile  sig_heatmap.txt.Aged.BM --outfile sig_heatmap.txt.Aged.BM.png --show_rownames --show_colnames --height 3000  --width 3000 --color.type 1  --fontsize 12 --fontsize_row 12 --fontsize_col 12  --is.log --legend
