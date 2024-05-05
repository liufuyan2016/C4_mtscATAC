####Quality check. -M is marker genes for plot and -F is fragment file from step1.

Rscript ArchR2_Quality.r  -M spleen.txt -H yes -m 2000 -F 01.P_mt -S Spleen -P mm10 -R 0.2 -O 01.P_mt

###merge samples. The Data is all the results directory from the step1 and sample list(see example) is sample information and marker_spleen.txt are genes used for annotation.

Rscript signac_merge_tsne.R -R Data -s sample.list -O Spleen_merge -M marker_spleen.txt

######the below scripts for cells annotation #########
##scCATCH software 
├── scCATCH_cell_anno_cus_marker.BM.r
├── scCATCH_cell_anno_cus_marker.Spleen.r
##scMCA software
├── scMCA_cell_anno.BM.r
├── scMCA_cell_anno.Spleen.r
###singleR software
├── SingleR_cell_anno.BM.r
├── SingleR_cell_anno_Spleen.r
