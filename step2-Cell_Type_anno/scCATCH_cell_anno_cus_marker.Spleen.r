library(scCATCH)
library(tibble)
library(Seurat)
RDS=readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/kangjingmin/project/01.Mt_sc/10.Article/04.Distribution/01.Old/03.Cluster_SignalC/01.P_TSNE_50/P_combined.harmony_GeneActivityLogNormalize.rds")
DefaultAssay(RDS) <- 'RNA'
obj <- createscCATCH(data = RDS[['RNA']]@data, cluster = as.character(RDS@meta.data$peaks_snn_res.0.4))
custom_marker <- data.frame(species = c("mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse","mouse"),

                            tissue = c(
"ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts","ts"),
                            cancer= c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal"),
                            condition = c("Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell","Normal cell"),
                            subtype1 = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"),
                            subtype2 = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"),
                            subtype3 = c("sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub"),
                            celltype = c("mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","mM12_Macro-Maf","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","Neutrophil","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM13_Macro-Ccl12","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM06_cDC1-Clec9a","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM08_Mono-Ly6c2","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","mM15_Macro-Vegfa","B","B","B","B","B","B","B","B","B","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","mM11_Macro-Mafb","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","M14_Macro-Mgl2","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","Mast","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","mM09_Mono-Nr4a1","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","pDC","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","mM07_cDC1-Ccl22","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","M05_cDC2-Itgax","T","T","T","T","T","T","T","T","T","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","mM10_Mono-Itgal","erythroblast","erythroblast","erythroblast","erythroblast","erythroblast","erythroblast","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a","M04_cDC2-Cd209a"),
gene=c("Maf","C1qa/b/c","Axl","Ccl12","Trem2","Tgfbr1","Cd81","Cd72","Cd63","Abhd12","Adgre1","Ms4a7","Ntpcr","Olfml3","Tmem119","Hpgds","Ets2","Mxd1","Klf2","Egr1","Csrnp1","Csf3r","Cxcr2","Cd24a","C5ar1","Trem1","Il1r2","S100a8","S100a9","Cxcl2","Il1rn","Clec4d","Clec4e","Upp1","Irf7","Mafb","C1qa/b/c","Ccl12","Sdc3","Lgmn","H2-Aa","Batf3","Cbfa2t3","Pa2g4","Xcr1","Itgb7","Arsb","Ckb","Fgd2","Naga","Pak1","Rab7b","Wdfy4","Ppm1m","Cd24a","Flt3","Cd83","Slamf7","Slamf8","Stat2","Mxd1","Bach1","Bcl3","Ly6c2","Mgst1","Ccl9","Smox","Ifit3","Ifi205","Ifit2","Isg20","Plaur","Tlr2","Cd14","Cd300lf","Tgm2","Cxcl10","Fn1","F13a1","Ccl2","Cebpb","Bhlhe40","Atf4","Tgif1","Atf3","Spp1","Vegfa","Mmp12","Adam8","Cd274","Cd63","Thbs1","C3ar1","Il1rn","Clec4d","Emp1","Arg1","Ero1l","Hilpda","Hmox1","Sgk1","Cd79a","Blk","Bcl11a","Ms4a1","Pax5","Cd19","H2-Aa","H2-Eb1","H2-Ab1","Mafb","Cebpb","C1qa/b","Axl","Spint1","C3ar1","F11r","Ccr5","Ccr1","Fcgr1","Fcgr4","Ly6i","Cfb","Mmp14","Clec5a","Ier3","Klf4","Mgl2","Lrp1","Fn1","Pltp","Axl","Ear2","Tnip3","Birc5","Gata2","Nfe2","Klf7","Cpa3","Ms4a2","Ifitm1","Hdc","Slc41a3","Csf1","Sytl3","Hgf","Il6","Il4","Ccl3","Cdh1","Cxcr2","Cd9","Cd69","Il18rap","Nr4a1","Rara","Klf13","Nfil3","Klf4","Bcl3","Cd300a","Il17ra","Tnfrsf21","Tnfsf13","Cd302","Cd14","Ifitm6","Gstm1","Idh1","Gsr","Adgre5","Fn1","F13a1","Anxa1","Ramp1","C3","Tcf4","Runx2","Bcl11a","Tsc22d1","Irf7","Sigech","Ly6d","Cox6a2","Smim5","Klk1","Rpgrip1","Lefty1","Upb1","Ccr9","Cd7","Sell","Cd164","Pltp","P2ry14","Ptprs","Relb","Etv3","Batf3","Aebp2","Nfkb2","Ccl22","Ccl5","Il15","Ccr7","Il15ra","Plxnc1","Prnp","Cd40","Birc2","Fscn1","Anxa3","Cacnb3","Nudt17","Socs2","Tspan3","Serpinb6b","Bcl3","Bhlhe40","Cbfa2t3","Itgax","Cd300a","Adrbk2","Il4i1","Mcemp1","Tmem176a","Slc38a2","Ramp1","Bcl11b","Cd3d","Cd4","Cd8b1","Cd8a","Nkg7","Cd28","Gzma","Cd3e","Nr4a1","Pou2f2","Klf4","Lyl1","tgal","Ace","Cd300a","Ceacam1","Spn","Tnfrsf1b","Il17ra","Lrp1","Adgre4","Treml4","Ear2","Stk10","Gypa","Gata1","Hbb−bs","Hbb−bt","Hba−a2","Hba−a1","Bcl3","Bhlhe40","Cbfa2t3","Cd209a","Cd300a","Cd83","Itgb7","Vrk1","Napsa","Ms4a4c"),
                            resource = c("Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment","Experiment"),
                            pmid = c("228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","sub","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882","228882"), stringsAsFactors = FALSE
)
index1="scCATCH_marker_P"
obj <- findmarkergene(obj,species = 'Mouse', tissue ='sp',marker = custom_marker ,if_use_custom_marker = TRUE, cell_min_pct = 0.25,logfc = 0.1, pvalue = 0.05)
clu_ann <- findcelltype(obj)
scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
sc_seurat_obj = AddMetaData(RDS, scCATCH_ann_data)

cellType=data.frame(ClusterID=sc_seurat_obj@meta.data$peaks_snn_res.0.4,celltype=sc_seurat_obj@meta.data$cell_type)


print(unique(cellType))
saveRDS(sc_seurat_obj, file=paste0(index1,"_combined.harmony_FindClusters.rds") )


index2="scCATCH_tissue_P"
obj <- findmarkergene(obj,species = 'Mouse', tissue ='Spleen',marker = cellmatch, cell_min_pct = 0.25,logfc = 0.1, pvalue = 0.05)
clu_ann <- findcelltype(obj)
scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
sc_seurat_obj = AddMetaData(RDS, scCATCH_ann_data)
sc_seurat_obj$seurat_clusters=sc_seurat_obj$peaks_snn_res.0.4
sc_seurat_obj@meta.data$cell_type=paste0(sc_seurat_obj@meta.data$seurat_clusters,":",sc_seurat_obj@meta.data$cell_type)
  table_cell.type=as.data.frame(table(as.character(sc_seurat_obj@meta.data$cell_type)))
  colnames(table_cell.type)=c("predicated.cell.type","number")
  table_cell.type$ratio=table_cell.type$number/colSums(table_cell.type[2])
  order_table_cell.type=table_cell.type[order(table_cell.type$number,decreasing= T),]
  write.table(order_table_cell.type,paste0(index2,"_Table2.celltypecount.csv"),sep= ",",quote = FALSE,row.names = FALSE)


saveRDS(sc_seurat_obj, file=paste0(index2,"_combined.harmony_FindClusters.rds") )


