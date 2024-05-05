#!/share/nas1/liux/Software/R.3.6.0/bin/Rscript

# usage function
usage <- function() {
        print("-------------------------------------------------------------------------------")
        print("Usage: Rscript ComplexHeatmap.r input_files outdir group order_names height width")
        print("1) input_files: files path split by \",\", such as [file1,file2,file3]")
        print("2) outdir: the dir for output")
	print("3) group: group files path split by \",\", such as [file1.group,file2.group,file3.group]")
	print("4) order_names: order names for each input_file, such as [G1,G2,G3]")
	print("5) height: height for pic")
	print("6) width: width for pic")
	print("7) ratio to make white color")
	print("8) make log")
        print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) < 7 ) {
        print(args)
        usage()
        stop("the length of args < 7")
}

height = as.numeric(args[5])
width = as.numeric(args[6])
ratio = as.numeric(args[7])
str(height)
str(width)

# load library
#require(edgeR)
#require(ggplot2)
require(circlize)
library("ComplexHeatmap")
library("Cairo")
library(RColorBrewer)
library(stringr)

order_names = strsplit(args[4],split = ",")[[1]]
#print(order_names)

# anno color
#anno_col = brewer.pal(9, "Set1")
anno_col = brewer.pal(12, "Paired")
files = strsplit(args[1],split = ",")[[1]]
groups = strsplit(args[3],split = ",")[[1]]

# color set
#my_col =(colorRampPalette(colors=c(rep("white",4),c("red",2))))(400)

ht_list1 <- 0
num <- 0
for (i in files){
  mat = read.delim(i, header=T, check.names=F, row.names=1)
  logN=as.numeric(args[8]) 
  if(logN==1){
   	mat_log = log10(mat+1)
   }else{
	mat_log = mat
   }
  # mat_log =mat
  cd<-c(as.matrix(mat_log))
  ht_list1 <-c(ht_list1,cd)
  num <- num+1
  if(num==1){
	mat1 <- sort(cd)
	combine_matrix <- mat_log
	c1<-ncol(combine_matrix)
	num_index1=c(1)
	num_index2=c(c1)
  }else{
	
	c1<-num_index2[num-1]+1
	c2<-ncol(mat_log)+num_index2[num-1]
	num_index2=c(num_index2,c2)
	num_index1=c(num_index1,c1)
        combine_matrix <-cbind(combine_matrix,mat_log)
  }
}
mat_scaled_all = t(apply(combine_matrix, 1, scale))






#color1=(colorRampPalette(colors="white"))(30)
#col2 =(colorRampPalette(colors=c(rep("white",4),"red")))(50)
#color2=(colorRampPalette(colors="red"))(30)
#my_col=c(color1,col2,color2)


ct=0
data_list = list()
ht = list()
for (i in files){
	ct = ct + 1
	mat = read.delim(i, header=T, check.names=F, row.names=1)
	group = read.delim(groups[ct],header=F,check.names=F)
	groups_types = group[,2][!duplicated(group[,2])]
	my_anno_col = list(GroupN = anno_col[1:length(groups_types)])
	print(my_anno_col)
	names(my_anno_col$GroupN) = groups_types
	my_anno_col2 = list(GroupT = anno_col[1:length(groups_types)])
        names(my_anno_col2$GroupT) = groups_types

	my_anno_col3 = list(Group = anno_col[1:length(groups_types)])
        names(my_anno_col3$Group) = groups_types
	print(1)

	group[,1] = as.character(group[,1])
	group[,2] = as.character(group[,2])
	minvalue<-min(na.omit(mat_scaled_all))
	mat_scaled_all[is.na(mat_scaled_all)]<- minvalue
	
	mat_scaled=mat_scaled_all[,num_index1[ct]:num_index2[ct]]
	colnames(mat_scaled) = colnames(mat)
	mat_scaled = as.matrix(mat_scaled[,group[,1]])
	data_list[[ct]] = mat_scaled
	#####for legend
	seqdata <- sort(mat_scaled_all)
	if(ct == 1){
		mat1 <- sort( mat_scaled)
		if(ratio==-1){
			white_col=ceiling(0.5 *length(mat_scaled_all))
		}else{
			max_value=mat1[ceiling(ratio*length(mat1))]
			white_col_array <- which(seqdata==max_value,arr.ind = TRUE)	
			white_col=white_col_array[1]
		}
       		#N <- ceiling(ratio *length(mat_scaled_all))
		white_colRatio <- white_col/length(mat_scaled_all)
		if(white_colRatio >0.85){
 		       white_col <- ceiling(0.8*length(mat_scaled_all))
		}
		
		piancha <- ceiling(0.15*length(mat_scaled_all))
		redcol<-length(mat_scaled_all)
		print (white_col)
		print (piancha)
		print (redcol)
		my_col=(colorRampPalette(colors=c("blue","white","red"))(redcol))
	}
	if (ct == 1){
		print (length(seqdata))
        print (length(my_col))
	
		groupN<-unique(group[,2])
		if(length(groupN)==1){
                       ha <- NULL
                }else{
                        ha = HeatmapAnnotation(GroupN = group[,2], col = my_anno_col,show_annotation_name =F)
                }
		print (my_anno_col)
			
		print(group[,2])
		ht[[ct]] = Heatmap(data_list[[ct]],  col=colorRamp2(seqdata, my_col ), name = order_names[ct],
		cluster_columns = F,   cluster_rows =F, show_row_dend = FALSE, show_row_names = T,show_column_names = F, column_title =  order_names[ct],show_column_dend=FALSE, top_annotation=ha,
		row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = "grey",lwd=0.5),heatmap_legend_param=list( border="WhiteSmoke"))
		ht_list = ht[[ct]]
	}else if(ct == length(groups)){
		groupN<-unique(group[,2])
                if(length(groupN)==1){
                       ha <- NULL
                }else{
			ha = HeatmapAnnotation(GroupT = group[,2], col = my_anno_col2,show_annotation_name =F)
		}
		
		ht[[ct]] = Heatmap(data_list[[ct]], col=colorRamp2(seqdata, my_col ), name = order_names[ct], show_row_names = T,
		show_column_names =F,cluster_columns = F,cluster_rows =F, row_title = order_names[ct],column_title =  order_names[ct],top_annotation=ha,
		row_names_gp = gpar(fontsize = 5),rect_gp = gpar(col = "grey",lwd=0.5),heatmap_legend_param=list( border="WhiteSmoke"))
		ht_list = ht_list + ht[[ct]]
	}else{
		 groupN<-unique(group[,2])
                if(length(groupN)==1){
                        ha <- NULL
                }else{
			ha = HeatmapAnnotation(Group = group[,2], col = my_anno_col3,show_annotation_name =F)
		}
		ht[[ct]] = Heatmap(data_list[[ct]],  col=colorRamp2(seqdata, my_col ), name = order_names[ct], show_row_names = F,
		show_column_names = F,cluster_rows =F, cluster_columns=F,column_title =  order_names[ct],top_annotation = ha,
		row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = "grey",lwd=0.5),heatmap_legend_param=list( border="WhiteSmoke"))
		ht_list = ht_list + ht[[ct]]
        }
}


print("Draw heatmap:")

CairoPNG(paste0(args[2],".ComplexHeatmap.png",sep = ""), height = height, width = width, res = 500)
draw(ht_list)
dev.off()
pdf(file=paste0(args[2],".ComplexHeatmap.pdf",sep = ""), height = height/500, width = width/500)
draw(ht_list)
dev.off()


#CairoTIFF(paste0(args[2],"/ComplexHeatmap.tiff",sep = ""), height = height, width = width, res = 300)
#draw(ht_list)
#dev.off()







