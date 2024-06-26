# load library
library('getopt');
require(ggplot2)
require(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggsignif)

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'value.col' , 'x', 1, "integer",
	'x.col' , 'y', 1, "integer",
	'group.col' , 'g', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'group.lab' , 'G', 1, "character",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
	'legend.xpos' , 'A', 1, "double",
	'legend.ypos' , 'B', 1, "double",
	'legend.col' , 'C', 1, "integer",
	'legend.size' , 'D', 1, "integer",
	'no.grid' , 'r', 0, "logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--value.col 	-x 	integer 	the col for x value [forced]
--x.col 	-y 	integer 	the col for x value [forced]
--group.col 	-g 	integer 	the col for group factor [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--group.lab 	-G 	character 	the lab for group factor [forced]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 14]
--legend.xpos 	-A 	double 		the x relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.ypos 	-B 	double 		the y relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.col 	-C 	integer 	the col number for legend disply [optional, default: NULL]
--legend.size 	-D 	integer 	the font size of text for legend [optional, default: 12]
--no.grid	-r 	NULL 		Do not drawing grid
--skip 		-k 	integer 	the number of line for skipping [optional, default: 0]
\n")
	q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$value.col) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$group.col) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }
if ( is.null(opt$group.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }
if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = NULL }
if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = NULL }
if ( is.null(opt$legend.col ) )		{ opt$legend.col = NULL }
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 12 }




#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, sep="\t", skip=opt$skip)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$value.col, opt$x.col, opt$group.col) ){
	cat("Final Error: max(value.col, opt$x.col, group.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df <- data.frame(y=data[,opt$value.col], x=as.factor(data[,opt$x.col]), group=as.factor(data[,opt$group.col]))



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
max1=max(df$y)
max2=max1+10
p=ggplot(data=df,aes(x=x,y=y,fill=group))+geom_boxplot(width=0.8,
               position = position_dodge(0.9))+geom_dotplot(binaxis = "y", stackdir = "center",
                position = position_dodge(1))+
  scale_fill_manual(values = c("#FBC29C","#A5DDC2"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method = "wilcox.test",
                     hide.ns = TRUE,
			label.y = c(90,95)
                   # label.y = c(max1,max2)
	)
  #geom_signif(annotations = c("",""),
   #           y_position = c(25,30),
   #           xmin = c(0.8,1.8),
   #           xmax = c(1.2,2.2),
   #           tip_length = c(0.02,0.4,0.02,0.2))



#p <- p + geom_density(aes(fill=group, colour=group), alpha = 0.2)

#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- p + labs(fill=opt$group.lab) + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
	axis.text.x = element_text(face="bold", size=opt$axis.size,angle=45,hjust=1,vjust=1),axis.text.y = element_text(face="bold", size=opt$axis.size))
# remove legend
#p <- p + theme(legend.position = "none")
# legend col
if( is.null(opt$legend.col) ){
	levels_num <- length( levels(df$group) )
	if( levels_num%%12==0 )
		opt$legend.col <- as.integer(levels_num / 12)
	else
		opt$legend.col <- as.integer(levels_num / 12) + 1
}
p <- p + guides(fill = guide_legend(ncol = opt$legend.col))
# grid and background
if ( !is.null(opt$no.grid) ) {
	p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
}



#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=opt$outfile, height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







