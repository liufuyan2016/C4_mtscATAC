# load library
library('getopt');
require(ggplot2)
library(plyr)
library(tidyverse)


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'group.col' , 'g', 1, "integer",
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'group.lab' , 'G', 1, "character",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'legend.xpos' , 'a', 1, "double",
	'legend.ypos' , 'b', 1, "double",
	'legend.col' , 'c', 1, "integer",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
	'legend.size' , 'd', 1, "integer",
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
--group.col 	-g 	integer 	the col for group factor [forced]
--x.col 	-x 	integer 	the col for x value [forced]
--y.col 	-y 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--group.lab 	-G 	character 	the lab for group factor [forced]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--legend.xpos 	-a 	double 		the x relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.ypos 	-b 	double 		the y relative position for legend, (0.0,1.0) [optional, default: NULL]
--legend.col 	-c 	integer 	the col number for legend disply [optional, default: NULL]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 14]
--legend.size 	-d 	integer 	the font size of text for legend [optional, default: 12]
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
if ( is.null(opt$group.col) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$y.col) )	{ print_usage(spec) }
if ( is.null(opt$group.lab) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$legend.xpos ) )	{ opt$legend.xpos = NULL }
if ( is.null(opt$legend.ypos ) )	{ opt$legend.ypos = NULL }
if ( is.null(opt$legend.col ) )		{ opt$legend.col = NULL }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 14 }
if ( is.null(opt$legend.size ) )	{ opt$legend.size = 12 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }



#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip,sep="\t")
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col, opt$group.col) ){
	cat("Final Error: max(x.col, y.col, group.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df1 <- data.frame(x=as.factor(data[,opt$x.col]), group=as.factor(data[,opt$group.col]), y=data[,opt$y.col])


library("RColorBrewer")

#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
col1=brewer.pal(n = 12, name = "Paired")
col2=brewer.pal(n = 8, name = "Dark2")
col=c(col1,col2[4],col2[7],col2[8],col1,col2[4],col2[7],col2[8],col1,col2[4],col2[7],col2[8])
df = ddply(df1,'x',transform,percent_y=y/sum(y)*100)
print (df)
p <- ggplot(df,aes(x=x, y=percent_y,fill=group))+geom_bar(stat="identity")+scale_fill_manual(values=col)+coord_flip() 

##, position="dodge", fill=group)

#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- p + labs(fill=opt$group.lab) + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
	axis.text = element_text(face="bold", color="black",size=opt$axis.size))
p <- p + theme(legend.title = element_text(face="bold", size=opt$legend.size),
	legend.text = element_text(size=opt$legend.size) )
# legend position
if( (!is.null(opt$legend.xpos)) && (!is.null(opt$legend.ypos)) ){
	p <- p + theme(legend.position=c(opt$legend.xpos, opt$legend.ypos))
}
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







