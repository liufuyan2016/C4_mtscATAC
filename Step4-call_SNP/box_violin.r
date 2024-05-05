#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');
require(ggplot2)
require(RColorBrewer)


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
        'group.col' , 'G', 1, "integer",
	'x.col' , 'g', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'lab.size' , 'l', 1, "integer",
	'axis.size' , 's', 1, "integer",
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
1) Rscript oneFactorBox.r --infile in_oneFactorBox.data --outfile out_oneFactorBox.png \\
	--value.col 2 --x.col 1 --x.lab \"x lab\" --y.lab \"y lab\" --title.lab \"title lab\" --skip 1

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--value.col 	-x 	integer 	the col for x value [forced]
--group.col     -g      integer         the col for group factor [forced]
--x.col 	-g 	integer 	the col for group factor [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 14]
--axis.size 	-s 	integer 	the font size of text for axis [optional, default: 11]
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
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 1500 }
if ( is.null(opt$width ) )		{ opt$width = 2000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 11 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 11 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }




#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile,sep = "\t", skip=opt$skip)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$value.col, opt$x.col) ){
	cat("Final Error: max(value.col, x.col) > the col of infile\n")
	print_usage(spec)
}
# create df
df <- data.frame(y=data[,opt$value.col], x=data[,opt$x.col],Replication=data[,opt$group.col])
df$x=factor(df$x,level=unique(df$x))


xxx=mean(df$y)
print (xxx)
#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
p <- ggplot(df, aes(x = x, y =y,fill=Replication)) +
geom_violin(position = position_dodge(0.9),alpha = 0.8,
              trim = TRUE,
              color = "black") +
  # 箱线图图层
  geom_boxplot(width = .1,outlier.shape = NA,show.legend = F,
               position = position_dodge(0.9),
               alpha = 0.5) +
  # 主题调整
  theme_bw(base_size =opt$axis.size) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),axis.text.y = element_text(size=opt$axis.size),
        legend.position = 'top',aspect.ratio = 0.4) +scale_fill_manual(values = c("#2A7080","#E64C35","#4B93C3"))+
 theme(legend.position = "bottom") + xlab(opt$x.lab) + ylab(opt$y.lab)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

# set lab and axis test size
# remove legend
p <- p + theme(legend.position = "none")
# grid and background
if ( !is.null(opt$no.grid) ) {
	p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
		panel.grid = element_blank())
}
theme_classic2<- function (base_size =opt$axis.size, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black",angle=45,hjust=1,vjust=1),
          axis.text.y = element_text(colour = "black" ),
          axis.title.x = element_text(colour = "black" ),
          axis.title.y = element_text( colour = "black",angle=90),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5))
}

p <-p+theme_classic2()


#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=opt$outfile, height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







