#!/usr/bin/perl -w
# Copyright (c) BGI 2022/6/28
# Writer:         liufuyan <liufuyan@genomics.cn>
# Program Date:   2022/6/28.
# Modifier:       liufuyan <liufuyan@genomics.cn>
# Last Modified:  2022/6/28.

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);


my $programe_dir=basename($0);
my $path=dirname($0);

my $ver    = "1.0";
my $Writer = "liufuyan <liufuyan\@genomics.cn>";
my $Data   = "2022/6/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($out,$MT,$input,$mix,$ID,$cell);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"input:s"=>\$input,
			"ID:s"=>\$ID,
			"mix:s"=>\$mix,
			"cell:s"=>\$cell,
			"MT:s"=>\$MT,		
			) || &help;
&help unless ($input && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: make stat for each cell MT fragments and MT ratio 
    Usage:
        -input      input all samples of C4_mtscATAC out dir
	-MT         MT id, eg. chrM 
        -o          outfile    must be given
        -h          Help document
	Usage End.
	exit;
}
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################
my %hash;
my $SP="SP";
`mkdir $out`;
my @aln=glob("$input/*/*/aln.bed");
for(my $i=0;$i<=$#aln;$i++){
	my ($sample)=$aln[$i]=~/$input\/(.*)\/.*\//;
        #next if($sample =~/JiaQuan/);
	my $d2c;
        if($aln[$i]=~/02.chromap/){
                $d2c="03.d2cfile";
        }else{
                $d2c="d2cfile";
        }

	 print "perl $Bin/classify_MT_to_cells.pl -aln  $aln[$i] -cell $input/$sample/$d2c/$sample.barcodeMerge.tsv  -MT $MT  -o  $out/$sample.MT.fragments.tsv\n";
	`perl $Bin/classify_MT_to_cells.pl -aln  $aln[$i] -cell $input/$sample/$d2c/$sample.barcodeMerge.tsv  -MT $MT  -o  $out/$sample.MT.fragments.tsv`;
			
	if(-s "$input/$sample/$d2c/$sample.fragments.tsv.gz"){	
		`rm $input/$sample/$d2c/fragments.tsv`;
		`cp $input/$sample/$d2c/$sample.fragments.tsv.gz $input/$sample/$d2c/fragments.tsv.gz`;   
    		`gzip -d $input/$sample/$d2c/fragments.tsv.gz`;
	}
       my $aln1="$input/$sample/$d2c/fragments.tsv";
       open (IN, "<$aln1")|| die "cannot open $aln1:$!";
       while (<IN>) {
	  chomp;
	  next if(/^$/);
	  my @tmp=split/\s+/;
	  if(defined $mix){
		($SP)=$tmp[0]=~/^(.*)\_/;
	  }
	  $hash{$SP}{$sample}{$tmp[3]}{"NUC"}+=$tmp[-1];
     }
     close IN;
    open (IN2, "<$out/$sample.MT.fragments.tsv")|| die "cannot open :$!";
    while (<IN2>){
        chomp;
        next if(/^$/);
        my @tmp=split/\s+/;
	if(defined $mix){
		($SP)=$tmp[0]=~/^(.*)\_/;
        }
	$hash{$SP}{$sample}{$tmp[3]}{"MT"}+=$tmp[-1];
    }        
}

open(SH,">$out/draw.sh");
foreach my $key(keys %hash){
   open(OUT,">$out/$key.MT.NUC.cell.stat.xls");
   foreach my $key2(keys %{$hash{$key}}){
         foreach my $key3(keys %{$hash{$key}{$key2}}){
              if(!defined $hash{$key}{$key2}{$key3}{"NUC"}){
	   	    $hash{$key}{$key2}{$key3}{"NUC"}=0;
   	       }
	      if(!defined $hash{$key}{$key2}{$key3}{"MT"}){
	   	    $hash{$key}{$key2}{$key3}{"MT"}=0;
	       }
	      my $ratio=$hash{$key}{$key2}{$key3}{"MT"}/($hash{$key}{$key2}{$key3}{"MT"}+$hash{$key}{$key2}{$key3}{"NUC"})*100;
	       my $MT_F;
		my  $NUC_F;
               if($hash{$key}{$key2}{$key3}{"MT"}==0){
			$MT_F=-1;
		}else{
			$MT_F=log($hash{$key}{$key2}{$key3}{"MT"}) / log(10);
		}
		if($hash{$key}{$key2}{$key3}{"NUC"} ==0){
			$NUC_F=-1;	
		}else{
			$NUC_F=log($hash{$key}{$key2}{$key3}{"NUC"}) / log(10);
		}
                print OUT "$key\t$key\_$key2\t$key3\t$ratio\t$MT_F\t$NUC_F\n"; 
            }
   }
   close OUT;
   print SH  "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls  --outfile $out/$key.MT.NUC.cell.stat.MT.png --value.col 4 --x.col 2 --x.lab \"Sample \" --y.lab \" MT ratio \" --title.lab \"  \"  --no.grid\n";
   print SH "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls  --outfile $out/$key.MT.NUC.cell.stat.MT.fragment.png --value.col 5 --x.col 2 --x.lab \"Sample\" --y.lab \"log10 mtDNA complexity\" --title.lab \"  \"  --no.grid\n";
  print SH "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls  --outfile $out/$key.MT.NUC.cell.stat.NUC.fragment.png --value.col 6 --x.col 2 --x.lab \"sample\" --y.lab \"log10 chromatin complexity\" --title.lab \"  \"  --no.grid\n";
}
`sh $out/draw.sh`;
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time : [",&sub_time($t),"]\n";
}
sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
	$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'f_dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/";
                chdir($pwd);
        }
		 elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}
