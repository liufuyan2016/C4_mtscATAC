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
my ($out,$MT,$topN,$input,$mix,$ID,$cell);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"input:s"=>\$input,
			"top:f"=>\$topN,
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
        function: change gff format to glean format
    Usage:
        -input      input dir
        -o          outfile    must be given
	-top        top N cells according to nuc fragments
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
my $SP="SP";
`mkdir $out`;
my %hash;
my @aln=glob("$input/*MT.NUC.cell.stat.xls");
if(-s $aln[0]){
	print "$aln[0]\n";
}else{
	push @aln,$input;
}
open(SH,">$out/draw.sh");
for(my $i=0;$i<=$#aln;$i++){
   my ($key)=$aln[$i]=~/([^\/]+).MT.NUC.cell.stat.xls/;
   my $tissue;
   my $lib;
   open(IN,"$aln[$i]");
   while(<IN>){
    chomp;
    my @tmp=split/\s+/;
    if($tmp[1]=~/SP_Y_P/){
	   $tissue="Young Spleen";$lib="mtscATAC";
    }elsif($tmp[1]=~/SP_P\d/){
    	    $tissue="Aged Spleen";$lib="mtscATAC";
    }elsif($tmp[1]=~/Aged/){
	   $tissue="Aged Spleen";$lib="scATAC";
    }elsif($tmp[1]=~/Young/){
           $tissue="Young Spleen";$lib="scATAC";
    }elsif($tmp[1]=~/^SP_GS/){
	    $tissue="Aged Bone Marrow";$lib="mtscATAC";
    }else{
	next;
    }
    push @tmp,($tissue,$lib);
    push @{$hash{$tmp[1]}},[@tmp];
  }
  close IN;
  open(OUT,">$out/$key.MT.NUC.cell.stat.xls.$topN");
  foreach my $key(keys %hash){
      my @array=@{$hash{$key}};
      @array=sort{$b->[5]<=>$a->[5]} @array;
      for(my $i=0;$i<=$#array;$i++){
	   if($i<$topN){
		print OUT join("\t",@{$array[$i]})."\n"
	   }
      }
   }
   close OUT;
  %hash=();
print SH  "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls.$topN --group.col 7   --outfile $out/$key.MT.NUC.cell.stat.MT.png --value.col 4 --x.col 2 --x.lab \"Sample \" --y.lab \" MT ratio \" --title.lab \"  \"  --no.grid\n";
   print SH "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls.$topN  --group.col 7 --outfile $out/$key.MT.NUC.cell.stat.MT.fragment.png --value.col 5 --x.col 2 --x.lab \"Sample\" --y.lab \"log10 mtDNA complexity\" --title.lab \"  \"  --no.grid\n";
  print SH "Rscript $Bin/box_violin.r  --infile $out/$key.MT.NUC.cell.stat.xls.$topN  --group.col 7 --outfile $out/$key.MT.NUC.cell.stat.NUC.fragment.png --value.col 6 --x.col 2 --x.lab \"sample\" --y.lab \"log10 chromatin complexity\" --title.lab \"  \"  --no.grid\n";
}
`sh  $out/draw.sh`;
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
