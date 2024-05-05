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
        function: mt distributions
    Usage:
        -input      input mtDNA stat
        -cell        cell annotation infor(must include header:cellID\tcellType) 
        -o          outfile    must be given
	-top        top N accoring to nuc fragment
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
my %hash2;
#open(SH,">$out/draw.sh");
open(OUT,">$out/draw.mt_ratio.txt");
open(OUT2,">$out/draw.mt_ratio_each_cell_type.txt");
open(IN,"$cell");
my $H=<IN>;
my @header=split/\s+/,$H;
my $cluster;
for(my $i=0;$i<=$#header;$i++){
  if($header[$i] eq "cellType"){
	$cluster=$i+1;
	last;
  }
}
while(<IN>){
   chomp;
   my @tmp=split/\t+/;
   $hash{$tmp[0]}=$tmp[$cluster];
}
close IN;
my %hashVS;
my %hash_type;
open(IN2,"$input");
while(<IN2>){
   chomp;
   my @tmp=split/\s+/;
   next if(!defined $hash{$tmp[2]});
   my $type=$hash{$tmp[2]};
   my $sample=$tmp[1];
   if($tmp[1]=~/_Y_P/){
	$tmp[1]="Young spleen";
   }elsif($tmp[1]=~/_P/){
	 $tmp[1]="Aged spleen";
   }elsif($tmp[1]=~/SP_GS/){
	  $tmp[1]="Aged bone marrow";
   }else{
	next;
   }
   $hashVS{$sample}=$tmp[1];
   push @{$hash_type{$sample}{$type}},$tmp[3];
   print OUT join("\t",@tmp)."\t$type\n";
}
close IN2;
foreach my $key(keys %hash_type){
       foreach my $key2(keys %{$hash_type{$key}}){
 		my @arr=@{$hash_type{$key}{$key2}};
                my $median=&median(@arr);
		print OUT2 "$key\t$key2\t$median\t$hashVS{$key}\n";
	}
}             
close OUT2;
print "/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Software/miniconda3/bin/Rscript  $Bin/twoFactorBox.r --infile  $out/draw.mt_ratio.txt --group.col 2 --outfile $out/statMT.png --value.col 4 --x.col 7  --x.lab \"Cell type\" --y.lab \"mtDNA fragments (%)\" --title.lab \"  \"  --no.grid  --group.lab  \"Sample\" --width 4000 --height 1600\n";
`/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Software/miniconda3/bin/Rscript  $Bin/twoFactorBox.r --infile  $out/draw.mt_ratio.txt --group.col 2 --outfile $out/statMT.png --value.col 4 --x.col 7  --x.lab "Cell type" --y.lab "mtDNA fragments (%)" --title.lab "  "  --no.grid  --group.lab  "Sample" --width 4000 --height 1600`;
print "/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Software/miniconda3/bin/Rscript  $Bin/twoFactorBox.r --infile  $out/draw.mt_ratio_each_cell_type.txt --group.col 4 --outfile $out/mt_ratio_each_cell_type.png --value.col 3 --x.col 2  --x.lab \"Cell type\" --y.lab \"mtDNA fragments (%)\" --title.lab \"  \"  --no.grid  --group.lab  \"Sample\" --width 4000 --height 1600\n";
`/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Software/miniconda3/bin/Rscript  $Bin/twoFactorBoxV2.r --infile  $out/draw.mt_ratio_each_cell_type.txt --group.col 4 --outfile $out/mt_ratio_each_cell_type.png --value.col 3 --x.col 2  --x.lab "Cell type" --y.lab "mtDNA fragments (%)" --title.lab "  "  --no.grid  --group.lab  "Sample" --width 4000 --height 1600`;


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub median{
   my (@numbers)=@_;
   my @sorted_numbers = sort {$a <=> $b} @numbers;

   my $count = @sorted_numbers;
   my $median;
   if ($count % 2 == 0) {
         $median = ($sorted_numbers[$count/2 - 1] + $sorted_numbers[$count/2]) / 2;
    } else {
         $median = $sorted_numbers[int($count/2)];
    }
    return $median;
}


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
