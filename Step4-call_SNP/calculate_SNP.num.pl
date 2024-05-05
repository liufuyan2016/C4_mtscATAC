#!/usr/bin/perl -w
# Copyright (c) BGI 2015/6/28
# Writer:         liufuyan <liufuyan@genomics.cn>
# Program Date:   2015/6/28.
# Modifier:       liufuyan <liufuyan@genomics.cn>
# Last Modified:  2015/6/28.

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);


my $programe_dir=basename($0);
my $path=dirname($0);

my $ver    = "1.0";
my $Writer = "liufuyan <liufuyan\@echo>";
my $Data   = "2015/6/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my (@in,$final,$file,$ann,$snp,$top,$pvalue,$od,@x,$cell,$index);
GetOptions(
			"h|?" =>\&help,
			"od:s"=>\$od,
			"i:s{,}"=>\@in,
			"cell:s"=>\$cell,
			"snp:s"=>\$snp,
			"final:s"=>\$final,
			"ann:s"=>\$ann,
			"top:s"=>\$top,
			 "index:s"=>\$index,	
			) || &help;
&help unless ($index && $od);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: 
    Usage:
        -i          infile     must be given
        -final       infile  final SNP results
	-index      index
	-cell       cell number 
	-ann        annotation SNP
	-top        top list to draw
	-snp        snp all
        -od         outdir    must be given
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
my %hash2;
my %hash3;
my @all1;
my @all2;
my %hash4;
`mkdir $od`;
my %hashselect;
if(defined $top){
	my $N=0;
	open (IN1, "<$top")|| die "cannot open $cell:$!";
	while (<IN1>){
 		 chomp;
  		 my @tmp=split/\t+/;
  		 $hashselect{$tmp[0]}=$N++;
      }
      close IN1;
}
my %hashfinal;
if(defined $final){
    open (IN1, "<$final")|| die "cannot open $cell:$!";
     while (<IN1>){
                 chomp;
                 my @tmp=split/\t+/;
                 $hashfinal{$tmp[0]}=1;
      }
      close IN1;
}
open(OUT8,">$od/$index.snp.all.stat");
open(OUT1,">$od/$index.snp.cell.txt");
open(OUT11,">$od/$index.snp.cell.txt.vaf");
open(OUT2,">$od/$index.snp.num.txt");
open(OUT3,">$od/$index.each_cell_snp.num.txt");
open(OUT88,">$od/$index.each_cell_snp.num.draw");

my %hash_ann_detail;
my %hashT;
my %hash_ann;
my %hashall;
if(defined $ann){
   open (IN1, "<$ann")|| die "cannot open $ann:$!";
   while (<IN1>){
       chomp;
       my @tmp=split/\t+/;
        my @infor=split/\:/,$tmp[3];

       $hash_ann{$tmp[0]}="$infor[0]\t$tmp[1]";
	my $SNP="$tmp[6]$tmp[7]>$tmp[8]";
       $hash_ann_detail{$tmp[0]}="$SNP\t$infor[0]\t$index";
   }
   close IN1;
}
open (IN1, "<$cell")|| die "cannot open $cell:$!";
while (<IN1>){
  chomp;
  my @tmp=split/\t+/;
  $hash{$tmp[0]}=$tmp[-1];
}
close IN1;
my %hash33;
my %hashall_cell;
my %hashall_cellselect;
my $Total;
open (IN2, "<$in[0]")|| die "cannot open $in[0]:$!";
my $header=<IN2>;
chomp $header;
my  $start;
my %hash_ann_stat;
my @cell=split/\s+/,$header;
while (<IN2>) {
	chomp;
	next if(/^$/);
	my @tmp=split/\s+/;
        if(defined $final){
		next if(!defined $hashfinal{$tmp[0]});
	}
        $Total=$#tmp;
	for(my $i=1;$i<=$#tmp;$i++){
		if($cell[0] eq "#ID"){
			$start=$i;
		}else{
			 $start=$i-1;
		}
		if($tmp[$i]==0){
			$hash3{$tmp[0]}{$hash{$cell[$start]}}{$cell[$start]}=0;
			$hash33{$tmp[0]}{$hash{$cell[$start]}}{$cell[$start]}=0;

			next;
		}else{
			$hashall{$tmp[0]}++;
			$hashall_cell{$cell[$start]}++;
			if(defined $hashselect{$tmp[0]}){
				 $hashall_cellselect{$cell[$start]}++;
			}
		}
		$hash2{$tmp[0]}{$hash{$cell[$start]}}++;
		if(defined $ann){
			$hash_ann_stat{$hash_ann{$tmp[0]}}{$hash{$cell[$start]}}++;
		}
		$hash3{$tmp[0]}{$hash{$cell[$start]}}{$cell[$start]}="1|$tmp[$i]";		
                $hash33{$tmp[0]}{$hash{$cell[$start]}}{$cell[$start]}=$tmp[$i];

		$hash4{$cell[$start]}++;
	}
  }
close IN2;
for(my $i=0;$i<=$#cell;$i++){
	next if(!defined $hash{$cell[$i]});
	$hashT{$hash{$cell[$i]}}++;
}
my %hash5;
my %hashVAF;
foreach my $key(sort keys %hash2){
	foreach my $key2(sort keys %{$hash2{$key}}){
		my $ratio=$hash2{$key}{$key2}/$hashT{$key2}*100;
		my ($snp)=$key=~/\d+(.*)/;		

		print OUT2 "$key2\t$key\t$hash2{$key}{$key2}\t$ratio\t$index\n";


                push @{$hash5{$key}},"$key2\t$key\t$hash2{$key}{$key2}\t$ratio\t$index";
		foreach my $key3(sort keys %{$hash3{$key}{$key2}}){
			my $type="\n";
			if(defined $ann){
				$type="\t$hash_ann_detail{$key}\n";
			}
			my($sample)=$key3=~/(.*)\_BC/;
			print OUT1 "$key3\t$key2\t$key\t$hash3{$key}{$key2}{$key3}\t$index\t$snp\t$sample$type";
			#if($hash33{$key}{$key2}{$key3}!=0){
			push @{$hashVAF{$sample}{$key2}},$hash33{$key}{$key2}{$key3};	
			#}
			print OUT11 "$key3\t$key2\t$key\t$hash33{$key}{$key2}{$key3}\t$index\t$snp\t$sample$type";
		}			
	}
}
close OUT2;
close OUT1;
close OUT11;
########median
open(OUT66,">$od/$index.snp.cell.txt.vaf.draw");
foreach my $key(keys %hashVAF){
   foreach my $key2(keys %{$hashVAF{$key}}){
	my @arr=@{$hashVAF{$key}{$key2}};
        my $average=&average(@arr);
	print OUT66 "$key\t$key2\t$average\t$index\n";
   }
}
close OUT66;

my %hashSNP_count;
foreach my $key(keys %hash4){
	 my ($sample)=$key=~/(.*)\_BC/;
	 #if($hash4{$key}!=0){
	 push @{$hashSNP_count{$sample}{$hash{$key}}},$hash4{$key};
	 #}
	 print OUT3 "$key\t$hash{$key}\t$hash4{$key}\t$sample\t$index\n";
}
close OUT3;

foreach my $key(keys %hashSNP_count){
    foreach my $key2(keys %{$hashSNP_count{$key}}){
	  my @arr=@{$hashSNP_count{$key}{$key2}};
          my $average=&average(@arr);
	  print OUT88 "$key\t$key2\t$average\t$index\n";
    }
}
close OUT88;



foreach my $key(keys %hashall){
    my $ratio=$hashall{$key}/$Total*100;
    print OUT8 "$key\t$hashall{$key}\t$ratio\n";
}
my %hashall_cell_stat;
my %hashall_cell_stat2;
foreach my $key(keys %hashall_cell){
	$hashall_cell_stat{$hashall_cell{$key}}++;
}

foreach my $key(keys %hashall_cell_stat){
    my $ratio=$hashall_cell_stat{$key}/$Total*100;
    print OUT8 "$key SNP\t$hashall_cell_stat{$key}\t$ratio\n";
}
close OUT8;


if(defined $top){
	foreach my $key(keys %hashall_cellselect){
                $hashall_cell_stat2{$hashall_cellselect{$key}}++;
	}
        open(OUT41,">$od/$index.snp.num.txt.top.stat");
	foreach my $key(keys %hashall_cell_stat2){
	    my $ratio=$hashall_cell_stat2{$key}/$Total*100;
    	    print OUT41 "$key SNP\t$hashall_cell_stat{$key}\t$ratio\n";
        }
        close OUT41;

	open(OUT4,">$od/$index.snp.num.txt.top");
	foreach my $key(sort {$hashselect{$a}<=>$hashselect{$b}} keys %hashselect){
		my @arr=@{$hash5{$key}};
		foreach my $key2(@arr){
			print OUT4 "$key2\n";
		}
	}
	close OUT4;
	`Rscript  $Bin/group_Bar2.r  --infile $od/$index.snp.num.txt.top --outfile $od/$index.snp.num.top.png --group.col 2 --x.col 1 --y.col 4 --group.lab "SNPs" --x.lab " " --y.lab "Percent of mutated cells(%)" --no.grid  --height 1400 --width 3000 --axis.size 9  --lab.size 9 --legend.size 7`;
}else{
`Rscript  $Bin/group_Bar2.r  --infile $od/$index.snp.num.txt --outfile $od/$index.snp.num.png --group.col 2 --x.col 1 --y.col 4 --group.lab "SNPs" --x.lab " " --y.lab "Percent of mutated cells(%)" --no.grid  --height 1400 --width 3000 --axis.size 9  --lab.size 9 --legend.size 7`;

}	

print "Rscript $Bin/oneFactorViolin.r --infile $od/$index.each_cell_snp.num.txt --outfile $od/$index.each_cell_snp.num.png --value.col 3 --x.col 2 --x.lab \"Cell Type\" --y.lab \"Number of variants\" --no.grid   --height 1400 --width 3000\n";


`Rscript $Bin/oneFactorViolin.r --infile $od/$index.each_cell_snp.num.txt --outfile $od/$index.each_cell_snp.num.png --value.col 3 --x.col 2 --x.lab "Cell Type" --y.lab "Number of variants" --no.grid   --height 1400 --width 3000`;
print "Rscript  $Bin/group_Bar2.r  --infile $od/$index.snp.num.txt --outfile $od/$index.snp.num.png --group.col 2 --x.col 1 --y.col 4 --group.lab \"SNPs\" --x.lab \" \" --y.lab \"Percent of mutated cells(%)\" --no.grid  --height 1400 --width 3000 --axis.size 9  --lab.size 9 --legend.size 7\n";






###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub average{
   my (@numbers)=@_;
  my $sum = 0;
  foreach my $num (@numbers) {
    $sum += $num;
  }
   my $count = scalar @numbers;
   my $average = $sum / $count;
   return $average;
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
