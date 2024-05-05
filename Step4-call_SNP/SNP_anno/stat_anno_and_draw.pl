#!/usr/bin/perl -w
# Copyright (c) ECO 2015/6/28
# Writer:         liufuyan <liufuyan@echo>
# Program Date:   2015/6/28.
# Modifier:       liufuyan <liufuyan@echo>
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
my (@in,$file,$ann,$snp,$top,$pvalue,$od,@x,$cell,$index);
GetOptions(
			"h|?" =>\&help,
			"od:s"=>\$od,
			"i:s{,}"=>\@in,
			"cell:s"=>\$cell,
			"snp:s"=>\$snp,
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
	-index      index
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
my %hashSNP;
my %hashT;
my %hash_mut;
my %hashVS;
my %hashfun;
my %hashG;
my %hashcell;
`mkdir $od`;
open(SH,">$od/$index.draw.sh");
open(OUT1,">$od/$index.snp.snp.xls");
open(OUT2,">$od/$index.snp.gene.xls");
open(OUT4,">$od/$index.snp.type.xls");
open(OUT3,">$od/$index.SNP_mutation_num.xls");
open(OUT5,">$od/$index.SNP_mutation_num.for_heatmap.txt");

open (IN2, "<$in[0]")|| die "cannot open $in[0]:$!";
while (<IN2>) {
	chomp;
	next if(/^$/);
	my @tmp=split/\t+/;
	#$hash_mut{$tmp[1]}{$tmp[0]}{$tmp[3]}++;
	$hashVS{$tmp[0]}=$tmp[-1];
        $hashcell{$tmp[1]}{$tmp[0]}=1;
	my @x=split/\|/,$tmp[3];
	$tmp[3]=$x[0];
	if($tmp[3]==1){
		$hash_mut{$tmp[1]}{$tmp[0]}{$tmp[3]}++;
		$hashfun{$tmp[-1]}{$tmp[1]}{$tmp[7]}++;
		$hashG{$tmp[-1]}{$tmp[1]}{$tmp[6]}++;
		$hashSNP{$tmp[-1]}{$tmp[1]}{$tmp[5]}++;
		$hashT{$tmp[-1]}{$tmp[1]}++;
		
	}
}
close IN2;
foreach my $key(sort keys %hashSNP){
	foreach my $key2( sort keys %{$hashSNP{$key}}){
		foreach my $key3( sort keys %{$hashSNP{$key}{$key2}}){
			print OUT1 "$key $key2\t$key3\t$hashSNP{$key}{$key2}{$key3}\n";
		}
		foreach my $key3( sort keys %{$hashG{$key}{$key2}}){
                        print OUT2 "$key $key2\t$key3\t$hashG{$key}{$key2}{$key3}\n";
                }
		foreach my $key3( sort keys %{$hashfun{$key}{$key2}}){
                        print OUT4 "$key $key2\t$key3\t$hashfun{$key}{$key2}{$key3}\n";
                }

	}
}

my %hash_cell_Type;
my %cell_type;
my %hash_mutN2;
my %hash_mutN;
my $N=0;
foreach my $key( keys %hashcell){
  foreach my $key2( keys %{$hashcell{$key}}){
	if(defined $hash_mut{$key}{$key2}{"1"}){		
		#$hash_mut{$key}{$key2}{"1"}=$hash_mut{$key}{$key2}{"1"}>10?"10+":$hash_mut{$key}{$key2}{"1"};
		$cell_type{$key}=1;	
		$hash_mutN{$key}{$hash_mut{$key}{$key2}{"1"}}{$hashVS{$key2}}++;
		$hash_mutN2{$hashVS{$key2}}{$hash_mut{$key}{$key2}{"1"}}{$key}++;
	}else{
		$hash_mutN{$key}{"0"}{$hashVS{$key2}}++;	
		$hash_mutN2{$hashVS{$key2}}{"0"}{$key}++;
	}
	 $hash_cell_Type{$hashVS{$key2}}{$key}++;
   }
}
print "$N\n";
foreach my $key(sort keys %hash_mutN){
    foreach my $key2( sort {$b<=>$a} keys %{$hash_mutN{$key}}){
		foreach my $key3( sort keys %{$hash_mutN{$key}{$key2}}){
			print OUT3 "$key\t$key2\t$hash_mutN{$key}{$key2}{$key3}\t$key $key3\n";
		}
    }
}

print OUT5 "ID";
foreach my $key(sort  keys %hash_mutN2){
    foreach my $key2( sort keys %{$hash_mutN2{$key}}){
	 foreach my $key3(sort  keys %cell_type){
		print OUT5 "\t$key3";
	}
	print OUT5 "\n";
	last;
     }
     last;
}
foreach my $key( sort{$b cmp $a} keys %hash_mutN2){
    foreach my $key2( sort{$a<=>$b} keys %{$hash_mutN2{$key}}){
		print OUT5 "$key\_$key2";
		foreach my $key3(sort  keys %cell_type){
			if(!defined $hash_mutN2{$key}{$key2}{$key3}){
				$hash_mutN2{$key}{$key2}{$key3}=0;
			}
			if(!defined $hash_cell_Type{$key}{$key3}){
				print OUT5 "\t0";
			}else{		
				my $ratio=$hash_mutN2{$key}{$key2}{$key3}/$hash_cell_Type{$key}{$key3};
				print OUT5 "\t$ratio";
			}
		}
		print OUT5 "\n";
                
    }
}


print SH "Rscript  $Bin/dodgedBar_percent.r --infile $od/$index.snp.snp.xls  --outfile $od/$index.snp.snp.png --group.col 2 --x.col 1   --y.col 3  --x.lab \"Cell Type\" --group.lab \"SNPs\" --y.lab \"Proportion\" --title.lab \" \" --legend.col 2  --width 6000 --no.grid\n";
print SH "Rscript  $Bin/dodgedBar_percent.r --infile $od/$index.snp.gene.xls  --outfile $od/$index.snp.gene.png --group.col 2 --x.col 1   --y.col 3  --x.lab \"Cell Type\" --group.lab \"mt gene\" --y.lab \"Proportion\" --title.lab \" \" --legend.col 2  --width 6000 --no.grid\n";
print SH "Rscript  $Bin/dodgedBar_percent.r --infile $od/$index.snp.type.xls  --outfile $od/$index.snp.type.png --group.col 2 --x.col 1   --y.col 3  --x.lab \"Cell Type\" --group.lab \"SNP annotation\" --y.lab \"Proportion\" --title.lab \" \" --legend.col 2  --width 6000 --no.grid\n";
print SH "Rscript  $Bin/dodgedBar_percent.r --infile $od/$index.SNP_mutation_num.xls  --outfile $od/$index.SNP_mutation_num.png  --group.col 2  --x.col 4   --y.col 3  --x.lab \"Cell Type\" --group.lab \"SNP number\" --y.lab \"Proportion\" --title.lab \" \" --legend.col 2  --width 6000 --no.grid\n";

`sh $od/$index.draw.sh`; 





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
