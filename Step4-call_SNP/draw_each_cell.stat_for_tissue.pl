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
        -cell        cell anno
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
open(SH,">$out/draw.sh");
open(OUT,">$out/draw.mt_ratio.txt");
open(OUT2,">$out/draw.mt_count.txt");
open(OUT3,">$out/draw.nuc.count.txt");

open(IN,"$cell");
while(<IN>){
   chomp;
   my @tmp=split/\s+/;
   my ($type)=$tmp[0]=~/\#(.*)/;
   $hash{$type}=$tmp[-1];
}
close IN;
open(IN2,"$input");
while(<IN2>){
   chomp;
   my @tmp=split/\s+/;
   next if(!defined $hash{$tmp[2]});
   my $x=$hash{$tmp[2]};
   if($tmp[1]=~/Y_GS/){
	$tmp[1]="Young Bone marrow";
    }else{
	$tmp[1]="Old Bone marrow";
    }
   next if($x eq "C1");
    print OUT join("\t",@tmp)."\t$x\n";
	 $x=~s/C//;
   
   push @{$hash2{$tmp[2]}{$x}},[@tmp];
   #print OUT join("\t",@tmp)."\n";
}
close IN2;
die;
#print Dumper %hash2;
foreach my $key(keys %hash2){
    print OUT "$key";
    print OUT2 "$key";
    print OUT3 "$key";
    foreach my $key2(sort {$a<=>$b} keys %{$hash2{$key}}){
          print OUT "\t$key2";
	  print OUT2 "\t$key2";
	  print OUT3 "\t$key2";
    }
            print OUT "\n";
        print OUT2 "\n";
        print OUT3 "\n";
    last;
}
	
foreach my $key(keys %hash2){
    print OUT "$key\t";
    print OUT2 "$key\t";
    print OUT3 "$key\t";
    foreach my $key2(sort {$a<=>$b} keys %{$hash2{$key}}){
	my @arr=@{$hash2{$key}{$key2}};
	print OUT "\t$arr[0][3]";
	print OUT2 "\t$arr[0][4]";
	print OUT3 "\t$arr[0][5]";
	}
	print OUT "\n";
	print OUT2 "\n";
	print OUT3 "\n";
}

close OUT;
close OUT2;
close OUT3;
#print SH  "Rscript $Bin/box_violin.r  --infile $out/draw.txt --group.col 7   --outfile $out/$key.MT.NUC.cell.stat.MT.png --value.col 4 --x.col 2 --x.lab \"Sample \" --y.lab \" MT ratio \" --title.lab \"  \"  --no.grid\n";
 #  print SH "Rscript $Bin/box_violin.r  --infile $out/draw.txt  --group.col 7 --outfile $out/$key.MT.NUC.cell.stat.MT.fragment.png --value.col 5 --x.col 2 --x.lab \"Sample\" --y.lab \"log10 mtDNA complexity\" --title.lab \"  \"  --no.grid\n";
 # print SH "Rscript $Bin/box_violin.r  --infile $out/draw.txt  --group.col 7 --outfile $out/$key.MT.NUC.cell.stat.NUC.fragment.png --value.col 6 --x.col 2 --x.lab \"sample\" --y.lab \"log10 chromatin complexity\" --title.lab \"  \"  --no.grid\n";


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
