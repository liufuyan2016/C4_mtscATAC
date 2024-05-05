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
use Bio::SeqIO;


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
my ($od,$SNP,$index,$barplot_SNP,$bam);
GetOptions(
			"h|?" =>\&help,
			"od:s"=>\$od,
			"barplot_SNP:s"=>\$barplot_SNP,
			"bam:s"=>\$bam,
			"SNP:s"=>\$SNP,		
			"index:s"=>\$index,
			) || &help;
&help unless ($SNP && $od);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: 
    Usage:
        -bam              bam file
        -SNP              raw SNP file
	-barplot_SNP      barplot_SNP file
	-index            sample index
        -od               out dir   must be given
        -h                Help document
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
`mkdir $od`;
open(SH,">$od/filter.sh");
open(ID,"$SNP");
while(<ID>){
   chomp;
   my @tmp=split/\s+/;
    my $ts=$tmp[0];
   $ts=~s/>/_/;
   `mkdir $od/$ts`;
   open(OUT,">$od/$ts/id");
   print OUT "$tmp[0]\n";
   close OUT;
   print SH "perl $Bin/filter_SNP_from_mgatk.pl -bam $bam  -SNP $od/$ts/id  -barplot_SNP $barplot_SNP -od $od/$ts/\n";
}
close ID;
close SH;
`sh $od/filter.sh`;
`cat $od/*/final_select.SNP.txt  >$od/final_select.SNP.txt`;
`perl $Bin/draw_fastq_quality.pl -i $od/final_select.SNP.txt  -index $index  -o $od/`;
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
