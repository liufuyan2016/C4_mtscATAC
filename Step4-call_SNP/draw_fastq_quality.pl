#!/usr/bin/perl -w
# Copyright (c) BGI 2015/6/28
# Writer:         liufuyan <liufuyan@biomarker.com.cn>
# Program Date:   2015/6/28.
# Modifier:       liufuyan <liufuyan@biomarker.com.cn>
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
my $Writer = "liufuyan <liufuyan\@genomcis.cn>";
my $Data   = "2015/6/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$out,$order,$sp,@sc,@bulk,$index,$fa);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"index:s"=>\$index,
			"i:s"=>\$in,
			"sc:s{,}"=>\@sc,
			"sp:s"=>\$sp,
			"order:s"=>\$order,
			"bulk:s{,}"=>\@bulk,			
			) || &help;
&help unless ($in && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function:
    Usage:
        -i          infile order file   must be given
	-index      index file
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
my %hash2;
my %hashT1;
my %hashT2;
my %hashTN1;
my %hashTN2;
my %hashDepth;
my %hashDepth2;
my %hashDD;
my %hashSS;
if(-s "$out/$index.mis_Q.txt"){
  `rm $out/$index.mis_Q.txt`;
}
open(IN,"$in");
while (<IN>){
      chomp;
      my @tmp=split/\s+/,$_;
      my $dir=$tmp[0];
      $dir=~s/>/_/;
     `cat $out/$dir/mis_Q.txt >>$out/$index.mis_Q.txt`;
}
close IN;
`Rscript $Bin/fastq_quality.r   --infile $out/$index.mis_Q.txt  --outfile $out/$index.mis_Q.png --value.col 2 --x.col 1 --x.lab " " --y.lab "Base quality score" --title.lab " "  --no.grid  --height 3755 --width 1000`;

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub reverse_complementary{
	my ($inputseq)=@_;
	my $sequence = reverse(substr($inputseq,0));
	$sequence =~ tr/ACGTUMRWSYKVHDBNacgtumrwsykvhdbnn/TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/d;
	$inputseq=substr($sequence, 0);
	return $inputseq;
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
