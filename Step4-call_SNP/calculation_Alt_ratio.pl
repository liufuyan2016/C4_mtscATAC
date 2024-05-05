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
my ($in,$out,$index,$fa);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"index:s"=>\$index,
			"i:s"=>\$in,
			"fa:s"=>\$fa,			
			) || &help;
&help unless ($in && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: produce bed from sam
    Usage:
        -i          infile fq1    must be given
	-fa          mt fa
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
open (OUT, ">$out")|| die "cannot open $out:$!";
open(IN2,"$fa");
<IN2>;
while (<IN2>){
  chomp;
  my @tmp=split//,$_;
  for(my $i=0;$i<=$#tmp;$i++){
	my $pos=$i+1;
  	$hash{$pos}=$tmp[$i];
  }
 
}
close IN2;

open (IN, "<$in")|| die "cannot open $in:$!";
my $H=<IN>;
chomp $H;
my @seq=split/\,/,$H;
my $geneID;
while (<IN>){
        chomp;
	next if(/^$/);
        my @tmp=split/\,/;
        #my $main=$hash{$tmp[0]};
	my $T=0;
        for(my $i=1;$i<=$#tmp;$i++){
		$T+=$tmp[$i];
	}
        for(my $i=1;$i<=$#tmp;$i++){
		if($T!=0){
                my $ratio=$tmp[$i]/$T;
		push @{$hash2{$tmp[0]}},[($ratio,$seq[$i],$tmp[$i])];
		}
        }       
}
close IN;	
foreach my $key(keys %hash2){
     my @arr=@{$hash2{$key}};
     @arr=sort{$b->[0] <=>$a->[0]} @arr;
     if($hash{$key} eq "$arr[0][1]"){
 	    print OUT "$key\t$index\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\n";
	    print OUT "$key\t$index\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\n";
	    print OUT "$key\t$index\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\n";
	    print OUT "$key\t$index\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\n";
    }elsif($hash{$key} eq "$arr[1][1]"){
	    print OUT "$key\t$index\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\n";
	    print OUT "$key\t$index\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\n";
  	    print OUT "$key\t$index\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\n";
	    print OUT "$key\t$index\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\n";

   }elsif($hash{$key} eq "$arr[2][1]"){
	  print OUT "$key\t$index\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\n";
            print OUT "$key\t$index\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\n";
            print OUT "$key\t$index\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\n";
            print OUT "$key\t$index\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\n";
  }elsif($hash{$key} eq "$arr[3][1]"){
          print OUT "$key\t$index\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\n";
            print OUT "$key\t$index\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\n";
            print OUT "$key\t$index\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\t$arr[2][1]\t$arr[2][0]\t$arr[2][2]\n";
            print OUT "$key\t$index\t$arr[3][1]\t$arr[3][0]\t$arr[3][2]\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\n";
  }else{
	if($hash{$key} eq "$arr[4][1]"){
          print OUT "$key\t$index\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\t$arr[0][1]\t$arr[0][0]\t$arr[0][2]\n";
            print OUT "$key\t$index\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\t$arr[1][1]\t$arr[1][0]\t$arr[1][2]\n";
            print OUT "$key\t$index\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\t$arr[3][1]\t$arr[2][0]\t$arr[2][2]\n";
            print OUT "$key\t$index\t$arr[4][1]\t$arr[4][0]\t$arr[4][2]\t$arr[2][1]\t$arr[3][0]\t$arr[3][2]\n";
  }
}

}
close OUT;	







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
{		#��ȡָ��Ŀ¼���ļ��ľ���·��
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
