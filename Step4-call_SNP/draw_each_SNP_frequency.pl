#!/usr/bin/perl -w
# Copyright (c) BGI  2015/6/28
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
my ($in,$out,$len,$order,@sc,@bulk,$index,$fa);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"index:s"=>\$index,
			"i:s"=>\$in,
			"sc:s{,}"=>\@sc,
		#	"sp:s"=>\$sp,
			"order"=>\$order,
			"len:i"=>\$len,
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
        -i           infile *.mt_SNP.txt or SNP order file   must be given
	-sc          C4_mtscATAC  mt frequency
	-order       use the SNP order
        -len         Mt length
	-index       index file
        -o           outfile    must be given
        -h           Help document
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
my %hashDepth;
my @draw_order;
my %hashDepthN;
my %hashDD;
my %hashSS;
open (OUT2, ">$out.all")|| die "cannot open $out:$!";
open (OUT, ">$out")|| die "cannot open $out:$!";
for(my $i=0;$i<=$#sc;$i++){
	open(IN2,"$sc[$i]");
	while (<IN2>){
  		chomp;
  		my @tmp=split/\s+/,$_;
		next if($tmp[0]>$len);
  		$hash{$tmp[0]}{$tmp[5]}+=$tmp[6];
		#$hashDepthN{$i}{$tmp[0]}++;
		#$hashDepth{$i}{$tmp[0]}+=($tmp[4]/$tmp[3]);
		$hashDepth{$i}{"$tmp[0]\t$tmp[5]"}=$tmp[-1];
		$hashDD{$tmp[0]}{$tmp[5]}+=$tmp[-1];
                $hashSS{$tmp[0]}=$tmp[2];

	}
	close IN2;
}
my $N=0;
my %hash_order;
my %hash_result;
open (IN, "<$in")|| die "cannot open $in:$!";
while (<IN>){
        chomp;
	next if(/^$/);
        my @tmp=split/\s+/;
	if(defined $order){
		push @draw_order,$tmp[0];
	}
	next if($tmp[0]!~/\d+[TCGA]>[TCGA]/);
	my ($pos,$seq,$SNP)=$tmp[0]=~/(\d+)([TCGA])>([TCGA])/;
	next if($pos>$len);
	$tmp[0]="$pos$seq>$SNP";
	my $sc_ratio=$hash{$pos}{$SNP}/($#sc+1);
	my $T;
	foreach my $key(keys %hashDepth){
		$T+=$hashDepth{$key}{"$pos\t$SNP"};
	}
	$hash_result{$tmp[0]}{"mtscATAC"}="$sc_ratio\t$T";
	$hash_order{$tmp[0]}=$sc_ratio;
}
close IN;
my %hashOrder;
my $M=0;
if(defined $order){
	@draw_order=reverse(@draw_order);
	for(my $i=0;$i<=$#draw_order;$i++){
		my $key=$draw_order[$i];
		 print OUT "$key\t$hash_result{$key}{\"mtscATAC\"}\tmtscATAC\n";
	}
}else{
	foreach my $key(sort{$hash_order{$a}<=>$hash_order{$b}} keys %hash_order){
		print OUT "$key\t$hash_result{$key}{\"mtscATAC\"}\tmtscATAC\n";
                $hashOrder{$key}=$M++;
	}
}

foreach my $key(keys %hash){
    #if($sp=~/human/i){
    #       next if($key>16569);
    #}
    foreach my $key2(keys %{$hash{$key}}){
           next if($key2 eq "N");
           my $sc_ratio=$hash{$key}{$key2}/($#sc+1);
           my $T1=$hashDD{$key}{$key2};
           print OUT2 "$key$hashSS{$key}>$key2\t$sc_ratio\t$T1\tmtscATAC\n";
    }
}
close OUT2;
open(OUT9,">$out.order");
foreach my $key(sort{$hashOrder{$b}<=>$hashOrder{$a}} keys %hashOrder){
     print OUT9  "$key\n";
}
close OUT9;

print "Rscript $Bin/group_Bar.r  --infile $out --outfile $out.png --group.col 4 --x.col 1 --y.col 2 --group.lab \" \" --x.lab \"SNPs\" --y.lab \"VAF\" --no.grid  --height 6500 --width 2000 --axis.size 9  --lab.size 12 --legend.size 7\n";
`Rscript $Bin/group_Bar.r  --infile $out --outfile $out.png --group.col 4 --x.col 1 --y.col 2 --group.lab " " --x.lab "SNPs " --y.lab "VAF" --no.grid  --height 4755  --width 2000 --axis.size 9  --lab.size 12 --legend.size 7`;
`Rscript $Bin/group_Bar.r  --infile $out --outfile $out.depth.png --group.col 4 --x.col 1 --y.col 3 --group.lab " " --x.lab "SNPs " --y.lab "Depth" --no.grid  --height 4755  --width 2000 --axis.size 9  --lab.size 12 --legend.size 7`;




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
