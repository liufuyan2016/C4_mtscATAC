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
my ($od,$SNP,$barplot_SNP,$bam);
GetOptions(
			"h|?" =>\&help,
			"od:s"=>\$od,
			"barplot_SNP:s"=>\$barplot_SNP,
			"bam:s"=>\$bam,
			"SNP:s"=>\$SNP,		
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
	-barplot_SNP     barplot_SNP file
        -od         out dir   must be given
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
my %hashAlt;
my %hashSNP;
my %hashStrand;
`mkdir $od`;
open(ID,"$SNP");
while(<ID>){
   chomp;
   my @tmp=split/\s+/;
   $hashSNP{$tmp[0]}=1;
}
close ID;
my %hashDepth;
my %hashVAF;
####Only use highest Alt for VAF lower to avoid the sequecing error
open(IN,"$barplot_SNP");
while(<IN>){
  chomp;
  my @tmp=split/\s+/;
  next if($tmp[-1] ne "mtscATAC");
  my ($pos)=$tmp[0]=~/(\d+)/;
  $tmp[1] = sprintf("%.8f", $tmp[1]);
  push @{$hash{$pos}},[($tmp[1],$tmp[2],$tmp[0])];
  $hashVAF{$tmp[0]}=$tmp[1];
}
close IN;
foreach my $key(keys %hash){
    my @arr=@{$hash{$key}};
    @arr=sort{$b->[0]<=>$a->[0]} @arr;
    if($arr[0][1]>30){
	if($arr[1][0] eq "0.00000000"){  
	     $hashAlt{$arr[0][2]}=99999;
        }elsif($arr[0][0]/$arr[1][0]>5){
		$hashAlt{$arr[0][2]}=$arr[0][0]/$arr[1][0];
	}		
    
     }
}
`mkdir $od/SAM`;
#######read merge sam file to judge strand information and base quality
my %hashQ;
my $H=`samtools view -H $bam`;
`samtools view $bam |grep -E  "NM:i:1|NM:i:2"  -w  >$od/SNP.sam`;
open(Q,">$od/mis_Q.txt");
my $AltFreq=0;
foreach my $key(keys %hashSNP){
  if(!defined  $hashAlt{$key}){
	$AltFreq++;
	 print "$key\tAl]t\n";
         next;
   }
  my ($d,$seq)=$key=~/(\d+)[ATCG]>([ATCG])/;
  my $id=$key;
  $id=~s/>/_/;
   open(OUT ,">$od/SAM/$id.sam");
   print OUT "$H";
   open(IN,"$od/SNP.sam");
   while(<IN>){
     chomp;
     my @tmp=split/\s+/;
     my $s1=$d-$tmp[3];
     my $len=length($tmp[9]);
     next if($s1<0||$s1>=$len);
     $s1=&find_snps($tmp[5],$tmp[3],$d); 
     if(!defined $s1){
	print "$_\n";
     }
     my $q1=substr($tmp[10],$s1,1);
     my $seq1=substr($tmp[9],$s1,1);
     $q1=ord($q1)-33;
     push @{$hashQ{$key}},$q1;
     print  Q "$key\t$q1\n";
     next if($seq1 ne $seq);
     my $strand=&Flagstat($tmp[1]);
     $hashStrand{$key}{$strand}++;
     print OUT  "$_\n";
   }
  close IN;
  close OUT;
  close Q;
  `samtools view -bS -o $od/SAM/$id.bam $od/SAM/$id.sam`;
  `samtools index $od/SAM/$id.bam`;
  `rm $od/SNP.sam`;

}
open(OUT ,">$od/final_select.SNP.txt");
my $LowQ=0;
my $strand=0;
foreach my $key(keys %hashStrand){
    my @Quality=@{$hashQ{$key}};
    my $sum=0;
    foreach my $num (@Quality) {
 	   $sum += $num;
    }
    my $avg = $sum / scalar(@Quality);
    my $dep=$hashStrand{$key}{"+"}+$hashStrand{$key}{"-"};
    my $S1=$hashStrand{$key}{"+"}/($hashStrand{$key}{"+"}+$hashStrand{$key}{"-"});
    my $S2=$hashStrand{$key}{"-"}/($hashStrand{$key}{"+"}+$hashStrand{$key}{"-"});
    if($avg<20){
	  $LowQ++;
	  print "LowQ\t$key\t$avg\n";
	  next;
     }elsif($S1>0.9||$S2<0.1){
		$strand++;
		print "strand\t$key\t$hashStrand{$key}{\"+\"}\+\|$hashStrand{$key}{\"-\"}\-\n";
		next;
      }elsif($dep <=30){
		print "strand\t$key\t$dep\n";
		 next;
      }
    print OUT  "$key\t$hashVAF{$key}\t$avg\t$S1\t$S2\t$hashAlt{$key}\t$dep\n";
}
close OUT;












###############Time
sub find_snps{
  my ($cigar_string, $read_start,$SNPpos) = @_;
  
  my $ref_pos = $read_start;
  my $query_pos = 0; 
  my @snp_positions; 
  my $dd=$SNPpos-$read_start; 
  my $SNP;
  while ($cigar_string =~ /(\d+)([mid])/ig) {
    my $op_len = $1;
    my $op_code = $2;
    print "$op_len\t$op_code\t$query_pos\n";
    if($op_code eq 'M') { 
      $ref_pos += $op_len;
      $query_pos += $op_len;
      if($ref_pos>$SNPpos) {
	    my $d2=$ref_pos-$SNPpos;
            $SNP=$query_pos-$d2;
	    last;
      }
   }elsif ($op_code eq 'I') { 
      $query_pos += $op_len;
   }elsif ($op_code eq 'D') {
      $ref_pos += $op_len;
   }
  }
   #my $seq="GTGTATAAGAGACAGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCC";
  #my $Q=substr($seq,$SNPpos,2);
  #print "$Q\n"; 
    return $SNP;
}


sub Flagstat{
 my ($flag) =@_;
my %flag_values = (
    '1'  => ['PAIRED',                 'paired-end read'],
    '2'  => ['PROPER_PAIR',            'properly paired read'],
    '4'  => ['UNMAPPED',               'read is unmapped'],
    '8'  => ['MATE_UNMAPPED',          'mate is unmapped'],
    '16' => ['REVERSE_STRAND',         'read is on the reverse strand'],
    '32' => ['MATE_REVERSE_STRAND',    'mate is on the reverse strand'],
    '64' => ['FIRST_IN_PAIR',          'read is the first read in a pair'],
    '128'=> ['SECOND_IN_PAIR',         'read is the second read in a pair'],
    '256'=> ['NOT_PRIMARY_ALIGNMENT', 'read is not primary alignment'],
    '512'=> ['QC_FAILED',              'read failed quality control'],
    '1024'=> ['DUPLICATE',             'read is a PCR duplicate'],
);
my $result;
 my @flags;
 foreach my $bit (reverse split //, sprintf("%012b", $flag)) {
     unshift @flags, $bit;
     }
     my $flag_string = join('', @flags);

     #print "Flag value: $flag\n";
     #print "Flag string: $flag_string\n";
     foreach my $bit (sort {$a <=> $b} keys %flag_values) {
         my ($name, $description) = @{$flag_values{$bit}};
             if ($flag & $bit) {
                     if($description=~/read is on the reverse strand/){
			  $result="-";
			   last;
		     }elsif($description=~/mate is on the reverse strand/){
			  $result="+";
			  last;
		     }
              }
      }
   return $result;
  }

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
