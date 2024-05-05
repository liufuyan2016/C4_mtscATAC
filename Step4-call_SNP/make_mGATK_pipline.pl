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
my ($out,$MT,$Dup,$SP,$input,$mix,$ID,$cell);
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"input:s"=>\$input,
			"mix"=>\$mix,
			"SP:s"=>\$SP,
			"Dup:s"=>\$Dup,
			) || &help;
&help unless ($input && $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: make mGTAK analysis 
    Usage:
        -input      input scATAC outdir
	-mix        for mix samples(human+mouse)
	-SP	     Species 
	-Dup	     keep Duplication for analysis
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
if(!defined $SP){
    $SP="mm10";
}
`mkdir $out`;
my @aln=glob("$input/*/aln.sam");

for(my $i=0;$i<=$#aln;$i++){
	my ($sample)=basename($input);
        print "$sample\n";
	my @xx=glob("$input/*/*.barcodeMerge.tsv");
	 my $merge=$xx[0];
        open(CELL,">$out/$sample.chrMt.cell.0");
	open (IN2, "<$merge")|| die "cannot open :$!"; 
	while(<IN2>) {
          chomp;
          next if(/^$/);
          my @tmp=split/\s+/;
	  $hash{$tmp[0]}=$tmp[1];
          print CELL "$tmp[1]\n";
        }
        close IN2;
        close CELL;
        my %hash_hg19=();
	my %hash_mm10=();
        my $erro_M_h;
        my $erro_h_M;
        my  $righ_hg19;
        my  $righ_mm10;
        `sort $out/$sample.chrMt.cell.0|uniq >$out/$sample.chrMt.cell`;
        open(OUT,">$out/$sample.chrMt.sam");
       if(defined $mix){
      		my $header=`head -500 $aln[$i]|grep "_chrM"|head -2`;
		 print OUT "$header";
                 print OUT "\@PG\tID:samtools\tPN:samtools\tVN:1.11CL\n";
                `grep -f $input/03.d2cfile/Human_barcode_list.txt $merge >$out/$sample.hg19.barcode`;
                `grep -f $input/03.d2cfile/Mouse_barcode_list.txt $merge >$out/$sample.mm10.barcode`;
	        open(H,"$out/$sample.hg19.barcode");    
                while(<H>){
                      chomp;
                      my @tmp=split/\s+/;
                      $hash_hg19{$tmp[0]}=$tmp[1]; 
                }
                close H;
                open(M,"$out/$sample.mm10.barcode");                    
                while(<M>){
                      chomp;
                      my @tmp=split/\s+/;
                      $hash_mm10{$tmp[0]}=$tmp[1];
                }
                close M;

      }else{
            	  my $header=`head -500 $aln[$i]|grep "chrM"|head -1`;
		  print OUT "$header";
                  print OUT "\@PG\tID:samtools\tPN:samtools\tVN:1.11CL:hg19_chrM mm10_chrM\n";
     }
       open (IN, "<$aln[$i]")|| die "cannot open :$!";
       while (<IN>) {
	  chomp;
	  next if(/^$/);
	  my @tmp=split/\s+/;
	  if(defined $tmp[11]){
          	my ($mis)=$tmp[11]=~/NM:i:(\d+)/; 
		next if($mis>2);
          }
	  if($tmp[2]=~/chrM/){
		my ($barcode)=$tmp[13]=~/CB\:Z\:(.*)/;	
		next if(!defined $hash{$barcode});
		$tmp[13]="CB\:Z\:$hash{$barcode}";
		print OUT join("\t",@tmp)."\n";
                 if(defined $mix){
                      if($tmp[2]=~/hg19/&& defined  $hash_mm10{$barcode}){
                            $erro_h_M++;
                       }elsif($tmp[2]=~/mm10/&& defined  $hash_hg19{$barcode}){
                            $erro_M_h++;
                       }
             	       if($tmp[2]=~/hg19/){
                            $righ_hg19++;
                       }else{
                            $righ_mm10++;
                       }
         	}              

	  }
	
      }
     close IN;
     close OUT;
      if(defined $mix){
                my $ratiohg19=$erro_h_M/$righ_hg19*100;
                my $ratiomm10=$erro_M_h/$righ_mm10*100;
                print "error hg19\t$erro_h_M\t$righ_hg19\t$ratiohg19\n";
                print "error mm10\t$erro_M_h\t$righ_mm10\t$ratiomm10\n";
		open(SH,">$out/mgatk.sh");
        }

     if(defined $mix){
     `grep "hg19_" $out/$sample.chrMt.sam |grep -v "mm10_"  >$out/$sample.hg19chrMt.sam`;
     `grep "mm10_" $out/$sample.chrMt.sam |grep -v "hg19_"  >$out/$sample.mm10chrMt.sam`;
     `samtools view -bS -o  $out/$sample.hg19Mt.bam $out/$sample.hg19chrMt.sam`;
     if(defined $Dup){
            print SH "cd $out &&source /home/kangjingmin/venv3/bin/activate && mgatk bcall -kd  -i  $sample.hg19Mt.bam   -n hg19 -o $sample.hg19 -c  20 -bt CB -b $sample.chrMt.cell  -g $Bin/hg19_mix.fasta";
     }else{
	  print SH "cd $out &&source /home/kangjingmin/venv3/bin/activate && mgatk bcall  -i  $sample.hg19Mt.bam   -n hg19 -o $sample.hg19 -c  20 -bt CB -b $sample.chrMt.cell  -g $Bin/hg19_mix.fasta";
     }
     print SH "&&samtools view -bS -o  $sample.mm10Mt.bam $sample.mm10chrMt.sam";
     if(defined $Dup){
	     print SH "source /home/kangjingmin/venv3/bin/activate && mgatk bcall  -kd -i  $sample.mm10Mt.bam   -n mm10  -o $sample.mm10 -c  20 -bt CB -b $sample.chrMt.cell  -g $Bin/mm10_mix.fasta";
      }else{
          print SH "&&cd $out &&source /home/kangjingmin/venv3/bin/activate && mgatk bcall  -i  $sample.hg19Mt.bam   -n hg19 -o $sample.hg19 -c  20 -bt CB -b $sample.chrMt.cell  -g $Bin/hg19_mix.fasta";
        }
    }else{
      `samtools view -bS -o  $out/$sample.${SP}Mt.bam  $out/$sample.chrMt.sam`;
       if(defined $Dup){
      		`cd $out &&source /home/kangjingmin/venv3/bin/activate && mgatk bcall -kd  -i $sample.${SP}Mt.bam -n mm10  -o $sample.$SP -c 90 -bt CB -b  $sample.chrMt.cell  -g $Bin/$SP\_add.fasta`;
	}else{
		`cd $out &&source /home/kangjingmin/venv3/bin/activate && mgatk bcall  -i $sample.${SP}Mt.bam -n mm10  -o $sample.$SP -c 90 -bt CB -b  $sample.chrMt.cell  -g $Bin/$SP\_add.fasta`;
	}
   }	
}

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
