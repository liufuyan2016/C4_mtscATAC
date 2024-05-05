#!/usr/bin/perl -w
# Copyright (c) BGI 2015/6/28

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
my $Writer = "liufuyan <liufuyan\@genomics.com>";
my $Data   = "2015/6/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$continue,$MT,$out,$cirRNA,$all,@AUC,$GSVA,$p,$complex,$tpm,$group,$N);
$p=1;
my $split_row;
my $order=0;
my $width=2000;
my $height=4000;
GetOptions(
			"h|?" =>\&help,
			"o:s"=>\$out,
			"i:s"=>\$in,
			 "N:i"=>\$N,
			"MT"=>\$MT,
			"GSVA"=>\$GSVA,
			"continue"=>\$continue,
			 "p:f"=>\$p,
			"cirRNA:s"=>\$cirRNA,
			 "AUC:s{,}"=>\@AUC,
			 "width:i"=>\$height,
			 "height:i"=>\$width,
			"split_row"=>\$split_row,
			"complex"=>\$complex,			
			 "cluster:f"=>\$order,
			"group:s"=>\$group,
			"tpm:s"=>\$tpm,
			"all"=>\$all,
			) || &help;
&help unless ($in && $N&& $out);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: only get Up DEGs and make drawing top N genes according to the fold; if give -GSVA is used for GSVA analysis and -i should be input U-test results.
    Usage:
        -i           infile DEG
	-all         for up and all DEGs. default is use up DEG, if use all ,please give this options
	-tpm	     expression  file
	-AUC	     AUC value and/or file, eg 0.85 
	-continue    for continue blue-red heatmap
	-split_row    split row according to the -i input  file 2 col  
	-p	     p value
	-MT 	     for C4_mtscATAC analysis
	-cirRNA      give VS file for cirRNA annotation 
	-GSVA        for GSVA picture
	-N	     top N number
	-group       group file
	-cluster       draw use input order,set 1 is for cluster
	-complex     for complex file
	-width	     width for picture,default is 4000
	-height      height for picture,default is 2000
        -o           outfile    must be given
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
if(-s "$out"){
	`rm -r $out`;
}
`mkdir $out`;
my $outdir=$out;
my $dirN=basename($out);
$out="$out/$dirN";

open (OUT, ">$out")|| die "cannot open $out:$!";
my %hash2;
my %hash;
my %hash3;
$hash2{"\#ID"}=1;
$hash3{"\#ID"}=1;

my %hash_group;
my %hash_group2;
my $Num=0;
my %hash_order;
open(G,"$group");
while(<G>){
  chomp;
  next if($_=~/Sample/);
  my @tmp=split/\s+/;
  if($#tmp>1){
         $hash_group2{$tmp[0]}=$tmp[2];
  }else{
         $hash_group2{$tmp[0]}=$tmp[1];
  }
  $hash_order{$tmp[1]}=$Num++;
  push @{$hash_group{$tmp[1]}},$tmp[0];
}
close G;
my %hash_AUC;
if(defined $AUC[0]){
	if(defined $AUC[1]){
		;
	}else{
		`perl $Bin/ROC/ROC_from_DEG_file.pl -i $in  -od $outdir/ROC`;
		$AUC[1]="$outdir/ROC/AUC.stat.xls";
	} 			
	open(IN2, "<$AUC[1]")|| die "cannot open $AUC[1]:$!";
	<IN2>;
	while (<IN2>){
		chomp;
		my @tmp=split/\s+/;
		if($tmp[1]<$AUC[0]){
			 $hash_AUC{$tmp[0]}=0;
		}else{
			$hash_AUC{$tmp[0]}=$tmp[1];
		}
	}
	close IN2;
}
my @fold_array;
open(SPLIT,">$out.gene_split");
open (IN, "<$in")|| die "cannot open $in:$!";
while (<IN>) 
{
	chomp;
	next if(/^$/);
	next if($_=~/\#/);
	my @tmp=split/\s+/;
	next if($tmp[0]=~/MSTRG\./);
	next if($tmp[0]=~/unconservative_/);
	next if($tmp[0]=~/_newGene/);
	next if($tmp[0]eq "x");
	if(defined $split_row){
		print SPLIT "$tmp[1]\n";
		$order=0;
	}
=cut
	if($tmp[-2] eq "Inf"){
		$tmp[-2]=&average(@tmp);
	}elsif($tmp[-2] eq "-Inf"){
		$tmp[-2]=&average(@tmp);
		$tmp[-2]=$tmp[-2]*-1;
        }
=cut
	if(defined $GSVA){
		next if($_=~/log2FC/);
		next if($tmp[-3]>$p);
		$hash{$tmp[0]}{"fold"}=$tmp[-3];
	        $hash{$tmp[0]}{"p"}=$tmp[-3];
		push @fold_array,[($tmp[0],$tmp[-3])];
	}else{
		if($#tmp<=4){
			if(defined $AUC[0]){
				$hash{$hash_AUC{$tmp[0]}}{$tmp[0]}{"fold"}=0.001;
                        	$hash{$hash_AUC{$tmp[0]}}{$tmp[0]}{"p"}=0.001;
			}else{
				$hash{$tmp[0]}{"fold"}=0.001;
                		$hash{$tmp[0]}{"p"}=0.001;
			}
			push @fold_array,[($tmp[0],0)];
		}else{
			my $fold;
			if(!defined $all){
			       next if($tmp[-1] ne  "up");
			        $fold=$tmp[-2];
			 }else{
			      $fold=abs($tmp[-2]);
			  }	
			 next if($tmp[-3]>$p);
			  push @fold_array,[($tmp[0],$tmp[-2])];
			 if(defined $AUC[0]){
				$hash{$hash_AUC{$tmp[0]}}{$tmp[0]}{"fold"}=$fold;
				$hash{$hash_AUC{$tmp[0]}}{$tmp[0]}{"p"}=$tmp[-3];
			 }else{
				$hash{$tmp[0]}{"fold"}=$fold;
				$hash{$tmp[0]}{"p"}=$tmp[-3];
			 }
		}
	}
}
close SPLIT;
my $num=0;
if(defined $AUC[0]){
	foreach  my $key(sort{$b<=>$a} keys %hash) {
		foreach  my $key2(sort {$hash{$key}{$b}{"fold"}<=>$hash{$key}{$a}{"fold"}}keys %{$hash{$key}}) {
        	$num++;
        	last if($num>$N);
		next if($hash_AUC{$key2}==0);
        	$hash2{$key2}=1;
		}
	}
	$num=0;
	foreach  my $key(sort{$b<=>$a} keys %hash) {
		foreach  my $key2(sort {$hash{$key}{$a}{"p"}<=>$hash{$key}{$b}{"p"}}keys %{$hash{$key}}) {
     		   $num++;
        	  last if($num>$N);
		  next if($hash_AUC{$key2}==0);
      		  $hash3{$key2}=1;
		}
	}
}else{
	foreach  my $key(sort {$hash{$b}{"fold"}<=>$hash{$a}{"fold"}}keys %hash) {
	$num++;
	last if($num>$N);
	$hash2{$key}=1;
	}
	$num=0;
	foreach  my $key(sort {$hash{$a}{"p"}<=>$hash{$b}{"p"}}keys %hash) {
        $num++;
        last if($num>$N);
        $hash3{$key}=1;
	}
}
my %hashVS;
if(defined $cirRNA){
	open(VS,$cirRNA);
}else{
	open(VS,"$Bin/Symbol_2_ID.list");
}
while(<VS>){
  chomp; 
  next if($_=~/\#/);
  my @tmp=split/\s+/;
  if(defined $cirRNA){
	next if($tmp[2]<90||(defined $tmp[3]&&$tmp[3]<90));
	my ($name)=$tmp[1]=~/^([^\|]+)/;
	$hashVS{$tmp[0]}=$name;
  }else{
  	$hashVS{$tmp[1]}=$tmp[0];
  }
}
close VS;
my $flag=0;
my %hash_express;
open (OUT2, ">$out.p")|| die "cannot open $out:$!";
open (IN2, "<$tpm")|| die "cannot open $tpm:$!";
while (<IN2>)
{
        chomp;
        next if(/^$/);
        my @tmp=split/\s+/;
	if($tmp[0]eq "ID"){
		$flag=1;
		$tmp[0]="#ID";
	}
	if($flag==1){
		pop @tmp;
		pop @tmp;
		pop @tmp;
	}		
         if(defined $hashVS{$tmp[0]}){
		my $key=$tmp[0];
		$tmp[0]=$hashVS{$tmp[0]};
		$hash_express{$key}= join("\t",@tmp);
	  }else{
		my $key=$tmp[0];
                $hash_express{$key}= join("\t",@tmp);
		}

}

close IN2;
print OUT "$hash_express{\"\#ID\"}\n";
print OUT2 "$hash_express{\"\#ID\"}\n";
my @fold_array_sort=sort{$b->[1]<=>$a->[1]} @fold_array;
for(my $i=0;$i<=$#fold_array_sort;$i++){
  if(defined $hash2{$fold_array_sort[$i][0]}){
		 print OUT "$hash_express{$fold_array_sort[$i][0]}\n";
   }
  if(defined $hash3{$fold_array_sort[$i][0]}){
     print OUT2 "$hash_express{$fold_array_sort[$i][0]}\n";
  }
}

my ($heatmap1,$groupall1,$lable1,$heatmap2,$groupall2, $lable2);
if(defined $complex){
   foreach my $key(sort{$hash_order{$a}<=>$hash_order{$b}} keys %hash_order){
	open(G2,">$out.group$key");
	open(G1,">$out.group$key.drawing");
	my @array=@{$hash_group{$key}};
	foreach my $list (@array){
		print G1 "$list\t$hash_group2{$list}\n";
	}
	close G1;
	 foreach my $list (@array){
                print G2 "$list\t$key\n";
        }
        close G2;
	$heatmap1.="$out.$key.txt,";
	$groupall1.="$out.group$key.drawing,";
	$lable1.="$key,";
        $heatmap2.="$out.$key.p.txt,";
        $groupall2.="$out.group$key.drawing,";
        $lable2.="$key,";
	print "perl $Bin/fpkm_prepare_accord_group.pl  -i $out  -group $out.group$key   -o  $out.$key.txt -order $key\n";
	`perl $Bin/fpkm_prepare_accord_group.pl  -i $out  -group $out.group$key   -o  $out.$key.txt -order $key `;
	`perl $Bin/fpkm_prepare_accord_group.pl  -i $out.p  -group $out.group$key   -o  $out.$key.p.txt -order $key `;
	}
}
my $log;
if(defined $GSVA){
	$log=-1;
}else{
	$log=1;
	if(defined $split_row){
		 $log.=" $out.gene_split  ";
	}
}
my $sampleN=keys %hash_group2;
$height=$sampleN*85+1400;
my $script;
if(defined $MT){
   	$script="ComplexHeatmap_Col.r";
}elsif(defined $split_row){
	$script="ComplexHeatmap_Col_similar_zscore.r4 ";
}elsif(defined $continue){
	$script="ComplexHeatmap_Col_similar_zscore.r3";
}else{
	if($order==0){
		  $script="ComplexHeatmap_Col_similar_zscore.r1";
	}else{
   		  $script="ComplexHeatmap_Col_similar_zscore.r0";
        }
}
print "Rscript $Bin/$script  $heatmap1   $out  $groupall1  $lable1 $width $height 0.95 $log\n";
if(defined $continue){
	`Rscript $Bin/$script  $heatmap1   $out.90  $groupall1  $lable1 $width $height  0.50 $log`;
}else{
	`Rscript $Bin/$script  $heatmap1   $out.90  $groupall1  $lable1 $width $height  0.90 $log`;
	`Rscript $Bin/$script  $heatmap1   $out.95  $groupall1  $lable1 $width $height 0.95 $log`;
	`Rscript $Bin/$script  $heatmap1   $out.85  $groupall1  $lable1 $width $height  0.85 $log`;
	`Rscript $Bin/$script  $heatmap1   $out.50  $groupall1  $lable1 $width $height  -1 $log`; 
	`Rscript $Bin/$script  $heatmap1   $out.80  $groupall1  $lable1 $width $height  0.8  $log`;
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
