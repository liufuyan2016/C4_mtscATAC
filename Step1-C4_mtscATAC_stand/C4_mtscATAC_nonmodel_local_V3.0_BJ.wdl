##################################################################################
# source: /jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Automated/USER/songyumo/pipeline/scATAC/src/C4scATAC_nonmodel_local_V2.1.wdl
###############################################################################################
workflow iDrop_scATAC_default_2 {
  Array[Array[String]] Data
  String Outdir
  String SampleID
  String ProjectID
  String reference
  String? readStructure
  String outdir=Outdir
  String ID=SampleID
  String PID = "st.q -P "+ProjectID+"_sentieon"

##########    genome infor   ############
 Map[String,Array[String]] refConfig = {
    "mt_mm10" : ["mm10","2871842584","Genome/mm10/mm10_mask_MT.fa","Genome/mm10/tss.bed","Genome/mm10/promoter_scATAC.bed","chrM","Genome/mm10/mm10_mask_MT.fa.index","Genome/mm10/mm10_chrom_size.tsv","None"],    
    "mt_hg19" : ["hg19","3137454505","Genome/hg19/hg19_mask_MT.fa","Genome/hg19/regions/tss.bed","Genome/hg19/promoter_scATAC.bed","chrM","Genome/hg19/hg19_mask_MT.fa.index","Genome/hg19/hg19_chrom_size.tsv","None"]
 }


########################################
 String refdir=refConfig[reference][2]
  String Rscript="Rscript"
  String root="/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Pipline/DNBelabC4scATAC/"##download from https://github.com/M-wen/C4_scATAC_analysis
  String tss=refConfig[reference][3]
  String promo=refConfig[reference][4]
  String chrmt=refConfig[reference][5]
  String refindex=refConfig[reference][6]
  String whitelist="whitelist.txt"
  String  chromeSize=refConfig[reference][7]
  String blacklist=refConfig[reference][8]
  String runID=ID
  String? defaultConfig="C4scATAClib_seqT1_R1_70_R2_50.json"
  String config=select_first([readStructure,defaultConfig])
  String lib="/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Pipline/DNBelabC4scATAC/script/lib.sh"##set your enviroment
  String species=refConfig[reference][0]
  String macs2sp=refConfig[reference][1]
  String d2cpara = if blacklist != "None" then "--bg "+chromeSize+" --ts "+tss+" --bl "+blacklist else "--bg "+chromeSize+" --ts "+tss
  Int ?ForceFrag

  String read_format = if config == "oldT1" then "bc:50:-1,r2:0:49" else if config == "500" then "bc:50:59,bc:66:-1,r2:0:49" else "bc:0:19,r1:20:-1"
  Array[String] barcode = if config == "oldT1" then fq2 else if config == "500" then fq2 else fq1


  call makedir {
    input:
    Dir=outdir,
  }
  scatter(fq in Data){
    String fq1=fq[0]
    String fq2=fq[1]
  }

 call mapping {
  input:
    ID=ID,
    fastq1=fq1,
    fastq2=fq2,
    root=root,
    outdir=makedir.Outdir,
    refdir=refdir,
    barcode=barcode,
    read_format=read_format,
    refindex=refindex,
    whitelist=whitelist,
 
  }
  call deconvolution{
    input:
    lib=lib,
    alnbed=mapping.alnbed,
    outdir=outdir,
    root=root,
    ID=ID,
    d2cpara=d2cpara,
    species=species,
    ForceFrag=ForceFrag,
    chrmt=chrmt,
}

  call qc{input:lib=lib,root=root,tss=tss,FragmentFile=deconvolution.FragmentFile,outdir=outdir,ID=ID,Rscript=Rscript}
  call peakcount{input:lib=lib,chrmt=chrmt,root=root,listtxt=deconvolution.listtxt,macs2sp=macs2sp,promo=promo,outdir=outdir,FragmentFile=deconvolution.FragmentFile,ID=ID}
  call report{input:ID=ID,fastq1=fq1,fastq2=fq2,root=root,refdir=refdir,outdir=outdir,Qcfie=qc.Qcfie,Qcfie2=peakcount.Qcfie2}
  call clean{input:outdir=Outdir,html=report.html}
  #call clean{input:outdir=Outdir,html=deconvolution.bapBam}

}

###### script  ######

task makedir {
  String Dir
  command {
              mkdir -p ${Dir}
        mkdir -p ${Dir}/d2cfile
        mkdir -p ${Dir}/out
        mkdir -p ${Dir}/out/Peak
        mkdir -p ${Dir}/out/Promoter
        mkdir -p ${Dir}/temp
        mkdir -p ${Dir}/report/div
              mkdir -p ${Dir}/report/base64
              mkdir -p ${Dir}/report/table
  }
  output {
    String Outdir="${Dir}"
  }
}

task mapping {
      Array[String] fastq1
    Array[String] fastq2
  String outdir
  String refdir
  String ID
  String root
  String refindex
  String read_format

  Int cpp=40
  Int mem=2
  Array[String] barcode

  String  whitelist
  command <<<
     ${root}/software/chromap/v0.2.3_r407/chromap --SAM --preset scATAC-seq  -x ${refindex} -r ${refdir} -1 ${sep=',' fastq1}  -2 ${sep=',' fastq2}  -o ${outdir}/temp/aln.sam --barcode ${sep=',' barcode} --barcode-whitelist ${whitelist} --read-format ${read_format} -t ${cpp}  2> ${outdir}/temp/alignment_report.tsv
    ${root}/software/chromap/v0.2.3_r407/chromap --preset scATAC-seq  --trim-adapters -x ${refindex} -r ${refdir} -1 ${sep=',' fastq1}  -2 ${sep=',' fastq2}  -o ${outdir}/temp/aln.bed --barcode ${sep=',' barcode} --barcode-whitelist ${whitelist} --read-format ${read_format} -t ${cpp}  2> ${outdir}/temp/alignment_report.tsv
   >>>
  runtime{
    backend:"Local"
    cpu:cpp
    memory:"${mem} GB"
  }
 output {
    String  alnbed="${outdir}/temp/aln.bed"
    String aln_report="${outdir}/temp/alignment_report.tsv"
  }
}






task deconvolution {
  String alnbed
  String outdir
  String root
  String ?lib
  String ?species
  String chrmt
  String ID
  String d2cpara
  Int ?ForceFrag
  Int ?mapq
  Int cpp=1
  Int mem=7
  command <<<
   grep -v "chrM" ${alnbed} > ${outdir}/temp/aln.rm_MT.bed 
   ${root}/software/d2c_v1.4.4/bin/d2c merge -i ${outdir}/temp/aln.rm_MT.bed  --mapq 30 --bf 0   -o ${outdir}/d2cfile -c 10 -n ${ID} --mc ${chrmt} --sat --bt1 CB ${d2cpara}
   cp ${outdir}/d2cfile/${ID}.CorrelationBarcodes.tsv.gz ${outdir}/report/plot_input2_Jaccard_Overlap_Knee.csv.gz
  cp ${outdir}/d2cfile/${ID}.barcodeCount.tsv ${outdir}/report/plot_input1_Bead_Barcode_Knee.csv
       tail -n +2 ${outdir}/d2cfile/${ID}.Metadata.tsv |awk '{print $1}' > ${outdir}/d2cfile/list.${ID}.txt

  >>>
  runtime{
    backend:"Local"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
                    String FragmentFile="${outdir}/d2cfile/${ID}.fragments.tsv.gz"
                    String listtxt= "${outdir}/d2cfile/list.${ID}.txt"
                    File FragmentFile1="${outdir}/d2cfile/${ID}.fragments.tsv.gz"
                    File listtxt1= "${outdir}/d2cfile/list.${ID}.txt"

  }
}



task qc {
        String outdir
        String Rscript
        String root
        String tss
        String FragmentFile
        String ID
        String ?lib
        Int cpp=1
        Int mem=1

  command <<<
    if [ -f ${default=abjdbashj lib} ]; then
           source ${lib}
    fi
        ${Rscript} ${root}/script/plot_TSSEnrichment_FragSize.R -T ${tss} -F ${FragmentFile} -G ${ID} -O ${outdir}
        nf=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2<147) print $0}' | wc -l)
        mn=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2>147 && $3-$2<294) print $0}' | wc -l)
        total=$(gunzip -c ${outdir}/d2cfile/${ID}.fragments.tsv.gz | wc -l)
        echo "$nf $total" | awk '{print "Fraction of nucleosome-free regions:"$1/$2*100"%"}' >> ${outdir}/report/5_2.library.QC.csv
        echo "$mn $total" | awk '{print "Fraction of fragments mono-nucleosome regions:"$1/$2*100"%"}' >> ${outdir}/report/5_2.library.QC.csv


 >>>
  runtime{
    backend:"Local"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String Qcfie = "${outdir}/report/5_2.library.QC.csv"
    File Qcfie1 = "${outdir}/report/5_2.library.QC.csv"

  }
}



task peakcount {
  	String root
        String outdir
        String ?lib
        String macs2sp
        String promo
        String ID
	String chrmt
        String FragmentFile
        String listtxt
        Int cpp=1
        Int mem=14

 command <<<
          if [ -f ${default=abjdbashj lib} ]; then
      source ${lib}
          fi
    ${root}/bin/macs2 callpeak -t ${FragmentFile} -f BED -g ${macs2sp} -n ${ID} -B -q 0.001 --nomodel --outdir ${outdir}/out
    ${root}/bin/Rscript ${root}/script/C4scATAC_Cluster_AnnotationV3.R -I ${outdir}/out/${ID}_peaks.narrowPeak -F ${FragmentFile} -C ${listtxt} -G ${promo} -Q ${outdir}/d2cfile/${ID}.Metadata.tsv -O  ${outdir} -MT ${chrmt}

  >>>
  runtime{
    backend:"Local"
    cpu:cpp
    memory:"${mem} GB"
  }
  output {
    String Qcfie2 = "${outdir}/out/${ID}_peaks.narrowPeak"
    File Qcfie21 = "${outdir}/out/${ID}_peaks.narrowPeak"
  }
}



task report {
  String ID
  String ?lib
  Array[String] fastq1
  Array[String] fastq2
  String outdir
            String Qcfie
  String Qcfie2
  String refdir
  String root
  Int cpp=1
  Int mem=5
  command <<<
    if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
    fi

    export LC_ALL=en_US.UTF-8
        echo "Sample ID:${ID}" > ${outdir}/report/2.sample.csv
        echo "FASTQ path:${sep=',' fastq1},${sep=',' fastq2}" >>  ${outdir}/report/2.sample.csv
        echo "Pipeline version:v3.0.0" >> ${outdir}/report/2.sample.csv
        echo "Reference path:${refdir}" >> ${outdir}/report/2.sample.csv
        grep "Number of reads:" ${outdir}/temp/alignment_report.tsv| sed 's/\.//g' | awk -F ":" '{printf "Total number of reads :%'"'"'0.0f\n",$2/2}' > ${outdir}/report/3.sequencing.csv
        grep "Number of barcodes in whitelist\|Number of reads:" ${outdir}/temp/alignment_report.tsv |sed 's/\.//g'|awk -F ':' '{print $2}'|awk 'NR==1{tmp=$1}NR>1{printf "Reads pairs with a valid barcode:%'"'"'0.2f%\n",$1*2/tmp*100}'  >> ${outdir}/report/3.sequencing.csv
        grep "Number of mapped reads:" ${outdir}/temp/alignment_report.tsv | sed 's/\.//g' |awk -F ":" '{print "Reads Mapped to Genome:"$2}' > ${outdir}/report/3.mapping.csv
        grep "Number of uniquely mapped reads:" ${outdir}/temp/alignment_report.tsv | sed 's/\.//g'| awk -F ":" '{printf "Uniquely mapped reads:%'"'"'0.0f\n",$2}' >> ${outdir}/report/3.mapping.csv
       #awk '{print "Mitochondria ratio:"$3}' ${outdir}/temp/aln.bed.mt.stat  >> ${outdir}/report/3.mapping.csv
       echo "Mitochondria ratio:0" >> ${outdir}/report/3.mapping.csv
       grep "bead_cutoff" ${outdir}/d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "Bead threshold:%'"'"'0.0f\n",$2}' > ${outdir}/report/4.cells.csv
        wc -l ${outdir}/d2cfile/${ID}.barcodeMerge.tsv > ${outdir}/report/wc_barcode.txt
        awk -F "," '{printf "Bead number:%'"'"'0.0f\n",$1}' ${outdir}/report/wc_barcode.txt >> ${outdir}/report/4.cells.csv
        grep "cor_cutoff" ${outdir}/d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "cor threshold:%0.5f\n",$2}' >> ${outdir}/report/4.cells.csv
        rm ${outdir}/report/wc_barcode.txt

          ${root}/bin/python ${root}/script/barcode.py --outPath ${outdir}
          ${root}/bin/python ${root}/script/jaccard.py --outPath ${outdir}
    ${root}/bin/python ${root}/script/svg_to_base64_string.py --outPath ${outdir}
          ${root}/bin/python ${root}/script/data-table2.py --outPath ${outdir}
    ${root}/bin/python ${root}/script/st.py --outPath ${outdir} --ID ${ID}
          ${root}/bin/python ${root}/script/generateV2.py --outPath ${outdir} --htmlTemplate ${root}/script/template-test-1-new-2.html --ID ${ID}


  >>>
  runtime{
    backend:"Local"
    cpu:cpp
    memory:"${mem} GB"
  }
  output{
    File html = "${outdir}/report/${ID}_scATAC_analysis_report.html"
    #File html = "${outdir}/report/${ID}_C4-scATAC_report.html"
                                                      }
}



task clean{
  String outdir
  String html
  command{
##    rm -rf ${outdir}/temp/*.bam
##    rm -rf ${outdir}/temp/*.bam.bai
  }
}



 
