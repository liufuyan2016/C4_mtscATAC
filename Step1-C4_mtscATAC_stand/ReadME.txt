Payatterntion：requriment
##download from https://github.com/M-wen/C4_scATAC_analysis scripts to replace path in the C4_mtscATAC_nonmodel_local_V3.0_BJ.wdl scripts.

35:  string root="/jdfsbjcas1/ST_BJ/PUB/User/kangjingmin/Pipline/DNBelabC4scATAC/"

##change the data.config input fastq and output directory and could run your data
/usr/bin/java -jar cromwell-35.jar run  C4_mtscATAC_nonmodel_local_V3.0_BJ.wdl  --inputs data.config
