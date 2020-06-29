#!/usr/bin/env bash

#change hostname to match your namenode
HDFS_HOST=localhost
#All paths are HDFS paths
REFERENCE_HEADER=/user/root/ref.header #this file must contain the header fields of the reference VCF
REFERENCE_PATH=/user/root/reference
TARGET_PATH=/user/root/target
INDEX_PATH=/user/root/reference
GENMAP_PATH=/user/root/genmap
OUT_PATH=/user/root/imputed

#parameters: chromosome executor_memory(GB) driver_memory(GB) window(interval length cM) overlap(cM) threads executor-cores
run() {
  MODE=cluster #client #Read https://spark.apache.org/docs/latest/submitting-applications.html
  MASTER=yarn #local[8] #spark://localhost:7077
  spark-submit --master yarn --deploy-mode ${MODE} \
  --executor-memory ${2}g \
  --driver-memory ${3}g \
  --executor-cores ${7} \
  --conf spark.dynamicAllocation.enabled=true \
  --conf spark.scheduler.mode=FAIR \
  --conf spark.shuffle.service.enabled=true \
  --conf spark.executor.memoryOverhead=1000 \
  --conf spark.port.maxRetries=100 \
  --class org.ngseq.sparkbeagle.Main ../target/sparkbeagle-0.5-jar-with-dependencies.jar \
   -rheader ${REFERENCE_HEADER} \
   -ref ${REFERENCE_PATH}/chr${1}.vcf.gz \
   -gt ${TARGET_PATH}/chr${1}.vcf.gz \
   -map ${GENMAP_PATH}/plink.chr${1}.GRCh37.map \
   -window ${4} -overlap ${5} \
   -index ${INDEX_PATH}/chr${1}.vcf.gz.tbi \
   -out ${OUT_PATH}/chr${1} \
   -nthreads ${6} \
   -chrom ${1} \
   -hdfs hdfs://${HDFS_HOST}:8020 &>imputation${1}.log
}

for i in {1..22}; do run "$i" 4 4 1 0.1 4 4& done




