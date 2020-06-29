##This folder contains scripts for performing imputation scalability tests on Apache Spark cluster

###Running the test
-Run preprocess.sh first to download 1000 Genomes reference panels and create target data. Downloads also genetic maps.
./preprocess.sh 

-Execute save2hdfs.sh to load data to HDFS
./save2hdfs.sh 

-Run the imputation
./run_test.sh

-After imputation is ready, get the imputed data to local filesystem
> hdfs dfs -getmerge /user/root/imputed imputed.vcf.gz
or whole directory as is
> hdfs dfs -get /user/root/imputed

