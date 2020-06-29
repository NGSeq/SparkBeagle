# SparkBeagle
SparkBeagle is distributed implementation of the widely adopted Beagle imputation tool to boost and to scale out genotype imputation from Whole-genome sequencing panels in the cloud environment. 


### Requirements
CentOS Linux 7(recommended), Java 8, Apache Maven 3.x, Hadoop 3.x.x, Spark 2.x.x, HDFS, tabix, bgzip

### Installing
>git clone https://github.com/NGSeq/SparkBeagle

>cd SparkBeagle

>mvn install package

Binary file is created including all needed library dependencies under
target/sparkbeagle-x.x-jar-with-dependencies.jar 

### Cluster setup
Automated installation with Hadoop distributions such as Hortonworks HDP or Cloudera CDH.

or
 
Manually with minimal configuration by following the instructions in https://github.com/NGSeq/SparkBeagle/test/hadoop-init.sh and running the script.

## Running the pipeline as root user
If root user does not have home folder in HDFS:

>su hdfs 

>hdfs dfs -mkdir -p /user/root 

>chown root:root /user/root

>exit

##### Go to test folder
>cd test

### Impute subset of chromosome 22 of 1000 Genomes data set

Data is found under data/chr22_subset folder

##### Put files to HDFS

>hdfs dfs -put data/chr22_subset chr22_subset

##### Submit Spark job

You might need to modify hdfs://localhost:8020 to match the namenode of your cluster. Check also the class path to jar file.
Read https://spark.apache.org/docs/latest/submitting-applications.html if you want to change Spark master or deploy mode.
>spark-submit --master yarn --deploy-mode cluster \
  --executor-memory 4g \
  --driver-memory 4g \
  --executor-cores 8 \
  --conf spark.dynamicAllocation.enabled=true \
  --conf spark.scheduler.mode=FAIR \
  --conf spark.shuffle.service.enabled=true \
  --conf spark.executor.memoryOverhead=1000 \
  --conf spark.port.maxRetries=100 \
  --class org.ngseq.sparkbeagle.Main ../target/sparkbeagle-0.5-jar-with-dependencies.jar \
   -rheader /user/root/chr22_subset/ref.header \
   -ref /user/root/chr22_subset/chr22ref.vcf.gz \
   -gt /user/root/chr22_subset/chr22target.vcf.gz \
   -map /user/root/chr22_subset/plink.chr22.GRCh37.map \
   -index /user/root/chr22_subset/chr22ref.vcf.gz.tbi \
   -out /user/root/chr22_subset/imputed \
   -nthreads 4 \
   -chrom 22 \
   -hdfs hdfs://localhost:8020

##### Get imputed data to local file
>hdfs dfs -getmerge /user/root/chr22_subset/imputed imputed_chr22_subset.vcf.gz

### Impute all autosomes of whole 1000 Genomes data set

##### Run preprocess.sh first to download 1000 Genomes reference panels and create target data. Downloads also genetic maps.
>./preprocess.sh 

##### Execute save2hdfs.sh to load data to HDFS
>./save2hdfs.sh 

##### Run the imputation
>./run_test.sh

##### After imputation, get imputed data to local filesystem
> hdfs dfs -getmerge /user/root/imputed imputed.vcf.gz

or whole directory as is

> hdfs dfs -get /user/root/imputed

### Running imputation in general with multiple chromosomes

##### Put your reference, target and genetic map data to separate folders

>mkdir reference target genmap 

Download genetic map files from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ and unzip to folder named "genmap".

##### Compress VCF data with bgzip and create tabix index (eg. chromosome 1 to 22)

>for i in {1..22};do bgzip reference/chr${i}.vcf &; done

>for i in {1..22};do bgzip target/chr${i}.vcf &; done

>for i in {1..22};do tabix -pvcf reference/chr${i}.vcf.gz &; done

##### Put files to separate HDFS folders

>hdfs dfs -put reference

>hdfs dfs -put target

>hdfs dfs -put genmap

##### Crop and save header text from VCF file to separate ref.header file and put to HDFS

>zcat reference/chr1.vcf.gz | head -n253 > ref.header

>hdfs dfs -put ref.header

##### Modify parameters in run.sh script

##### Execute

>./run.sh

##### Get imputed data to local file

>hdfs dfs -getmerge imputed imputed.vcf.gz

