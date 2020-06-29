#!/usr/bin/env bash

mkdir -p data/reference
mkdir -p data/target
mkdir -p data/genmap

dl_and_preprocess_panels(){
      #Download 1000 Genomes phase 3 data
      wget -O data/reference/chr${1}.vcf.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
      #Initialize target file with modified header (should not contain same samples as reference)
      zcat data/reference/chr${1}.vcf.gz | head -n253  | sed -e 's/HG/TG/g' -e 's/NA/TG/g' | bgzip -c > data/target/chr${1}.vcf.gz
      zcat data/reference/chr${1}.vcf.gz | tail -n+254 > data/tmp_chr${1}.vcf
      #Filter every fifth marker to target data
      cat data/tmp_chr${1}.vcf | awk 'NR%5==1' | bgzip >> data/target/chr${1}.vcf.gz
      #Strip VCF header to file
      HFILE=data/ref.header
      if ! [ -f "$HFILE" ]; then
           zcat data/reference/chr${1}.vcf.gz | head -n253 > data/ref.header
      fi
      #Create tabix index for reference
      tabix -pvcf data/reference/chr${1}.vcf.gz
      rm -f data/tmp_chr${1}.vcf
}

#Run chromosomes in parallel
for chr in {1..22}; do dl_and_preprocess_panels "$chr" & done
#Download and decompress genetic maps
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip -O data/genmap/map.zip
unzip data/genmap/map.zip -d data/genmap/
rm -f data/genmap/map.zip
