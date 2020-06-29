#!/usr/bin/env bash

#Remove test data from local FS and HDFS
rm -rf data/reference/
rm -rf data/target/
rm -rf data/genmap/
rm -rf data/ref.header

hdfs dfs -rm -r -f -skipTrash data/reference/  &
hdfs dfs -rm -r -f -skipTrash data/target/ &
hdfs dfs -rm -r -f -skipTrash data/genmap/ &
hdfs dfs -rm -r -f -skipTrash data/ref.header &