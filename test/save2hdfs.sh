#!/usr/bin/env bash

# Create home folder for HDFS root user
#su hdfs
#hdfs dfs -mkdir -p /user/root
#hdfs dfs -chown root:root /user/root
#exit

#change hostname to match your namenode
HDFS_HOST=localhost

hdfs dfs -put data/reference/ hdfs://${HDFS_HOST}:8020/user/root/ &
hdfs dfs -put data/target/ hdfs://${HDFS_HOST}:8020/user/root/ &
hdfs dfs -put data/genmap/ hdfs://${HDFS_HOST}:8020/user/root/ &
hdfs dfs -put data/ref.header hdfs://${HDFS_HOST}:8020/user/root/ &
