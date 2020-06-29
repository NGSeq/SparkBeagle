#!/usr/bin/env bash

###### USE ONLY FOR TESTING PURPOSES. USE AT YOUR OWN RISK. ######

#Assume we have three nodes having hostnaems node-1,node-2,node-3 and we are logged in as root user
#Hadoop is installed in folder /opt/hadoop and Spark in /opt/spark

#####Get Hadoop and Spark
wget https://archive.apache.org/dist/hadoop/core/hadoop-2.7.3/hadoop-2.7.3.tar.gz -O /opt/hadoop-2.7.3.tar.gz
tar -xf /opt/hadoop-2.7.3.tar.gz
mv /opt/hadoop-2.7.3 /opt/hadoop

wget wget http://archive.apache.org/dist/spark/spark-2.3.2/spark-2.3.2-bin-hadoop2.7.tgz -O /opt/spark-2.3.2-bin-hadoop2.7.tgz
tar -xf /opt/spark-2.3.2-bin-hadoop2.7.tgz
mv /opt/spark-2.3.2-bin-hadoop2.7 /opt/spark

#Modify hosts, HDFS dirs, memory options in /opt/hadoop/etc/hadoop/core-site.xml, /opt/hadoop/etc/hadoop/hdfs-site.xml, /opt/hadoop/etc/hadoop/yarn-site.xml and /opt/spark/conf/spark-defaults.conf
#add all nodes to /etc/hosts and /opt/hadoop/etc/hadoop/workers file

## Minimal configurations

##### /opt/hadoop/etc/hadoop/core-site.xml
#<configuration>
#	<property>
#            <name>fs.default.name</name>
#            <value>hdfs://node-1:8020</value>
#        </property>
#</configuration>

##### /opt/hadoop/etc/hadoop/hdfs-site.xml
#<configuration>
#    <property>
#            <name>dfs.namenode.name.dir</name>
#            <value>/hadoop/dfs/nn</value>
#    </property>
#
#    <property>
#            <name>dfs.datanode.data.dir</name>
#            <value>/hadoop/dfs/dn</value>
#    </property>
#
#    <property>
#            <name>dfs.replication</name>
#            <value>3</value>
#    </property>
#</configuration>

##### /opt/hadoop/etc/hadoop/yarn-site.xml
#
#<configuration>
#
#<!-- Site specific YARN configuration properties -->
#    <property>
#            <name>yarn.acl.enable</name>
#            <value>0</value>
#    </property>
#
#    <property>
#            <name>yarn.resourcemanager.hostname</name>
#            <value>node-1</value>
#    </property>
#
#    <property>
#            <name>yarn.nodemanager.aux-services</name>
#            <value>mapreduce_shuffle</value>
#    </property>
#    <property>
#            <name>yarn.nodemanager.resource.memory-mb</name>
#            <value>5000</value>
#    </property>
#
#    <property>
#            <name>yarn.scheduler.maximum-allocation-mb</name>
#            <value>30000</value>
#    </property>
#
#    <property>
#            <name>yarn.scheduler.minimum-allocation-mb</name>
#            <value>1024</value>
#    </property>
#
#    <property>
#            <name>yarn.nodemanager.vmem-check-enabled</name>
#            <value>false</value>
#    </property>
#</configuration>

#### /opt/hadoop/etc/hadoop/workers
#node-1
#node-2
#node-3

#### /opt/spark/conf/spark-defaults.conf
#spark.eventLog.enabled           true
#spark.eventLog.dir               hdfs://localhost:8020/tmp/applicationHistory
#spark.history.fs.logDirectory	  /tmp/applicationHistory

## Environment variables, add to .bashrc
export HDFS_NAMENODE_USER="root"
export HDFS_DATANODE_USER="root"
export HDFS_SECONDARYNAMENODE_USER="root"
export YARN_RESOURCEMANAGER_USER="root"
export YARN_NODEMANAGER_USER="root"

export YARN_CONF_DIR=/opt/hadoop/etc/hadoop/
export PATH=/opt/hadoop/bin:/opt/hadoop/sbin:/opt/spark/bin:$PATH
#export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk #Modify and uncomment if java is not found

#Generate ssh key
ssh-keygen
#copy ssh key to all nodes
ssh-copy-id -i $HOME/.ssh/id_rsa.pub root@node-1
ssh-copy-id -i $HOME/.ssh/id_rsa.pub root@node-2
ssh-copy-id -i $HOME/.ssh/id_rsa.pub root@node-3

distribute(){

    #####Install java-1.8.0-openjdk
    ssh -tt node-$i yum install -y java-1.8.0-openjdk
    #Copy Hadoop and Spark
    ssh -tt node-$i mkdir /hadoop
    scp -r /opt/hadoop node-$i:/
    scp -r /opt/spark node-$i:/
    #distribute hosts file
    scp /etc/hosts node-$i:/etc/hosts

    #Copy only configurations if changed
    #scp /opt/hadoop/etc/hadoop/* node-$i:/opt/hadoop/etc/hadoop/;

    #If more disk space is needed, add volumes to all nodes and mount to /hadoop (default hdfs dir will be /hadoop/dfs assigned in /opt/hadoop/etc/hadoop/hdfs-site.xml)
    #ssh -tt node-$i mkfs.ext4 /dev/vdb
    #ssh -tt node-$i mount /dev/vdb /hadoop

    #Remove all HDFS data if needed
    #ssh -tt node-$i rm -r -f /hadoop/*
}
#change iterator to match hostnames
for i in {1..3}; do distribute "$i" & done

#start services
hdfs namenode -format
/opt/hadoop/sbin/start-all.sh
hdfs dfs -mkdir /tmp
hdfs dfs -mkdir /user
hdfs dfs -mkdir /user/root