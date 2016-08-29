#!/bin/bash
echo \# Download choosen fastq files from 1000genome project
echo \# run with: swarm -f /data/wangj36/1000genomes_phase3_exome/download_data.swarm -g 4 -t 8 
filename="Filterd_with_POP.index"
while read  -r line || [ -n "$line" ]; do
	#root=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
	root=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/
	ary=($line)
	if [ ${ary[7]} = EUR ] && [ ${ary[6]} != FIN ] && [ ${ary[6]} = TSI ]
	then
		#echo cmd: wget $root$ary -P /data/wangj36/1000genomes_phase3_exome
		#echo ${ary[0]} ${ary[7]} ${ary[6]}
		if [ ! -d /data/wangj36/1000genomes_phase3_exome/${ary[6]} ]; then
			mkdir /data/wangj36/1000genomes_phase3_exome/${ary[6]}
		fi
		echo wget -c $root$ary -P /data/wangj36/1000genomes_phase3_exome/${ary[6]}
		#if [ ! -f $aim_file]; then
		#	echo wget -c $root$ary -P /data/wangj36/1000genomes_phase3_exome/${ary[6]}
		#fi
	fi

done<"$filename"
