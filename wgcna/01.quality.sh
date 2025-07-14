#!/bin/bash

#***************
# Quality check
#
# FastQC
#***************

# Directories
DATA=/Volumes/HD2Backup/Coan/dbGaP-17758/sra
VOLUME=/Volumes/HD3/Brohl_chr

cd $VOLUME
mkdir Data/fastqc
cd Data/fastqc

for file in ls ${DATA}/*.fastq.gz
do
  cp $file .
done

N=1
MAXJOBS=8

function multiple_threads {
    if [ "$N" -lt "$MAXJOBS" ] ; then
        N=$((N + 1))
	else
        wait
    fi
}

echo "********************************"
echo "STARTING QC ANALYSIS WITH FASTQC"
echo "********************************"
date

for FILE in *.fastq.gz #Change when necessary
do
  echo -e "Iniciando arquivo $FILE"
  docker run --rm -v $(PWD):/data biocontainers/fastqc:0.11.5 fastqc $FILE &
  multiple_threads
done

rm *.fastq.gz