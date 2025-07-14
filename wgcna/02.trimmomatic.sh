#!/bin/bash

#************************************************
# TRIMMOMATIC DOCKER VERSION
#
# http://www.usadellab.org/cms/?page=trimmomatic
#
# Run inside DATA directory
#************************************************

# Directories
DATA=/Volumes/HD2Backup/Coan/dbGaP-17758/sra
VOLUME=/Volumes/HD3/Brohl_chr

N=1
MAXJOBS=8

function multiple_threads {
    if [ "$N" -lt "$MAXJOBS" ] ; then
        N=$((N + 1))
	else
        wait
    fi
}

cd $VOLUME
mkdir Data/trimmed-reads
cd Data/trimmed-reads

for file in ls ${DATA}/*.fastq.gz
do
  cp $file . &
  multiple_threads
done

echo "**********************************"
echo "STARTING TRIMMING WITH TRIMMOMATIC"
echo "**********************************"
date

for FILE_R1 in *_1.fastq.gz
do
    FILE_R2=`echo -e $FILE_R1 | sed 's/_1.fastq.gz/_2.fastq.gz/'` #Input reverse
    FILE_OUT_R1=`echo -e $FILE_R1 | sed 's/_1.fastq.gz/_1.trim.gz/'` #Output forward paired
    FILE_OUT_R2=`echo -e $FILE_R2 | sed 's/_2.fastq.gz/_2.trim.gz/'` #Output reverse paired

    FILE_OUT_R1_UNP=`echo -e $FILE_OUT_R1 | sed 's/_1.trim.gz/_1.trim-unpaired.gz/'` #Output forward unpaired
    FILE_OUT_R2_UNP=`echo -e $FILE_OUT_R2 | sed 's/_2.trim.gz/_2.trim-unpaired.gz/'` #Output reverse unpaired

    docker run -w /home -v $(PWD):/home --rm fjukstad/trimmomatic PE \
    -threads 2 $FILE_R1 $FILE_R2 $FILE_OUT_R1 $FILE_OUT_R1_UNP $FILE_OUT_R2 $FILE_OUT_R2_UNP \
    ILLUMINACLIP:/tools/trimmomatic/adapters/TruSeq3-PE-2.fa:4:30:10 SLIDINGWINDOW:4:20 MINLEN:50 AVGQUAL:28 &
    
    multiple_threads #Coloca 4 trabalhos de uma só vez
done

rm *.fastq.gz
mkdir unpaired-reads
mv *.trim-unpaired.gz unpaired-reads
