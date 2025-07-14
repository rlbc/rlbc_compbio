#!/bin/bash

#********************************************************
# DOCKER VERSION
#
# This script makes salmon quantification on transcripts.
# Salmon 0.12
#********************************************************

REFFASTA="gencode.v30.pc_lncRNA.fa" #Reference fasta file; previously created a hard link for the file
INDEXLIB="gencode.v30.transcriptsLib"
READSFOLDER="/Volumes/HD3/Brohl_chr/Data/trimmed-reads/"

N=1
MAXJOBS=3

function multiple_threads {
    if [ "$N" -lt "$MAXJOBS" ] ; then
        N=$((N + 1))
	else
        wait
    fi
}

echo "***********************"
echo "STARTING BUILDING INDEX"
echo "***********************"
date

mkdir /Volumes/HD3/Brohl_chr/Data/quantification
cd /Volumes/HD3/Brohl_chr/Data/quantification
ln /Volumes/HD3/Referencias/gencode.v30.pc_lncRNA.fa

ls -1 gencode.v30.transcriptsLib* > /dev/null 2>&1
if [ "$?" = "0" ]; then
    echo "Index exists!"
else
    docker run -w /mnt -v $(PWD):/mnt --rm combinelab/salmon salmon index -t $REFFASTA -i $INDEXLIB
fi

echo "************************************"
echo "STARTING QUANTIFICATION WITH SALMON"
echo "************************************"

for FILE_R1 in $READSFOLDER*_1.trim.gz
do
    FILE_R2=`echo -e $FILE_R1 | sed 's/_1.trim.gz/_2.trim.gz/'` #Input reverse
    ln $FILE_R1 $(PWD)
    ln $FILE_R2 $(PWD)
    FILE_R1=`echo -e $FILE_R1 | cut -d "/" -f 7`
    FILE_R2=`echo -e $FILE_R2 | cut -d "/" -f 7`
    date
    echo -e $FILE_R1
    echo -e $FILE_R2
    FILE_OUT=`echo -e $FILE_R1 | sed 's/_1.trim.gz/_quant/'`
    docker run -w /mnt -v $(PWD):/mnt --rm combinelab/salmon salmon quant -i $INDEXLIB -l A -1 $FILE_R1 -2 $FILE_R2 -p 20 --gcBias --validateMappings -o quants/${FILE_OUT}
    wait
    rm $FILE_R1 $FILE_R2
done
