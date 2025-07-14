#!/bin/bash

#****************************************************
# NORMAL VERSION
#
# Run this script if a second alignment is necessary
# due to strand information and after SeqMonk
#****************************************************

#samtools 1.7
#Using htslib 1.7

REFFOLDER="/Volumes/HD3/Referencias/" 
ANASILYSFOLDER="/Volumes/HD3/Brohl_chr/Data/"
READSFOLDER="/Volumes/HD3/Brohl_chr/Data/trimmed-reads/"
RESULTSFOLDER="/Volumes/HD3/Brohl_chr/Data/alignments/"

REFHG38="hg38_chromosomes.fa" #Reference fasta file
INDEXLIB="hg38_chromosomes" #Library name for hisat2-build
STRAND= #Library strand (Single-end: F -transcript- or R -reverse complement; Pair-end: FR or RF)
HISAT2="/Users/rafaelcoan/Programas/hisat2-2.1.0/"

N=1
MAXJOBS=10

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

cd $REFFOLDER
ls -1 $(pwd)/*.ht2 > /dev/null 2>&1
if [ "$?" = "0" ]; then
    echo "Index exists!"
else
    docker run -w /mnt -v $(PWD):/mnt --rm kathrinklee/rna-seq-pipeline-hisat2 hisat2-build -p 20 $REFHG38 $INDEXLIB
fi
cd $ANASILYSFOLDER
mkdir alignments
cd $RESULTSFOLDER

echo "******************************"
echo "STARTING ALIGNMENT WITH HISAT2"
echo "******************************"

for FILE_R1 in ${READSFOLDER}*_1.trim.gz
do
    date
    FILE_R2=`echo -e $FILE_R1 | sed 's/_1.trim.gz/_2.trim.gz/'` #Input reverse
    FILE_SAM=`echo -e $FILE_R1 | sed 's/_1.trim.gz/_trim.sam/'`
    FILE_BAM=`echo -e $FILE_SAM | sed 's/sam/st.bam/'`
    echo $FILE_BAM

    ${HISAT2}hisat2 -p 20 --dta -x ${REFFOLDER}${INDEXLIB} -1 $FILE_R1 -2 $FILE_R2 -S ${FILE_SAM} --new-summary --summary-file ${FILE_SAM}.log
    wait
    
    samtools view -@ 15 -bSh -f 3 -q 30 $FILE_SAM | samtools sort -@ 15 - > $FILE_BAM
    wait
    rm $FILE_SAM
    mv $FILE_BAM $RESULTSFOLDER
done

echo "********************"
echo "SAMTOOLS CORRECTIONS"
echo "********************"

for FILE in *.bam
do
    date
    echo $FILE
    FILE_OUT=`echo -e $FILE | sed 's/bam/bam.bai/'`

    samtools index -@ 20 $FILE $FILE_OUT &
    wait
done
