#!/bin/bash

#****************************************************
# SYSTEM VERSION
#
# Assembly of novel transcripts in StringTie 
# 
#****************************************************

REFGTF="/Volumes/HD3/Referencias/gencode.v30.annotation.gtf" #Reference fasta file
STRINGTIE="/Users/rafaelcoan/Programas/stringtie-1.3.4d.OSX_x86_64/"
ALIGN_FOLDER="/Volumes/HD3/Brohl_chr/Data/alignments/"

N=1
MAXJOBS=10

function multiple_threads {
    if [ "$N" -lt "$MAXJOBS" ] ; then
        N=$((N + 1))
	else
        wait
    fi
}

echo "*********************************"
echo "STARTING ASSEMBLY WITH STRINGTIE"
echo "*********************************"

cd $ALIGN_FOLDER

for FILE in *.bam
do
    date
    echo $FILE
    FILE_OUT=`echo -e $FILE | sed 's/bam/gtf/'`
    FILE_LABEL=`echo -e $FILE | sed 's/_trim.st.bam/_STRG/'`
    ${STRINGTIE}stringtie -p 20 -G $REFGTF -o $FILE_OUT -l $FILE_LABEL -A ${FILE_OUT}.abund.tab -C ${FILE_OUT}.cov.gtf $FILE
    wait
done

mkdir assembly
mv *.gtf assembly
mv *.tab assembly
mv assembly ..

date
echo "Done."