# one-liners
My bioinformatics one liners.

## Linux log to file and stdout
`<command> 2>&1 | tee outputfile`

## tsv to csv
`sed 's/\t/,/g' input_file > output_file`

`sed -i 's/\t/,/g' input_file` (inline)

## grep multiple entries
`grep 'word1\|word2' file`

## Check if 2 folders are the same
`diff -r -q /path/dir1 /path/dir2`

## R change ENSGXYZ.X to ENSGXYZ
`sub('\\.[0-9]*$', '', item)`

## fasta files
(https://www.biostars.org/p/17680/)

Add something to the end of headers
`sed 's/>.*/&WHATEVERYOUWANT/' file.fa > outfile.fa`

Extract IDs

`grep -o -E "^>\w+" file.fasta | tr -d ">"`

Linearize sequence and extract IDs

`while read line;do if [ "${line:0:1}" == ">" ]; then echo -e "\n"$line; else echo $line | tr -d '\n' ; fi; done < input.fasta > output.fasta`

`grep -A1 'Q15049' output.fasta`

Alternative:

`awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' input.fasta > output.fasta`

Get the number for reads in fastq

`echo $(( $(wc -l <reads.fq) /4 ))`

## Bed files

Get midpoint windows (get the midpoint of the interval and then add a fixed number of bases)

`awk -vOFS="\t" -vEXT=4 'width=$3-$2 {if(width % 2 != 0) {width+=1} ; mid=$2+width/2; print $1,mid-EXT,mid+EXT,$4}' a.bed`

Calculate on columns

`awk '{sum+=($3-$2)+1} END {print sum}' file.bed`

## Other files

Split a file based on column value (in this case for tab-delimited)

`awk -F '\t' '{print>$13}' myPeaks.bed`
