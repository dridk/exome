#!/bin/bash 

FOLDER=$1 
SHEET=$2

echo "conversion bcl2fastq"
echo $FOLDER
echo $SHEET 

mkdir -p $FOLDER/output
bcl2fastq --create-fastq-for-index-reads -R $FOLDER -o $FOLDER/output --sample-sheet $SHEET  


mkdir -p $FOLDER/output/merged
cd $FOLDER/output

NAMES=$(ls *.fastq.gz|cut -d"_" -f1|uniq)

echo "fusion des reads par lane"

for i in $NAMES
do
	R1=$i*R1*.fastq.gz
	R2=$i*R2*.fastq.gz
	
	zcat $R1|gzip > merged/$i.R1.fastq.gz &
	zcat $R2|gzip > merged/$i.R2.fastq.gz &
	
done 
