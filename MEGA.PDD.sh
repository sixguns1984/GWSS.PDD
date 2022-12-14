################################
# Pipeline for MEGA PD clincial progression analysis
# Authors: Bioinformatics Team @ Scherzer's lab 
# Email: gliu@mgh.harvard.edu
# Version:0.1
################################
#!/bin/sh

USAGEMSG="usage: $(basename $0) [-i input file Subject][-l input file LONGITUDINAL][-g genotype][-s study sample][-o output file][-f flag]

################################
# Pipeline for genome-wide genotyping imputation Chr1-22
# Authors: Bioinformatics Team @ Scherzer's lab 
# Email: gliu@mgh.harvard.edu
# Version:0.1
################################

*   -i input file with Clinical Individual subjects
*   -l input file with Clinical Longitudinal record 
*	-g input file genotype from PLINK/IMPUTE2
*	-s study file with sample
*	-o output file name	
*   -f flag

Pre-requirement: R/3.1.0 "lmer", "survival" package
"

[ $# -lt 5 ] && echo "$USAGEMSG" >&2 && exit 1

MAF="0.05"

while getopts "i:l:g:s:o:f:" opt
do
	case $opt in
    (i) SUBINPUT=$OPTARG;;
    (l) LONGINPUT=$OPTARG;;
	(g) GENOTYPEINPUT=$OPTARG;;
	(s) MEGASAMPLE=$OPTARG;;
	(o) OUTPUT=$OPTARG;;
	(f) FLAG=$OPTARG;;
	(*) echo "$0: error - unrecognized option $1" >&2; exit 1;;
	esac
done

if [ -f $SUBINPUT ]
then
echo "Find the input Clinical Subject file"
else
echo "no such input Clinical Subject file";exit 0;
fi

if [ -f $LONGINPUT ]
then
echo "Find the input Clinical Longitudinal file"
else
echo "no such input Clinical Longitudinal file";exit 0;
fi

if [ -f $GENOTYPEINPUT ]
then
echo "Find the input study genotype file"
else
echo "no such input study genotype file";exit 0;
fi

if [ -f $MEGASAMPLE ]
then
echo "Find the input study sample file"
else
echo "no such input study sample file";exit 0;
fi

mkdir $FLAG.split
cd    $FLAG.split

# split genotype file into small parts
split -l 20000 $GENOTYPEINPUT

cd ..

for i in $FLAG.split/*; do Rscript MEGA.PDD.PCA.R $SUBINPUT $LONGINPUT $MEGASAMPLE $i $OUTPUT; done

rm -r $FLAG.split
