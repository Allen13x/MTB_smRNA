#!/bin/bash

mkdir Predictions
mkdir Raw

### Run IntaRNA sTar prediction for query sRNA

IntaRNA3 -q query/q.fasta -t target/t.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar.csv


## Run IntaRNA3 prediction for query sRNA

IntaRNA3 -q query/q.fasta -t target/t.fa --threads 4 --outmode=C --out Predictions/res_inta.csv

## Run IntaRNA duplex prediction for query sRNA

IntaRNA3 -q query/q.fasta -t target/t.fa --personality=IntaRNAduplex --threads 4 --outmode=C --out Predictions/res_duplex.csv

### Run miranda prediction for query sRNA

miranda query/q.fasta target/t.fa -en -4 -out Raw/miranda

{ grep "Seq1,Seq2" Raw/miranda | head -n 1 | sed 's/,/;/g'; grep ">>" Raw/miranda | sed 's/\t/;/g'| sed 's/ /_/g'| sed 's/>//g'; } > Predictions/res_miranda.csv

### Run RNAplex prediction for query sRNA

RNAplex -q query/q.fasta -t target/t.fa > Raw/plex

i=1
while (($i < "$(wc -l Raw/plex | awk '{print $1}')")) 
do
	j=$i+2
	string=$(cat Raw/plex | awk "NR== $j")
 	if [[ $string == *"&"* ]]
		then
			cat Raw/plex | awk "NR >=$i && NR <= $j" | awk '{print}' ORS=' ' | sed 's/  */ /g' | sed 's/ /\t/g' | awk '{print $1, $2, $4, $6, $7}' | sed 's/ /;/g' | sed 's/[>()]//g'>> Predictions/res_plex.csv 
			i=$i+3
		else
			i=$i+1
	fi	
done

### Run sRNARFTarget prediction for query sRNA
git clone https://github.com/BioinformaticsLabAtMUN/sRNARFTarget.git

cd sRNARFTarget
nextflow run sRNARFTarget.nf --m ../target/t.fa --s ../query/q.fasta

cp sRNARFTargetResult/Prediction_probabilities.csv ../Predictions/Prediction_probabilities.csv
