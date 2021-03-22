#!/bin/bash

mkdir Predictions
mkdir Raw

### Run IntaRNA sTar prediction for ryhB sRNA

IntaRNA3 -q query/r.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar_r.csv


### Run IntaRNA sTar prediction for sgrS sRNA

IntaRNA3 -q query/s.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar_s.csv


### Run miranda prediction for ryhB sRNA

miranda query/r.fasta target/NC_000913.fa -en -4 -out Raw/miranda_r

{ grep "Seq1,Seq2" Raw/miranda_r | head -n 1 | sed 's/,/;/g'; grep ">>" Raw/miranda_r | sed 's/\t/;/g'| sed 's/ /_/g'; } > Predictions/miranda_r.csv

### Run miranda prediction for sgrS sRNA

miranda query/s.fasta target/NC_000913.fa -en -4 -out Raw/miranda_s

{ grep "Seq1,Seq2" Raw/miranda_s | head -n 1 | sed 's/,/;/g'; grep ">>" Raw/miranda_s | sed 's/\t/;/g'| sed 's/ /_/g'; } > Predictions/miranda_s.csv

