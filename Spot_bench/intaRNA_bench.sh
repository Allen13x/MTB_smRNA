#!/bin/bash

### Run IntaRNA sTar prediction for ryhB sRNA

IntaRNA3 -q query/r.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar_r.csv


### Run IntaRNA sTar prediction for sgrS sRNA

IntaRNA3 -q query/s.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar_s.csv


### Run miranda prediction for ryhB sRNA

miranda query/r.fasta target/NC_000913.fa -en -4 -out Predictions/miranda_r

{ grep "Seq1,Seq2" Predictions/miranda_r | head -n 1; grep ">>" Predictions/miranda_r; } >> miranda_r

### Run miranda prediction for sgrS sRNA

miranda query/r.fasta target/NC_000913.fa -en -4 -out Predictions/miranda_s

{ grep "Seq1,Seq2" Predictions/miranda_s | head -n 1; grep ">>" Predictions/miranda_s; } >> miranda_s
