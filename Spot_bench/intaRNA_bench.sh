#!/bin/bash

### Run prediction for ryhB sRNA

IntaRNA3 -q query/r.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar.csv


### Run prediction for sgrS sRNA

IntaRNA3 -q query/s.fasta -t target/NC_000913.fa --threads 4 --personality=IntaRNAsTar --outMode=C --out Predictions/res_tar_s.csv
