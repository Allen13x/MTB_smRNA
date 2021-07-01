mkdir Predictions

## IntaRNA Predictions
IntaRNA3 -q query/q.fasta -t target/t.fa --threads 4 --outmode=C --out Predictions/res_inta.csv

## sRNARFTarget Predictions

git clone https://github.com/BioinformaticsLabAtMUN/sRNARFTarget.git

cd sRNARFTarget
nextflow run sRNARFTarget.nf --m ../target/t.fa --s ../query/q.fasta

cp sRNARFTargetResult/Prediction_probabilities.csv ../Predictions/Prediction_probabilities.csv
