# sRNA targeting host-pathogen

Code Repository for the project of Bioinformatic characterization of Mycobacterium tuberculosis smRNA in Host-Pathogen interaction

No data will be shared, just codes!

The overall workflow is divided in sub-folders ordered by number


0.smRNA_Candidates -> 1.Target_benchamrking -> 2.Matrix_workflow -> 3.Host-Pathogen_Interactions

Each folder containts the req.txt file to create a conda environment with all the package needed to run the scripts (except R)

## 0.smRNA_Candidates

Starting from the fastq files, the pipeline.sh script will produce the file counts/mtb_counts.tsv file with the counts over the candidates smRNA and mRNA from Mycobacterium Tuberculosis (H37Rv)

The general scheme is remove the reads which allign to Human Genome and then use all the remaining reads which map on H37Rv genome to produce the count table

The count.R file will extract the smRNA which pass the count filters (Candidates present in at least 2 samples of the group ATB, or 3 samples between ATB and LTBI and non present in the CTRL group or present in more than 30 samples between ATB and LTBI). A sample table must be provided by the user.


## 1.Target_benchmarking

The script srna_bench.sh will take as input the query small rnas from a fasta file in query/q.fasta and the target fasta from target/t.fa

The predictions from each algorithm will be stored in the Predictions/ folder.

The Benchmark.R script will compare the AUC of each algorithm and joint predictions and produce two figures in Fig/ folder. (True Positive and False Positive targets must be provied in the script)

## 2.Target_scan

### req1.txt

A nextflow pipeline to produce the training set for the Tensorflow model is contained in the script Matrix_maker.nf. The input file must have the format:

smRNA;mRNA;NCBI Reference Sequence;TAXID;commonID;smRNA sequence;mRNA sequence

Use the following code to launch the pipeline
```
nextflow run Matrix_maker.nf --train --list INPUT --outdir Train_dataset
```
Replace INPUT with the path to the input file

### req2.txt

The Nextflow pipeline will produce a Train_dataset folder with the energy variation matrixes inside the two folders t and f (Labels for the model training)


The codes inside tf.py will create, train and store the model based on the content of Train_dataset folder.

## 3.Host-Pathogen_Interactions

### req1.txt

Put the final candidates smRNAs fasta file in the query folder as q.fasta and the human target genome in the the folder target as t.fa.

Run the script Target_prediction.sh. The results will be available inside the folder Predictions.

The target.R script will create the joint prediction dataset and produce the file with the list of top 100 targets for each candidate.

Create a folder in which store the single fasta of each smRNA and mRNA.

Create an input file with this format:

mRNA;i;smRNA;path to mRNA .fasta file; path to smRNA .fasta file

"i" is just a custom label.

and use the header: mrna;i;srna;pm;ps

Use the following code to launch the pipeline
```
nextflow run ../2.matrix_workflow/Matrix_maker.nf --list INPUT
```
This will produce a results/ folder containing the energy variation matrixes for all the couples. 

### req2.txt

The code inside the model_filter.py will use the trained model from point 2 to discriminate true and false positive targets (True:1; False:0).

The code inside TopGO.R will perform the TopGO analysis over the true targets (interactive mode - e.g RStudio)








