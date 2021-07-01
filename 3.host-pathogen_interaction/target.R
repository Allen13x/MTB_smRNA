
# R <- Version 4.0.5

# Settings ----------------------------------------------------------------


#Import Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

##### Import Data #####


######## Predictions
srnar <- read.delim('Predictions/Prediction_probabilities.csv',header=T,
                    stringsAsFactors = F,sep='\t') %>% 
  rename(id1=mRNA_ID,id2=sRNA_ID,E=Prediction_Probability) %>% 
  mutate(E=1-E,alg='srnar')

inta <-  aggregate(E~id1+id2,
                   read.delim('Predictions/res_inta.csv',header = T,sep = ';',
                              stringsAsFactors = F) %>% 
                     select(id1,id2,E),min) %>% mutate(alg='inta')




db = data.frame()
for (i in list(srnar)) {
  for (j in list(inta)) {
    
    db_int <- rbind(i %>% mutate(E=rank(E)),
                    j %>% mutate(E=rank(E))) %>%mutate(alg1=i %>% select(alg) %>% 
                                                         distinct() %>% pull(),
                                                       alg2=j %>% select(alg) %>% 
                                                         distinct() %>% pull()) %>% 
      add_count(id1,id2) %>% unite(col = 'alg',alg1,alg2,sep='-')
    
    
    aggregate(E~id1+id2+alg,
              db_int %>% filter(n>1),median) -> db_int
    
    db=rbind(db,db_int)
  }
  
}





db %>% group_by(id2) %>% 
  filter(rank(E)<= 100) %>% ungroup() %>% 
  select(id1,id2) %>% rename(srna=id2,mrna=id1) %>% 
  mutate(tm=paste('/beegfs/scratch/ric.cirillo/dimarco.federico/Human/CNN/couples/',srna,'/mrna/',mrna,'.fasta',sep=''),
         ps=paste('/beegfs/scratch/ric.cirillo/dimarco.federico/Human/CNN/couples/',srna,'/srna/',srna,'.fasta',sep='')) %>% 
  gather(m,pm,tm) %>% select(-m) %>%  
  distinct() %>% arrange(srna,mrna) %>% head()
write_delim('path.csv',col_names = T,delim=',')


