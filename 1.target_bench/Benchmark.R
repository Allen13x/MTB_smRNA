
# R <- Version 4.0.2

# Settings ----------------------------------------------------------------


#Import Libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

#Custom functions
simple_auc <- function(rel, fp){
  dFPR <- c(diff(fp), 0)
  dTPR <- c(diff(rel), 0)
  sum(rel * dFPR) + sum(dTPR * dFPR)/2
}


##### Import Data #####


######## Predictions


###sRNARFTarget
srnar <- read.delim('Predictions/Prediction_probabilities.csv',header=T,
                    stringsAsFactors = F,sep='\t') %>% 
  rename(id1=mRNA_ID,id2=sRNA_ID,E=Prediction_Probability) %>% 
  mutate(E=1-E,alg='srnar')
###IntaRNA3
inta <-  aggregate(E~id1+id2,
                   read.delim('Predictions/res_inta.csv',header = T,sep = ';',
                              stringsAsFactors = F) %>% 
                     select(id1,id2,E),min) %>% mutate(alg='inta')
###IntaRNADuplex
duplex<- aggregate(E~id1+id2,
          read.delim('Predictions/res_duplex.csv',header = T,sep = ';',
                     stringsAsFactors = F) %>% 
            select(id1,id2,E),min) %>% mutate(alg='duplex')

###IntaRNAsTar
tar<- aggregate(E~id1+id2,
                read.delim('Predictions/res_tar.csv',header = T,sep = ';',
                   stringsAsFactors = F) %>% 
                  select(id1,id2,E),min) %>% mutate(alg='sTar')

###MiRanda
mir<- aggregate(E~id1+id2,read.delim('Predictions/miranda.csv',header=T,
                   sep=';',stringsAsFactors = F) %>% 
  rename(id2=Seq1,id1=Seq2,E=Max.Energy) %>% 
  select(id1,id2,E),min) %>% mutate(alg='mir')


###RNAPlex
plex <- aggregate(E~id1+id2,
                  read.delim('Predictions/res_plex.csv',header = F,
                   col.names =c('id1','id2','range1','range2','E'),
                   sep=';',stringsAsFactors = F) %>% 
                    select(id1,id2,E) %>% 
                    mutate(E=as.numeric(E)) %>% filter(id1!='null'),min) %>% 
  mutate(alg='plex')

######## Validated Target

ver_1 <- c()  ## Insert list validated target srna-1


######## Validated non-Target
fp_1 <- c()   ## Insert list validated non-target srna-1



######## Validated Target



ver_2 <- c()  ## Insert list validated target srna-2


######## Validated non-Target
fp_2 <- c()   ## Insert list validated non-target srna-2


######## Join predictions dataset
db = rbind(tar,inta,srnar,plex,duplex,mir)




######## Verified dataset target (replace srna-1/2 with smRNA-id )
verified <- data.frame(v=c(rep('t',length(ver_1)),rep('f',length(fp_1)),rep('t',length(ver_2)),rep('f',length(fp_2))),
                       srna=c(rep('srna-1',length(ver_1)+length(fp_1)),rep('srna-2',length(ver_2)+length(fp_2))),
                       target=c(ver_1,fp_1,ver_2,fp_2))


#### Check TP and FP rates - overall sRNAs#####
alg <- db %>% select(alg) %>% distinct() %>% pull()
query<- db %>% select(id2) %>% distinct() %>% pull()


d<- data.frame(t=c('a'),srn=c('a'),rel=c('1'),fp=c('1'),a=c('a'))

for (a in alg){
for (t in seq(0,9000,100)){

    vt<-verified %>% filter(v=='t')
    vf<-verified %>% filter(v=='f')
    
    n<- length(vt %>% pull(srna))
    f<- length(vf %>% pull(srna))
    m<-db %>% filter(alg==a) %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vt,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    z<-db %>% filter(alg==a) %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vf,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    d<-rbind(d,c(t,'srna',m/n,z/f,a))
    

}
  d<- rbind(d,c(t+100,'srna',1,1,a))
  }





d=d[-1,]


# Ranking algorithms
auc=c()

for (i in d %>% select(a) %>% distinct() %>% pull(a)){
  
  auc[i] = with(d %>% filter(a==i)%>%mutate_at(vars(fp,rel),as.numeric),simple_auc(rel,fp))
  
}

auc[order(-auc)]

##### ROC figures


d  %>%  mutate_at(.vars = vars(rel,fp,t),as.numeric) %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'orange',high = 'red' )+
  ggtitle('Algorithm ROC Curve')+
  geom_abline(slope = 1,intercept = 0)+
  geom_label(data = data.frame(auc) %>% rownames_to_column('a'),
              aes(label=paste('AUC=',round(auc,2)),y=0.10,x=0.75))+
  facet_wrap(~a,labeller = labeller(a=c('duplex' = 'IntaRNA-Duplex','inta'='IntaRNA3','mir'='miRanda','plex'='RNAPlex',
                                        'srnar'='RNARFTarget','sTar'='IntaRna-sTar')))



ggsave('Fig/roc_joint.png', dpi=300)






#### Joint method -----------------------------------------------------


### Create joint df

list_db <- list(tar, inta, srnar) # instert prediction for desired algorithms

db = data.frame()
for (i in list_db) {
  for (j in list_db) {

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


#### Check TP and FP rates - overall sRNAs#####

alg <- db %>% select(alg) %>% distinct() %>% pull()
query<- db %>% select(id2) %>% distinct() %>% pull()

d<- data.frame(t=c('a'),srn=c('a'),rel=c('1'),fp=c('1'),a=c('a'))

for (a in alg){
  for (t in seq(0,9000,100)){
    
    vt<-verified %>% filter(v=='t')
    vf<-verified %>% filter(v=='f')
    
    n<- length(vt %>% pull(srna))
    f<- length(vf %>% pull(srna))
    m<-db %>% filter(alg==a) %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vt,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    z<-db %>% filter(alg==a) %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vf,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    d<-rbind(d,c(t,'srna',m/n,z/f,a))
    
    
  }
  d<- rbind(d,c(t+100,'srna',1,1,a))
  }

d=d[-1,]
d  %>%  mutate_at(.vars = vars(rel,fp,t),as.numeric)-> d

# Ranking algorithms
auc=c()

for (i in d %>% select(a) %>% distinct() %>% pull(a)){
  
  auc[i] = with(d %>% filter(a==i)%>%mutate_at(vars(fp,rel),as.numeric),simple_auc(rel,fp))
  
}

auc[order(-auc)]
##### ROC figures - single sRNAs

d  %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'orange',high = 'red' )+
  ggtitle('Roc curve - Joint')+
  geom_abline(slope = 1,intercept = 0)+
  geom_label(data = data.frame(auc) %>% rownames_to_column('a'),
             aes(label=paste('AUC=',round(auc,2)),y=0.10,x=0.75))+
  facet_wrap(~a)


ggsave('Fig/roc_mixalg.png', dpi=300)






