
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
srnar <- read.delim('Predictions/Prediction_probabilities.csv',header=T,
                    stringsAsFactors = F,sep='\t') %>% 
  rename(id1=mRNA_ID,id2=sRNA_ID,E=Prediction_Probability) %>% 
  mutate(E=1-E,alg='srnar')

inta <-  aggregate(E~id1+id2,
                   read.delim('Predictions/res_inta.csv',header = T,sep = ';',
                              stringsAsFactors = F) %>% 
                     select(id1,id2,E),min) %>% mutate(alg='inta')

duplex<- aggregate(E~id1+id2,
          read.delim('Predictions/res_duplex.csv',header = T,sep = ';',
                     stringsAsFactors = F) %>% 
            select(id1,id2,E),min) %>% mutate(alg='duplex')

tar<- aggregate(E~id1+id2,
                read.delim('Predictions/res_tar.csv',header = T,sep = ';',
                   stringsAsFactors = F) %>% 
                  select(id1,id2,E),min) %>% mutate(alg='sTar')

mir<- aggregate(E~id1+id2,read.delim('Predictions/miranda.csv',header=T,
                   sep=';',stringsAsFactors = F) %>% 
  rename(id2=Seq1,id1=Seq2,E=Max.Energy) %>% 
  select(id1,id2,E),min) %>% mutate(alg='mir')

plex <- aggregate(E~id1+id2,
                  read.delim('Predictions/plex1.csv',header = F,
                   col.names =c('id1','id2','range1','range2','E'),
                   sep=';',stringsAsFactors = F) %>% 
                    select(id1,id2,E) %>% 
                    mutate(E=as.numeric(E)) %>% filter(id1!='null'),min) %>% 
  mutate(alg='plex')

######## Validated Target

ver_r <- c('acnB','erpA','fepB','fur','sdhC',
           'sdhD','marA','fumA','sodB','msrB',
           'shiA','cirA','iscS','nirB','hdeA',
           'cysE','zapB','frdA')

# ver_r<- c('b0118','b0156','b0592','b0683','b0721',
#           'b0722','b1531','b1612','b1656','b1778',
#           'b1981','b2155','b2530','b3365','b3510',
#           'b3607','b3928','b4154')

######## Validated non-Target
fp_r <- c('bfr','sufB','fhuF','sufA','fhuA',
          'fthA','ygdQ')

# fp_r <- c('b3336','b1683','b4367','b1684',
#           'b0150','b1905','b2832')


######## Validated Target

ver_s <- c('ptsG','purR','manX','manY',
           'folE','asd','yigL','adiY')

# ver_s <- c('b1101','b1658','b1817','b1818',
#          'b2153','b3433','b3826','b4116')

######## Validated non-Target

fp_s <- c('zinT','cusF','ykgM','ykgO','znuA',
          'yebA','purR','yeeD','znuC',
          'ydjN','malK','mglB')
# 
# fp_s <- c('b1973','b0573','b0296','b4506','b1857','b1856',
#         'b1178','b2012','b1858','b1729','b4035','b2150')



######## Join predictions dataset
db = rbind(tar,inta,srnar,plex,duplex,mir)


#db = rbind(tar_r,tar_s)

######## Verified dataset target
verified <- data.frame(v=c(rep('t',18),rep('f',7),rep('t',8),rep('f',12)),
                       srna=c(rep('ryhB',25),rep('sgrS',20)),
                       target=c(ver_r,fp_r,ver_s,fp_s))

# verified_1 <- verified


#### Check TP and FP rates - single sRNAs#####
alg <- db %>% select(alg) %>% distinct() %>% pull()
query<- db %>% select(id2) %>% distinct() %>% pull()

d<- data.frame(t=c('a'),srn=c('a'),rel=c('1'),fp=c('1'),a=c('a'))

for (a in alg){
  for (i in query) {
    n<- length(verified %>% filter(srna==i,v=='t') %>% pull(srna))
    f<- length(verified %>% filter(srna==i,v=='f') %>% pull(srna))
for (t in seq(0,2300,100)){
 
  m<-db %>% filter(id2==i,alg==a) %>%
    mutate(r=rank(E))%>% 
    filter(id1%in%(verified %>% filter(srna==i,v=='t') %>% 
                     pull(target)),r<t) %>% 
    summarise(n=n()) %>% pull(n)
  
  z<-db %>% filter(id2==i,alg==a) %>%
    mutate(r=rank(E))%>% 
    filter(id1%in%(verified %>% filter(srna==i,v=='f') %>% 
                     pull(target)),r<t) %>% 
    summarise(n=n()) %>% pull(n)
  
  d<-rbind(d,c(t,i,m/n,z/f,a))
  
}
}
}


d<- d[-1,]


##### ROC figures - single sRNAs

d %>%  mutate_at(.vars = vars(rel,fp,t),as.numeric) %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',
                       low = 'red',high = 'blue' )+
  ggtitle('Roc curve - IntaRNA_sTar single sRNA')+
  geom_abline(slope = 1,intercept = 0)+
  facet_grid(a~srn)


ggsave('Fig/roc_Sgrs_facet.png', dpi=300)



#### Check TP and FP rates - overall sRNAs#####



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

##### ROC figures - single sRNAs


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



ggsave('Fig/roc_Sgrs.png', dpi=300)






# Intersection method -----------------------------------------------------

# # Double check - rank by sTar
# db_c<-inner_join(tar,plex,by=c('id1','id2')) %>% 
#   mutate(alg='mix_sTar',E=E.x) %>% select(-E.x,-E.y,-alg.x,-alg.y)
# # Double check - rank by plex
# db_p<-inner_join(tar,plex,by=c('id1','id2')) %>% 
#   mutate(alg='mix_plex',E=E.y) %>% select(-E.x,-E.y,-alg.x,-alg.y)

######## Create mixed models
db = data.frame()
for (i in list(tar, inta, srnar)) {
  for (j in list(tar,inta,srnar)) {

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





db_int<- rbind(tar, inta) %>% group_by(alg) %>% 
  mutate(E=rank(E)) %>% ungroup() %>% 
  add_count(id1,id2)


db_b<-aggregate(E~id1+id2, db_int%>% filter(n>1),median) %>% tibble() %>% mutate(alg='mix_2')

db_t<-aggregate(E~id1+id2,db_int %>% filter(n>2),median) %>% tibble() %>% mutate(alg='mix_3')
# Aggregate algorithms predictions
db<-rbind(inta,tar,plex,duplex,mir,srnar)
  
aggregate(E~id1+id2+alg,db,min) %>% tibble()->db




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
##### ROC figures - single sRNAs

d  %>% 
  separate(a,c('alg1','alg2'),sep='-') %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'orange',high = 'green' )+
  ggtitle('Roc curve - E.coli - Mixed')+
  geom_abline(slope = 1,intercept = 0)+
  facet_grid(alg1~alg2)


ggsave('Fig/roc_mixalg.png', dpi=300)


# Ranking algorithms
auc=c()

for (i in d %>% select(a) %>% distinct() %>% pull(a)){

auc[i] = with(d %>% filter(a==i)%>%mutate_at(vars(fp,rel),as.numeric),simple_auc(rel,fp))

}

auc[order(-auc)]



