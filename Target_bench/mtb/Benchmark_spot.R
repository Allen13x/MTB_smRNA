
# R <- Version 4.0.2

# Settings ----------------------------------------------------------------


#Import Libraries
library(dplyr)
library(tidyr)
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
  mutate(E=1-E,alg='srnar') %>% separate(id1,c('id1'),sep = '\\|')

inta <-  aggregate(E~id1+id2,
                   read.delim('Predictions/res_inta.csv',header = T,sep = ';',
                              stringsAsFactors = F) %>% 
                     select(id1,id2,E),min) %>% mutate(alg='inta') %>% separate(id1,c('id1'),sep = '\\|')

duplex<- aggregate(E~id1+id2,
          read.delim('Predictions/res_duplex.csv',header = T,sep = ';',
                     stringsAsFactors = F) %>% 
            select(id1,id2,E),min) %>% mutate(alg='duplex') %>% separate(id1,c('id1'),sep = '\\|')

tar<- aggregate(E~id1+id2,
                read.delim('Predictions/res_tar.csv',header = T,sep = ';',
                   stringsAsFactors = F) %>% 
                  select(id1,id2,E),min) %>% mutate(alg='sTar') %>% separate(id1,c('id1'),sep = '\\|')

mir<- aggregate(E~id1+id2,read.delim('Predictions/miranda.csv',header=T,
                   sep=';',stringsAsFactors = F) %>% 
  rename(id2=Seq1,id1=Seq2,E=Max.Energy) %>% 
  select(id1,id2,E),min) %>% mutate(alg='mir')

plex <- aggregate(E~id1+id2,
                  read.delim('Predictions/plex.csv',header = F,
                   col.names =c('id1','id2','range1','range2','E'),
                   sep=';',stringsAsFactors = F) %>% 
                    select(id1,id2,E) %>% 
                    mutate(E=as.numeric(E)) %>% filter(id1!='null'),min) %>% 
  mutate(alg='plex') %>% separate(id1,c('id1'),sep = '\\|')

######## Validated Target
v<- c(9,7,4)


ver<-c('Rv0058','Rv3864','Rv3865','Rv3870','Rv3872',
       'Rv3883c','Rv0172','Rv0174','Rv0285', ##<----B11
       'Rv2093c','Rv2618','Rv0092','Rv2221c','Rv0378',
       'Rv2013','Rv2997',#<----Mcr7
       'Rv1876','Rv3106','Rv0009','Rv0069c' #<----MrsI
       ) 

# ver_r <- c('dnaB','espE','espF','eccCa1','PE35',
#            'mycP1','mce1D','mce1F','PE5')


######## Validated non-Target
fp <- c('	Rv1876','Rv3841')

# fp <- c('bfrA','bfrB','sufR')

# fp_r <- c('b3336','b1683','b4367','b1684',
#           'b0150','b1905','b2832')





######## Join predictions dataset
db = rbind(tar_r,tar_s,plex)


#db = rbind(tar_r,tar_s)

######## Verified dataset target
verified <- data.frame(v=c(rep('t',length(ver)),rep('f',3*length(fp))),
                       srna=c(rep('B11',v[1]),rep('Mcr7',v[2]),rep('MrsI',v[3]),rep(v,each=2)),
                       target=c(ver,fp,fp,fp))

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
    

}}





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
  scale_color_gradient(name='Threshold - Top Targets',low = 'red',high = 'blue' )+
  ggtitle('Roc curve - IntaRNA_sTar_joint')+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~a)



ggsave('Fig/roc_Sgrs.png', dpi=300)






# Intersection method -----------------------------------------------------

# # Double check - rank by sTar
# db_c<-inner_join(tar,plex,by=c('id1','id2')) %>% 
#   mutate(alg='mix_sTar',E=E.x) %>% select(-E.x,-E.y,-alg.x,-alg.y)
# # Double check - rank by plex
# db_p<-inner_join(tar,plex,by=c('id1','id2')) %>% 
#   mutate(alg='mix_plex',E=E.y) %>% select(-E.x,-E.y,-alg.x,-alg.y)

# Double check - rank by min

db = data.frame()
for (i in list(tar, srnar,inta)) {
  for (j in list(tar,srnar,inta)) {
    
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

# Aggregate algorithms predictions
db<-rbind(inta,tar,plex,duplex,db_b,srnar)
  
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
d  %>%  mutate_at(.vars = vars(rel,fp,t),as.numeric) -> d
##### ROC figures - single sRNAs

d  %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'orange',high = 'green' )+
  ggtitle('Roc curve - IntaRNA_sTar_joint')+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~a)



d %>% 
  separate(a,c('alg1','alg2'),sep='-') %>% 
  ggplot(aes(y=rel,x=t))+
  geom_path()+
  xlab('Thresholds')+ylab('True Positive Rates')+
  theme_bw()+
  facet_grid(alg1~alg2)
  
  
ggsave('Fig/TP_genome.png', dpi=300)


# Ranking algorithms
auc=c()

for (i in d %>% select(a) %>% distinct() %>% pull(a)){

auc[i] = with(d %>% filter(a==i),simple_auc(rel,fp))

}

auc[order(-auc)]



# 

d %>% 
  filter(t<9100) %>% 
  group_by(a) %>% 
  summarise(m=summary(lm(rel~t))$coefficients[2,1]) %>% arrange(-m)

