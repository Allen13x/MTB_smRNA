
# R <- Version 4.0.2

#Import Libraries
library(dplyr)
library(tidyr)
library(ggplot2)

##### Import Data #####

#### IntaRNA_sTAR predictions

### sRNA ryhB
######## Predictions

tar_r<- read.delim('Predictions/res_tar.csv',header = T,sep = ';',stringsAsFactors = F) %>% select(id1,id2,E) %>% mutate(alg='sTar')

mir_r<- read.delim('Predictions/miranda_r.csv',header=T,sep=';',stringsAsFactors = F) %>% rename(id2=Seq1,id1=Seq2,E=Max.Energy) %>% select(id1,id2,E) %>% mutate(alg='mir')

######## Validated Target

ver_r<- c('b0118','b0156','b0592','b0683','b0721',
          'b0722','b1531','b1612','b1656','b1778',
          'b1981','b2155','b2530','b3365','b3510',
          'b3607','b3928','b4154')

######## Validated non-Target

fp_r <- c('b3336','b1683','b4367','b1684',
          'b0150','b1905','b2832')

### sRNA sgrS
######## Predictions
tar_s<- read.delim('Predictions/res_tar_s.csv',header = T,sep = ';',stringsAsFactors = F)%>% mutate(id2='sgrS')%>% select(id1,id2,E) %>% mutate(alg='sTar')

mir_s<- read.delim('Predictions/miranda_s.csv',header=T,sep=';',stringsAsFactors = F) %>% rename(id2=Seq1,id1=Seq2,E=Max.Energy) %>% select(id1,id2,E)%>% mutate(alg='mir')

######## Validated Target

ver_s <- c('b1101','b1658','b1817','b1818',
         'b2153','b3433','b3826','b4116')

######## Validated non-Target

fp_s <- c('b1973','b0573','b0296','b4506','b1857','b1856',
        'b1178','b2012','b1858','b1729','b4035','b2150')



######## Join predictions dataset
tar= rbind(tar_r,tar_s)



######## Verified dataset target
verified <- data.frame(v=c(rep('t',18),rep('f',7),rep('t',8),rep('f',12)),
                       srna=c(rep('ryhB',25),rep('sgrS',20)),
                       target=c(ver_r,fp_r,ver_s,fp_s))




#### Check TP and FP rates - single sRNAs#####


d<- data.frame(t=c('a'),srn=c('a'),rel=c('1'),fp=c('1'))


for (t in seq(0,2300,100)){
  for (i in c('ryhB','sgrS')) {


  n<- length(verified %>% filter(srna==i,v=='t') %>% pull(srna))
  f<- length(verified %>% filter(srna==i,v=='f') %>% pull(srna))
  m<-tar %>% filter(id2==i) %>%mutate(r=rank(E))%>% filter(id1%in%(verified %>% filter(srna==i,v=='t') %>% pull(target)),r<t) %>% summarise(n=n()) %>% pull(n)
  z<-tar %>% filter(id2==i) %>%mutate(r=rank(E))%>% filter(id1%in%(verified %>% filter(srna==i,v=='f') %>% pull(target)),r<t) %>% summarise(n=n()) %>% pull(n)
  d<-rbind(d,c(t,i,m/n,z/f))
  
}
}



d<- d[2:49,]


##### ROC figures - single sRNAs

d %>% mutate_at(.vars = vars(rel,fp,t),as.numeric) %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'red',high = 'blue' )+
  ggtitle('Roc curve - IntaRNA_sTar single sRNA')+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~srn)


ggsave('Fig/roc_Sgrs_facet.png', dpi=300)



#### Check TP and FP rates - overall sRNAs#####



d<- data.frame(t=c('a'),srn=c('a'),rel=c('1'),fp=c('1'))


for (t in seq(0,2300,100)){

    vt<-verified %>% filter(v=='t')
    vf<-verified %>% filter(v=='f')
    
    n<- length(vt %>% pull(srna))
    f<- length(vf %>% pull(srna))
    m<-tar  %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vt,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    z<-tar %>%mutate(r=rank(E))%>% filter(r<t) %>%
      inner_join(vf,by=c('id2'='srna','id1'='target')) %>% 
      summarise(n=n()) %>% pull(n)
    d<-rbind(d,c(t,'srna',m/n,z/f))
    

}





d=d[2:25,]

##### ROC figures - single sRNAs


d %>% mutate_at(.vars = vars(rel,fp,t),as.numeric) %>% 
  ggplot(aes(x=fp,y=rel))+
  ylim(0,1)+xlim(0,1)+
  geom_point(aes(color=t),lwd=3)+
  geom_path()+
  theme_classic()+
  ylab('True Positive')+
  xlab('False Positive')+
  scale_color_gradient(name='Threshold - Top Targets',low = 'red',high = 'blue' )+
  ggtitle('Roc curve - IntaRNA_sTar_joint')+
  geom_abline(slope = 1,intercept = 0)



ggsave('Fig/roc_Sgrs.png', dpi=300)
