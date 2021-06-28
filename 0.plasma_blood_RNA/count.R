library(tidyverse)
library(DESeq2)
library(readxl)
library(xlsx)
library(ggpubr)
library(pheatmap)
tpm3 <- function(counts,len) {
	  x <- counts/(len/1000)
  return(t(t(x)*1e6/colSums(x)))
}

length <- read_delim('counts/mtb_count.tsv',delim='\t',comment = '#')$Length
mtb_gene <- read_delim('counts/mtb_count.tsv',delim='\t',comment = '#')$Geneid


read_delim('counts/mtb_count.tsv',delim='\t',comment = '#') %>% head()
  select(-Start,-End,-Strand,-Chr) %>% 
	    gather(s,c,-Geneid,-Length) %>% 
	      mutate(s=str_replace(s,'-','_')) %>% 
	        extract(s,c('s','dr'),regex='(.*)_(.*)$') %>%select(-dr) %>%
		  extract(s,c('dr','s'),regex='(.*)/(.*)$')%>%select(-dr) %>%
		    arrange(s) %>% 
		      spread(s,c) %>% 
		        column_to_rownames('Geneid') %>% 
			  replace(is.na(.), 0) -> mtb_counts




		  mtb_counts<-mtb_counts[mtb_keep,]


		  mtb_l<-mtb_counts$Length
		  mtb_counts<- mtb_counts[,-1]


		  mtb_tpm<-tpm3(mtb_counts,mtb_l)


mtb_tpm %>% 
	  group_by(Geneid,cat) %>%
	    summarise(m=sum(c>0,na.rm=T)) %>% spread(cat,m)%>%filter((CTRL==0 | ATB+LTBI>10),(ATB+LTBI>2 | ATB> 1))  %>%
	      filter(str_detect(Geneid,'candidate'))%>%write_delim('Final_candidates_counts-tsv',delim='\t')
