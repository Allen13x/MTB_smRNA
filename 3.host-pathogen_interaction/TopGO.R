# R <- version 4.0.3

########library
library(topGO)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)

########load file 
expredata <- read_delim('filtered_couples.csv',delim=',')



########Pull genes list
transcripts <- expredata %>%separate(mrna,c('mrna'),sep='\\.') %>% pull(mrna)


#########Transcript to genes

myMart<- useMart('ensembl',dataset='hsapiens_gene_ensembl')

annotations <- getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id','entrezgene_id','external_gene_name','description'),
                     values =  transcripts,
                     filters = 'ensembl_transcript_id',mart = myMart)

genes<-annotations$ensembl_gene_id

########BP Mapping (on biological process ontology)
##GO2gene mapping
BPmapping <- annFUN.org(whichOnto = 'BP',mapping = 'org.Hs.eg.db',ID = 'ensembl')

##gene2GO mapping
BPmapping_gene2go <- inverseList(BPmapping)


########Universe (term of comparisons)
geneUniverse <- unique(unlist(BPmapping))
geneList <- factor(as.integer(geneUniverse%in%genes))
names(geneList)<-geneUniverse


########Top Go main object

GOdata_BP <- new('topGOdata',ontology='BP',
                 allGenes = geneList,
                 nodeSize = 5,
                 annot = annFUN.org,
                 mapping = 'org.Hs.eg.db',
                 ID = 'ensembl')




resultFisher.classic_BP <- runTest(GOdata_BP, algorithm = 'classic',statistic = 'fisher')

resultFisher.elim_BP <- runTest(GOdata_BP, algorithm = 'elim',statistic = 'fisher')



allRes_BP2 <- GenTable(GOdata_BP,
                       elimFisher = resultFisher.elim_BP,
                       classicFisher = resultFisher.classic_BP,
                       orderBy = 'elimFisher',
                       #ranksOf = 'elimFisher',
                       topNodes=5000,
                       numChar=100)



##### Top 20 go terms

fisherBP <- allRes_BP2 %>%  mutate(fractionGenes = Significant/Annotated*100) %>% dplyr::select(Term, Annotated,fractionGenes,Significant,elimFisher )

elim_BP_top <- fisherBP %>%  filter(elimFisher<0.01)




elim_BP_top %>% mutate(elimFisher=as.numeric(elimFisher)) %>% arrange(elimFisher) %>%head(20) %>% 
  ggplot(aes(x=Significant,y=reorder(Term, -elimFisher)))+
  geom_point(aes(color=elimFisher, size=fractionGenes))+
  scale_x_log10()+
  theme_bw()+
  xlab('Significant Genes')+
  #scale_color_viridis_c()
  scale_color_gradient(low='purple4',high='deepskyblue1')



##### GraphNodes



printGraph(object = GOdata_BP,result = resultFisher.elim_BP,firstSigNodes=5,pdfSW=T)