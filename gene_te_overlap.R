library(topGO)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)


############################################
#        function for topGO enrichment     #
############################################

get_enriched_terms<-function(gene_list, mappings, return_sample_GOData=FALSE){
  # use the gene 2 GOterms mapping provided for D. incarnata
  geneID2GO<-mappings
  # the input genes form the input, use these to annotate all genes, 1 is present in input list, 0 is absent
  geneSel<-gene_list
  geneSel<-factor(as.integer(names(geneID2GO) %in% geneSel))
  names(geneSel)<-names(geneID2GO)
  
  # set up the topGO object
  sampleGOdata <- new("topGOdata",
                      ontology = "BP",
                      allGenes = geneSel, 
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)
  
  # run three tests, fisher, Kol-Smirn, and Kol-Smirn with elimination
  resultFisher <- runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
  
  # generate summary tane and return it
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100,
                     numChar=1000 )
  #allRes<-GenTable(sampleGOdata, Fis = resultFisher, topNodes = 20)
  
  if (return_sample_GOData == TRUE){
    return(list(result=allRes, goData=sampleGOdata))
  } else {
    return(allRes)
  }
}

##############################################################
#       function to get genes underlying enriched terms      #
##############################################################

get_de_genes_in_term<-function(degs, go_term, go_data){
  genes_in_term<-genesInTerm(go_data, go_term)[[1]]
  degs_in_term<-intersect(genes_in_term, degs)
  return(degs_in_term)
}

###############################################
#       function to filter topGO object        #
###############################################

filter_topGO<-function(topgo_object){
  return(topgo_object$result %>% filter(classicFisher < 0.05))
}


############################################
#        load in the GO ID mappings        #
############################################

# for the GO term enrichment tests
mp_impolita<-readMappings("impolita_topGO_annotation.txt")
mp_vieillardii<-readMappings("vieillardii_topGO_annotation.txt")
mp_pancheri<-readMappings("pancheri_topGO_annotation.txt")
mp_revolutissima<-readMappings("revolutissima_topGO_annotation.txt")
mp_yahouensis<-readMappings("yahouensis_topGO_annotation.txt")



#############################################################
#    read in annotation of genes within 5KB of a TE         #
#############################################################


pancheri.cactaTIR<-read_delim("pancheri/pancheri.gene_cactaTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="cactaTIR")
pancheri.copiaLTR<-read_delim("pancheri/pancheri.gene_copiaLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="copiaLTR")
pancheri.gypsyLTR<-read_delim("pancheri/pancheri.gene_gypsyLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="gypsyLTR")
pancheri.harbingerTIR<-read_delim("pancheri/pancheri.gene_harbingerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="harbingerTIR")
pancheri.helitron<-read_delim("pancheri/pancheri.gene_helitron_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="helitron")
pancheri.marinerTIR<-read_delim("pancheri/pancheri.gene_marinerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="marinerTIR")
pancheri.mutatorTIR<-read_delim("pancheri/pancheri.gene_mutatorTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="mutatorTIR")

revolutissima.cactaTIR<-read_delim("revolutissima/revolutissima.gene_cactaTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="cactaTIR")
revolutissima.copiaLTR<-read_delim("revolutissima/revolutissima.gene_copiaLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="copiaLTR")
revolutissima.gypsyLTR<-read_delim("revolutissima/revolutissima.gene_gypsyLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="gypsyLTR")
revolutissima.harbingerTIR<-read_delim("revolutissima/revolutissima.gene_harbingerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="harbingerTIR")
revolutissima.helitron<-read_delim("revolutissima/revolutissima.gene_helitron_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="helitron")
revolutissima.marinerTIR<-read_delim("revolutissima/revolutissima.gene_marinerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="marinerTIR")
revolutissima.mutatorTIR<-read_delim("revolutissima/revolutissima.gene_mutatorTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="mutatorTIR")

vieillardii.cactaTIR<-read_delim("viellardiei/vieillardii.gene_cactaTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="cactaTIR")
vieillardii.copiaLTR<-read_delim("viellardiei/vieillardii.gene_copiaLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="copiaLTR")
vieillardii.gypsyLTR<-read_delim("viellardiei/vieillardii.gene_gypsyLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="gypsyLTR")
vieillardii.harbingerTIR<-read_delim("viellardiei/vieillardii.gene_harbingerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="harbingerTIR")
vieillardii.helitron<-read_delim("viellardiei/vieillardii.gene_helitron_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="helitron")
vieillardii.marinerTIR<-read_delim("viellardiei/vieillardii.gene_marinerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="marinerTIR")
vieillardii.mutatorTIR<-read_delim("viellardiei/vieillardii.gene_mutatorTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="mutatorTIR")

yahouensis.cactaTIR<-read_delim("yahouensis/yahouensis.gene_cactaTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="cactaTIR")
yahouensis.copiaLTR<-read_delim("yahouensis/yahouensis.gene_copiaLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="copiaLTR")
yahouensis.gypsyLTR<-read_delim("yahouensis/yahouensis.gene_gypsyLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="gypsyLTR")
yahouensis.harbingerTIR<-read_delim("yahouensis/yahouensis.gene_harbingerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="harbingerTIR")
yahouensis.helitron<-read_delim("yahouensis/yahouensis.gene_helitron_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="helitron")
yahouensis.marinerTIR<-read_delim("yahouensis/yahouensis.gene_marinerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="marinerTIR")
yahouensis.mutatorTIR<-read_delim("yahouensis/yahouensis.gene_mutatorTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="mutatorTIR")


impolita.cactaTIR<-read_delim("impolita/impolita.gene_cactaTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="cactaTIR")
impolita.copiaLTR<-read_delim("impolita/impolita.gene_copiaLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="copiaLTR")
impolita.gypsyLTR<-read_delim("impolita/impolita.gene_gypsyLTR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="gypsyLTR")
impolita.harbingerTIR<-read_delim("impolita/impolita.gene_harbingerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="harbingerTIR")
impolita.helitron<-read_delim("impolita/impolita.gene_helitron_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="helitron")
impolita.marinerTIR<-read_delim("impolita/impolita.gene_marinerTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="marinerTIR")
impolita.mutatorTIR<-read_delim("impolita/impolita.gene_mutatorTIR_window1000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="mutatorTIR")

########################################################################
#       get the go enrichments and objects of closely related genes    #
########################################################################


pancheri.cactaTIR_go<-get_enriched_terms(pancheri.cactaTIR$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.copiaLTR_go<-get_enriched_terms(pancheri.copiaLTR$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.gypsyLTR_go<-get_enriched_terms(pancheri.gypsyLTR$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.harbingerTIR_go<-get_enriched_terms(pancheri.harbingerTIR$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.helitron_go<-get_enriched_terms(pancheri.helitron$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.marinerTIR_go<-get_enriched_terms(pancheri.marinerTIR$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.mutatorTIR_go<-get_enriched_terms(pancheri.mutatorTIR$X1, mp_pancheri, return_sample_GOData=TRUE)


revolutissima.cactaTIR_go<-get_enriched_terms(revolutissima.cactaTIR$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.copiaLTR_go<-get_enriched_terms(revolutissima.copiaLTR$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.gypsyLTR_go<-get_enriched_terms(revolutissima.gypsyLTR$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.harbingerTIR_go<-get_enriched_terms(revolutissima.harbingerTIR$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.helitron_go<-get_enriched_terms(revolutissima.helitron$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.marinerTIR_go<-get_enriched_terms(revolutissima.marinerTIR$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.mutatorTIR_go<-get_enriched_terms(revolutissima.mutatorTIR$X1, mp_revolutissima, return_sample_GOData=TRUE)

vieillardii.cactaTIR_go<-get_enriched_terms(vieillardii.cactaTIR$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.copiaLTR_go<-get_enriched_terms(vieillardii.copiaLTR$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.gypsyLTR_go<-get_enriched_terms(vieillardii.gypsyLTR$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.harbingerTIR_go<-get_enriched_terms(vieillardii.harbingerTIR$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.helitron_go<-get_enriched_terms(vieillardii.helitron$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.marinerTIR_go<-get_enriched_terms(vieillardii.marinerTIR$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.mutatorTIR_go<-get_enriched_terms(vieillardii.mutatorTIR$X1, mp_vieillardii, return_sample_GOData=TRUE)


yahouensis.copiaLTR_go<-get_enriched_terms(yahouensis.copiaLTR$X1, mp_yahouensis, return_sample_GOData=TRUE)
yahouensis.gypsyLTR_go<-get_enriched_terms(yahouensis.gypsyLTR$X1, mp_yahouensis, return_sample_GOData=TRUE)


impolita.cactaTIR_go<-get_enriched_terms(impolita.cactaTIR$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.copiaLTR_go<-get_enriched_terms(impolita.copiaLTR$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.gypsyLTR_go<-get_enriched_terms(impolita.gypsyLTR$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.harbingerTIR_go<-get_enriched_terms(impolita.harbingerTIR$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.marinerTIR_go<-get_enriched_terms(impolita.marinerTIR$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.mutatorTIR_go<-get_enriched_terms(impolita.mutatorTIR$X1, mp_impolita, return_sample_GOData=TRUE)







# yahou-panch sister
intersect(filter_topGO(yahouensis.gypsyLTR_go)$Term, filter_topGO(pancheri.gypsyLTR_go)$Term)

# yahou-impo same soil (volcanic)
intersect(filter_topGO(yahouensis.gypsyLTR_go)$Term, filter_topGO(impolita.gypsyLTR_go)$Term)

# impo-revo sister
intersect(impolita.gypsyLTR_go$result$Term, revolutissima.gypsyLTR_go$result$Term)

# panch-revo same soil (ultramafic)
intersect(pancheri.gypsyLTR_go$result$Term, revolutissima.gypsyLTR_go$result$Term)


# how old are the TEs close to genes here?



genesInTerm(yahouensis.gypsyLTR5000_go$goData, c("GO:0042127", "GO:0001522"))
get_de_genes_in_term(yahouensis.gypsyLTR5000$X1, "GO:0042127", yahouensis.gypsyLTR5000_go$goData)

get_de_genes_in_term(pancheri.gypsyLTR5000$X1, "GO:0015689", pancheri.gypsyLTR5000_go$goData)
get_de_genes_in_term(revolutissima.gypsyLTR5000$X1, "GO:0015689", revolutissima.gypsyLTR5000_go$goData)






vieillardii.gene_copiaLTR<-read.table("to_local/vieillardii.gene_copiaLTR") %>% mutate(species="vieillardii", class="copia")
revolutissima.gene_copiaLTR<-read.table("to_local/revolutissima.gene_copiaLTR") %>% mutate(species="revolutissima", class="copia")
impolita.gene_gypsyLTR<-read.table("to_local/impolita.gene_gypsyLTR") %>% mutate(species="impolita", class="gypsy")
pancheri.gene_gypsyLTR<-read.table("to_local/pancheri.gene_gypsyLTR") %>% mutate(species="pancheri", class="gypsy")
yahouensis.gene_copiaLTR<-read.table("to_local/yahouensis.gene_copiaLTR") %>% mutate(species="yahouensis", class="copia")
yahouensis.gene_gypsyLTR<-read.table("to_local/yahouensis.gene_gypsyLTR") %>% mutate(species="yahouensis", class="gypsy")
pancheri.gene_copiaLTR<-read.table("to_local/pancheri.gene_copiaLTR") %>% mutate(species="pancheri", class="copia")
impolita.gene_copiaLTR<-read.table("to_local/impolita.gene_copiaLTR") %>% mutate(species="impolita", class="copia")
revolutissima.gene_gypsyLTR<-read.table("to_local/revolutissima.gene_gypsyLTR") %>% mutate(species="revolutissima", class="gypsy")
vieillardii.gene_gypsyLTR<-read.table("to_local/vieillardii.gene_gypsyLTR") %>% mutate(species="vieillardii", class="gypsy")


all_ltr<-rbind(vieillardii.gene_copiaLTR,
      revolutissima.gene_copiaLTR,
      impolita.gene_gypsyLTR,
      pancheri.gene_gypsyLTR,
      yahouensis.gene_copiaLTR,
      yahouensis.gene_gypsyLTR,
      pancheri.gene_copiaLTR,
      impolita.gene_copiaLTR,
      revolutissima.gene_gypsyLTR,
      vieillardii.gene_gypsyLTR)

all_ltr %>% group_by(class) %>% summarise_at("species", sum) 

all_ltr %>% group_by(species, class) %>%
  summarise(count=n())

all_ltr$species <- factor(all_ltr$species, levels=c("yahouensis", "pancheri", "impolita", "revolutissima", "vieillardii"))


all_ltr %>% 
  ggplot(aes(x = V1, fill = class)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.25) +
  facet_wrap(~ species, ncol=2)








