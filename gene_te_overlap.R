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


pancheri.cactaTIR5000<-read_delim("pancheri/pancheri.gene_cactaTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="cactaTIR")
pancheri.copiaLTR5000<-read_delim("pancheri/pancheri.gene_copiaLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="copiaLTR")
pancheri.gypsyLTR5000<-read_delim("pancheri/pancheri.gene_gypsyLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="gypsyLTR")
pancheri.harbingerTIR5000<-read_delim("pancheri/pancheri.gene_harbingerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="harbingerTIR")
pancheri.helitron5000<-read_delim("pancheri/pancheri.gene_helitron_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="helitron")
pancheri.marinerTIR5000<-read_delim("pancheri/pancheri.gene_marinerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="marinerTIR")
pancheri.mutatorTIR5000<-read_delim("pancheri/pancheri.gene_mutatorTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="pancheri", te_class="mutatorTIR")

revolutissima.cactaTIR5000<-read_delim("revolutissima/revolutissima.gene_cactaTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="cactaTIR")
revolutissima.copiaLTR5000<-read_delim("revolutissima/revolutissima.gene_copiaLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="copiaLTR")
revolutissima.gypsyLTR5000<-read_delim("revolutissima/revolutissima.gene_gypsyLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="gypsyLTR")
revolutissima.harbingerTIR5000<-read_delim("revolutissima/revolutissima.gene_harbingerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="harbingerTIR")
revolutissima.helitron5000<-read_delim("revolutissima/revolutissima.gene_helitron_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="helitron")
revolutissima.marinerTIR5000<-read_delim("revolutissima/revolutissima.gene_marinerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="marinerTIR")
revolutissima.mutatorTIR5000<-read_delim("revolutissima/revolutissima.gene_mutatorTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="revolutissima", te_class="mutatorTIR")

vieillardii.cactaTIR5000<-read_delim("viellardiei/vieillardii.gene_cactaTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="cactaTIR")
vieillardii.copiaLTR5000<-read_delim("viellardiei/vieillardii.gene_copiaLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="copiaLTR")
vieillardii.gypsyLTR5000<-read_delim("viellardiei/vieillardii.gene_gypsyLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="gypsyLTR")
vieillardii.harbingerTIR5000<-read_delim("viellardiei/vieillardii.gene_harbingerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="harbingerTIR")
vieillardii.helitron5000<-read_delim("viellardiei/vieillardii.gene_helitron_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="helitron")
vieillardii.marinerTIR5000<-read_delim("viellardiei/vieillardii.gene_marinerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="marinerTIR")
vieillardii.mutatorTIR5000<-read_delim("viellardiei/vieillardii.gene_mutatorTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="vieillardii", te_class="mutatorTIR")

yahouensis.cactaTIR5000<-read_delim("yahouensis/yahouensis.gene_cactaTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="cactaTIR")
yahouensis.copiaLTR5000<-read_delim("yahouensis/yahouensis.gene_copiaLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="copiaLTR")
yahouensis.gypsyLTR5000<-read_delim("yahouensis/yahouensis.gene_gypsyLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="gypsyLTR")
yahouensis.harbingerTIR5000<-read_delim("yahouensis/yahouensis.gene_harbingerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="harbingerTIR")
yahouensis.helitron5000<-read_delim("yahouensis/yahouensis.gene_helitron_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="helitron")
yahouensis.marinerTIR5000<-read_delim("yahouensis/yahouensis.gene_marinerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="marinerTIR")
yahouensis.mutatorTIR5000<-read_delim("yahouensis/yahouensis.gene_mutatorTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="yahouensis", te_class="mutatorTIR")


impolita.cactaTIR5000<-read_delim("impolita/impolita.gene_cactaTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="cactaTIR")
impolita.copiaLTR5000<-read_delim("impolita/impolita.gene_copiaLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="copiaLTR")
impolita.gypsyLTR5000<-read_delim("impolita/impolita.gene_gypsyLTR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="gypsyLTR")
impolita.harbingerTIR5000<-read_delim("impolita/impolita.gene_harbingerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="harbingerTIR")
impolita.helitron5000<-read_delim("impolita/impolita.gene_helitron_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="helitron")
impolita.marinerTIR5000<-read_delim("impolita/impolita.gene_marinerTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="marinerTIR")
impolita.mutatorTIR5000<-read_delim("impolita/impolita.gene_mutatorTIR_window5000.annotation", col_names = FALSE) %>% mutate(species="impolita", te_class="mutatorTIR")

########################################################################
#       get the go enrichments and objects of closely related genes    #
########################################################################


pancheri.cactaTIR5000_go<-get_enriched_terms(pancheri.cactaTIR5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.copiaLTR5000_go<-get_enriched_terms(pancheri.copiaLTR5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.gypsyLTR5000_go<-get_enriched_terms(pancheri.gypsyLTR5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.harbingerTIR5000_go<-get_enriched_terms(pancheri.harbingerTIR5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.helitron5000_go<-get_enriched_terms(pancheri.helitron5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.marinerTIR5000_go<-get_enriched_terms(pancheri.marinerTIR5000$X1, mp_pancheri, return_sample_GOData=TRUE)
pancheri.mutatorTIR5000_go<-get_enriched_terms(pancheri.mutatorTIR5000$X1, mp_pancheri, return_sample_GOData=TRUE)


revolutissima.cactaTIR5000_go<-get_enriched_terms(revolutissima.cactaTIR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.copiaLTR5000_go<-get_enriched_terms(revolutissima.copiaLTR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.gypsyLTR5000_go<-get_enriched_terms(revolutissima.gypsyLTR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.harbingerTIR5000_go<-get_enriched_terms(revolutissima.harbingerTIR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.helitron5000_go<-get_enriched_terms(revolutissima.helitron5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.marinerTIR5000_go<-get_enriched_terms(revolutissima.marinerTIR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)
revolutissima.mutatorTIR5000_go<-get_enriched_terms(revolutissima.mutatorTIR5000$X1, mp_revolutissima, return_sample_GOData=TRUE)

vieillardii.cactaTIR5000_go<-get_enriched_terms(vieillardii.cactaTIR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.copiaLTR5000_go<-get_enriched_terms(vieillardii.copiaLTR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.gypsyLTR5000_go<-get_enriched_terms(vieillardii.gypsyLTR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.harbingerTIR5000_go<-get_enriched_terms(vieillardii.harbingerTIR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.helitron5000_go<-get_enriched_terms(vieillardii.helitron5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.marinerTIR5000_go<-get_enriched_terms(vieillardii.marinerTIR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)
vieillardii.mutatorTIR5000_go<-get_enriched_terms(vieillardii.mutatorTIR5000$X1, mp_vieillardii, return_sample_GOData=TRUE)


yahouensis.copiaLTR5000_go<-get_enriched_terms(yahouensis.copiaLTR5000$X1, mp_yahouensis, return_sample_GOData=TRUE)
yahouensis.gypsyLTR5000_go<-get_enriched_terms(yahouensis.gypsyLTR5000$X1, mp_yahouensis, return_sample_GOData=TRUE)


impolita.cactaTIR5000_go<-get_enriched_terms(impolita.cactaTIR5000$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.copiaLTR5000_go<-get_enriched_terms(impolita.copiaLTR5000$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.gypsyLTR5000_go<-get_enriched_terms(impolita.gypsyLTR5000$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.harbingerTIR5000_go<-get_enriched_terms(impolita.harbingerTIR5000$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.marinerTIR5000_go<-get_enriched_terms(impolita.marinerTIR5000$X1, mp_impolita, return_sample_GOData=TRUE)
impolita.mutatorTIR5000_go<-get_enriched_terms(impolita.mutatorTIR5000$X1, mp_impolita, return_sample_GOData=TRUE)


yahouensis.gypsyLTR5000_go$result %>% filter(classicFisher < 0.05)

yahouensis.gypsyLTR5000_go$goData


genesInTerm(yahouensis.gypsyLTR5000_go$goData, c("GO:0042127", "GO:0001522"))

get_de_genes_in_term(yahouensis.gypsyLTR5000$X1, "GO:0042127", yahouensis.gypsyLTR5000_go$goData)

# yahou-panch sister
intersect(filter_topGO(yahouensis.gypsyLTR5000_go)$Term, filter_topGO(pancheri.gypsyLTR5000_go)$Term)

# yahou-impo same soil (volcanic)
intersect(filter_topGO(yahouensis.gypsyLTR5000_go)$Term, filter_topGO(impolita.gypsyLTR5000_go)$Term)

# impo-revo sister
intersect(impolita.gypsyLTR5000_go$result$Term, revolutissima.gypsyLTR5000_go$result$Term)

# panch-revo same soil (ultramafic)
intersect(pancheri.gypsyLTR5000_go$result$Term, revolutissima.gypsyLTR5000_go$result$Term)


get_de_genes_in_term(pancheri.gypsyLTR5000$X1, "GO:0015689", pancheri.gypsyLTR5000_go$goData)
get_de_genes_in_term(revolutissima.gypsyLTR5000$X1, "GO:0015689", revolutissima.gypsyLTR5000_go$goData)





