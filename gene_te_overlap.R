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

species_te_df<-rbind(pancheri.cactaTIR,
                     pancheri.copiaLTR,
                     pancheri.gypsyLTR,
                     pancheri.harbingerTIR,
                     pancheri.helitron,
                     pancheri.marinerTIR,
                     pancheri.mutatorTIR,
                     revolutissima.cactaTIR,
                     revolutissima.copiaLTR,
                     revolutissima.gypsyLTR,
                     revolutissima.harbingerTIR,
                     revolutissima.helitron,
                     revolutissima.marinerTIR,
                     revolutissima.mutatorTIR,
                     vieillardii.cactaTIR,
                     vieillardii.copiaLTR,
                     vieillardii.gypsyLTR,
                     vieillardii.harbingerTIR,
                     vieillardii.helitron,
                     vieillardii.marinerTIR,
                     vieillardii.mutatorTIR,
                     yahouensis.cactaTIR,
                     yahouensis.copiaLTR,
                     yahouensis.gypsyLTR,
                     yahouensis.harbingerTIR,
                     yahouensis.helitron,
                     yahouensis.marinerTIR,
                     yahouensis.mutatorTIR,
                     impolita.cactaTIR,
                     impolita.copiaLTR,
                     impolita.gypsyLTR,
                     impolita.harbingerTIR,
                     impolita.helitron,
                     impolita.marinerTIR,
                     impolita.mutatorTIR)







# a function which takes a gene ID <> GO terms mapping object per species, filters genes of interest, and counts how many time each GO term occurs
get_go_counts <- function(mp, data_frame, selection_statement){
  query_df <- data_frame %>% filter(rlang::eval_tidy(rlang::parse_expr(selection_statement)))
  mp %>% keep(names(.) %in% query_df$X1) %>% unlist() %>% table() %>% data.frame()
}



pancheri_table<-get_go_counts(mp_pancheri, species_te_df, "species == 'pancheri' & te_class == 'gypsyLTR'")
revolutissima_table<-get_go_counts(mp_revolutissima, species_te_df, "species == 'revolutissima' & te_class == 'gypsyLTR'")
vieillardii_table<-get_go_counts(mp_vieillardii, species_te_df, "species == 'vieillardii' & te_class == 'gypsyLTR'")
yahouensis_table<-get_go_counts(mp_yahouensis, species_te_df, "species == 'yahouensis' & te_class == 'gypsyLTR'")
impolita_table<-get_go_counts(mp_impolita, species_te_df, "species == 'impolita' & te_class == 'gypsyLTR'")


test<-purrr::reduce(list(pancheri_table,
                   revolutissima_table,
                   vieillardii_table,
                   yahouensis_table,
                   impolita_table), dplyr::inner_join, by = ".") %>%
  set_colnames(c("GO", "pancheri", "revolutissima", "vieillardii", "yahouensis", "impolita"))
               
              
test$stdev<- test %>% dplyr::select(pancheri, revolutissima, vieillardii, yahouensis, impolita) %>% apply(1, sd)

test$stdev %>% hist()

test %>% filter(stdev > 100) %>% dplyr::select(GO)

test %>% filter(stdev > 100) %>% dplyr::select(pancheri, revolutissima, vieillardii, yahouensis, impolita) %>% boxplot()






mp_impolita[grep("GO:0002215", mp_impolita)] %>% names()


mp_impolita[["g13983.t1"]]

mp_impolita
mp_vieillardii
mp_pancheri
mp_revolutissima
mp_yahouensis


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



testing4<-intersect(impolita.gypsyLTR_go$result$GO.ID, revolutissima.gypsyLTR_go$result$GO.ID)

get_de_genes_in_term(impolita.gypsyLTR$X1, "GO:0016192", impolita.gypsyLTR_go$goData) %>% data.frame()

sapply(testing4, function(x) get_de_genes_in_term(impolita.gypsyLTR$X1, x, impolita.gypsyLTR_go$goData) %>% data.frame())





# how old are the TEs close to genes here?



get_de_genes_in_term(pancheri.gypsyLTR$X1, "GO:0048544", pancheri.gypsyLTR_go$goData) %>% data.frame()
get_de_genes_in_term(revolutissima.gypsyLTR$X1, "GO:0048544", revolutissima.gypsyLTR_go$goData) %>% data.frame()

genes2test<-get_de_genes_in_term(impolita.gypsyLTR$X1, "GO:0016192", impolita.gypsyLTR_go$goData) %>% data.frame() %>% pull() %>% str_split_i("\\.", 1)


testing2<-read.table("impolita_gene_te_dists") %>% set_colnames(c("geneid", "basepairs", "classification", "insertion"))

testing2 %>% filter(colour_by == "yes") %>% pull(basepairs) %>% hist()



testing2 %>% 
  #mutate(colour_by=ifelse(geneid %in% genes2test, "yes", "no")) %>% 
  #filter(classification %in% c("LTR/Copia", "LTR/Gypsy")) %>% 
  arrange(desc(colour_by)) %>% 
  filter(abs(basepairs) < 25000 & classification %in% c("LTR/Copia", "LTR/Gypsy")) %>% 
  drop_na() %>%
  ggplot(aes(x=basepairs, y=insertion, colour=colour_by, alpha=colour_by)) +
  geom_point(size=1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~classification) + 
  scale_colour_manual(values = c("red", "grey89")) +
  scale_alpha_manual(values=c(1, 0.1))



ggplot(testing3, aes(x=basepairs, y=insertion)) +
  geom_point(size=0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

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


all_ltr %>% group_by(species, class) %>%
  summarise(count=n())

all_ltr$species <- factor(all_ltr$species, levels=c("yahouensis", "pancheri", "impolita", "revolutissima", "vieillardii"))



pdf("fuckssake.pdf", height=6, width=6)
all_ltr %>% 
  ggplot(aes(x = V1, fill = class)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.25) +
  facet_wrap(~ species, ncol=2)
dev.off()

rbind(read.table("pancheri.gene_gypsyLTR_window1000_insertiondates") %>% mutate(species="pancheri_1000", class="gypsy"), 
      read.table("revolutissima.gene_gypsyLTR_window1000_insertiondates") %>% mutate(species="revolutissima_1000", class="gypsy")) %>%
  ggplot(aes(x = V1, fill = class)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.25) +
  facet_wrap(~ species, ncol=2)






