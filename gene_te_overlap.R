library(topGO)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(ape)
library(ggtree)
library(ggplot2)
library(egg)

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
mp_sandwicensis<-readMappings("sandwicensis_topGO_annotation.txt")


names(mp_impolita) %<>% str_split_i(".t1", 1)
names(mp_vieillardii) %<>% str_split_i(".t1", 1)
names(mp_pancheri) %<>% str_split_i(".t1", 1)
names(mp_revolutissima) %<>% str_split_i(".t1", 1)
names(mp_yahouensis) %<>% str_split_i(".t1", 1)
names(mp_sandwicensis) %<>% str_split_i(".t1", 1)



#######################################################################################
#     for each species, read in GO mappings and gene-TE proximity table               #
#     inner join them to get the transcript ID (arbitrary but need the tranbscript)   #
#      get enriched terms of the transcripts which have a TE in them (0 distance)
########################################################################################


  



#sandwicensis_GO<-read_delim("to_r/sandwicensis_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#sandwicensis.gene_te<-read_delim("to_r/sandwicensis.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
sandwicensis.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/sandwicensis.gene_te_dists", col_names = c("annotation", "gene", "gene_dist", "te_length", "insertion_date"))

#vieillardii_GO<-read_delim("to_r/vieillardii_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#vieillardii.gene_te<-read_delim("to_r/vieillardii.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
vieillardii.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/vieillardii.gene_te_dists", col_names = c("annotation", "gene", "gene_dist", "te_length", "insertion_date"))

#impolita_GO<-read_delim("to_r/impolita_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#impolita.gene_te<-read_delim("to_r/impolita.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
impolita.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/impolita.gene_te_dists", col_names = c("annotation", "gene", "gene_dist", "te_length", "insertion_date"))

#yahouensis_GO<-read_delim("to_r/yahouensis_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#yahouensis.gene_te<-read_delim("to_r/yahouensis.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
yahouensis.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/yahouensis.gene_te_dists", col_names = c("annotation", "delete", "gene", "gene_dist", "te_length", "insertion_date")) %>% dplyr::select(-delete)

#revolutissima_GO<-read_delim("to_r/revolutissima_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#revolutissima.gene_te<-read_delim("to_r/revolutissima.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
revolutissima.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/revolutissima.gene_te_dists", col_names = c("annotation", "gene", "gene_dist", "te_length", "insertion_date"))

#pancheri_GO<-read_delim("to_r/pancheri_GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
#pancheri.gene_te<-read_delim("to_r/pancheri.gene_te_dists_annotation", col_names = c("gene", "dist", "te", "ins", "none"))
pancheri.gene_te<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/gene_te_overlaps/pancheri.gene_te_dists", col_names = c("annotation", "gene", "gene_dist", "te_length", "insertion_date"))

all.gene_te<-rbind(sandwicensis.gene_te %>% mutate(species="D. sandwicensis"),
                   vieillardii.gene_te %>% mutate(species="D. vieillardii"),
                   impolita.gene_te %>% mutate(species="D. impolita"),
                   yahouensis.gene_te %>% mutate(species="D. yahouensis"),
                   revolutissima.gene_te %>% mutate(species="D. revolutissima"),
                   pancheri.gene_te %>% mutate(species="D. pancheri"))


pdf("test3.pdf", height=13, width=16)
all.gene_te %>% filter(annotation %in% c("LTR/Gypsy", "LTR/Copia")) %>% 
  ggplot(aes(x=te_length, fill=annotation)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25)  + 
  facet_wrap(~species) +
  geom_vline(xintercept=c(5300), linetype="dotted", colour="red", size=1) +
  geom_vline(xintercept=c(8500), linetype="dotted", colour="red", size=1) +
  geom_vline(xintercept=c(11550), linetype="dotted", colour="red", size=1) 
dev.off()


all.gene_te %>% filter(te_length > 8000 & te_length < 8700) %>%
  filter(abs(as.numeric(gene_dist)) < 10000 & annotation == "LTR/Gypsy") %>%
  ggplot(aes(x=gene_dist, fill=species)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) + facet_wrap(~species)
  

all.gene_te %>% filter(te_length > 11300 & te_length < 11700) %>%
  mutate(insertion_date=as.numeric(insertion_date)) %>%
  filter(annotation == "LTR/Gypsy") %>%
  ggplot(aes(x=insertion_date, fill=species)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) + facet_wrap(~species) +
  geom_vline(xintercept=c(0.973), linetype="dotted", colour="red", size=1) 


all.gene_te %>% filter(te_length > 8000 & te_length < 8700) %>%
  mutate(insertion_date=as.numeric(insertion_date)) %>%
  filter(annotation == "LTR/Gypsy") %>%
  ggplot(aes(x=insertion_date, fill=species)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) + facet_wrap(~species) +
  geom_vline(xintercept=c(0.96), linetype="dotted", colour="red", size=1) 





all.gene_te %>% ggplot(aes(x=te_length, fill=species)) + geom_histogram() + facet_wrap(~species)


##################################################
#                read in species tree            #
################################################## 

species_tree<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/lib1234_speciestree.editedTiplabs.nwk")
species_tree.rooted <- root(species_tree, which(species_tree$tip.label == "D.sandwicensis"))
tip<-c("D.olen", "D.fasciculosa", "D.macrocarpa", "D.ferrea")
species_tree.rooted<-drop.tip(species_tree.rooted, tip)
species_tree.rooted$species <- species_tree.rooted$tip.label
species_tree.rooted$tip.label<-str_replace(species_tree.rooted$tip.label, "D.", "D. ")


colours_tips <- case_when(species_tree.rooted$tip.label == "D. sandwicensis" ~ "D. sandwicensis",
                          species_tree.rooted$tip.label == "D. vieillardii" ~ "D. vieillardii",
                          species_tree.rooted$tip.label == "D. pancheri" ~ "D. pancheri",
                          species_tree.rooted$tip.label == "D. revolutissima" ~ "D. revolutissima",
                          species_tree.rooted$tip.label == "D. impolita" ~ "D. impolita",
                          species_tree.rooted$tip.label == "D. yahouensis" ~ "D. yahouensis",
                          !(species_tree.rooted$tip.label %in% c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. revolutissima", "D. impolita", "D. yahouensis")) ~ "No data")



dd <- data.frame(taxa=species_tree.rooted$tip.label, tipcols=colours_tips)
p<-ggtree(species_tree.rooted, size=1)
p <- p %<+% dd + geom_tippoint(aes(color=tipcols), size=10)
p2<-p + geom_tiplab(size=11, aes(color=tipcols), offset=0.002, show.legend=FALSE) + 
  scale_color_manual(values=c("goldenrod", "darkblue", "cornflowerblue", "darkolivegreen4", "#B25D91FF", "brown3"), limits = c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. revolutissima", "D. impolita", "D. yahouensis"), na.value = "grey77") + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none") +
  expand_limits(x = 0.07)






#########################################################################
#.    make hist of gene te distances across species with phylo tree.    #
#########################################################################

ho<-rbind(sandwicensis.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. sandwicensis"),
vieillardii.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. vieillardii"),
impolita.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. impolita"),
yahouensis.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. yahouensis"),
revolutissima.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. revolutissima"),
pancheri.gene_te %>% dplyr::select(gene_dist) %>% mutate(species="D. pancheri"))
ho$species <- factor(ho$species, levels=c("D. yahouensis", "D. pancheri",
                                            "D. impolita", "D. revolutissima",
                                            "D. vieillardii", "D. sandwicensis"))
ho <- ho %>%
  filter(abs(as.numeric(gene_dist)) < 10000) %>%
  ggplot(aes(x=gene_dist, fill=species)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) +
  facet_wrap(~species, ncol=2) +
  geom_vline(xintercept=c(0), linetype="dotted", colour="red", size=1) + 
  scale_fill_manual(values=c("goldenrod", "darkblue", "cornflowerblue", "darkolivegreen4", "#B25D91FF", "brown3"), limits = c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. revolutissima", "D. impolita", "D. yahouensis")) + 
  theme(legend.text = element_text(size=25),
        legend.position="none",
        strip.text = element_text(size = 25, face="italic"),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=25),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=25)) +
  xlab("Distance from gene (bp") +
  ylab("Number of LTR Retrotransposons")


pdf("test.pdf", width=20, height=15)
grid.arrange(p2, ho, ncol=2)
dev.off()

#######################################################################################
#.    plot paired scatter for TEs 0 distance and TEs > 0 distance for all species.    #
#######################################################################################


sample_data<-data.frame(counts=c(sandwicensis.gene_te %>% filter(gene_dist == 0) %>% nrow(), 
           sandwicensis.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow(),
           vieillardii.gene_te %>% filter(gene_dist == 0) %>% nrow(),
           vieillardii.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow(),
           impolita.gene_te %>% filter(gene_dist == 0) %>% nrow(),
           impolita.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow(),
           yahouensis.gene_te %>% filter(gene_dist == 0) %>% nrow(),
           yahouensis.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow(),
           revolutissima.gene_te %>% filter(gene_dist == 0) %>% nrow(),
           revolutissima.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow(),
           pancheri.gene_te %>% filter(gene_dist == 0) %>% nrow(),
           pancheri.gene_te %>% filter(gene_dist > 0 & gene_dist < 5000) %>% nrow()),
           category=rep(c("0", "1 - 5000"), 6),
           species=c("D. sandwicensis", "D. sandwicensis",
                     "D. vieillardii", "D. vieillardii",
                     "D. impolita", "D. impolita",
                     "D. yahouensis", "D. yahouensis",
                     "D. revolutissima", "D. revolutissima",
                     "D. pancheri", "D. pancheri")) %>%
  set_colnames(c("counts", "category", "species"))


sample_data$species <- factor(sample_data$species, levels=c("D. yahouensis",
                                                               "D. pancheri",
                                                               "D. impolita",
                                                               "D. revolutissima", 
                                                               "D. vieillardii",
                                                               "D. sandwicensis"))



distance_point<-ggplot(sample_data, aes(category, counts, fill=category)) + 
  geom_line(aes(group=species)) +
  geom_point(aes(colour=species), size=9) +
  theme(legend.text = element_text(size=25, face = "italic"),
        legend.title = element_blank(),
        strip.text = element_text(size = 25, face="italic"),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=25),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=25)) +
  scale_fill_discrete(guide="none") + 
  scale_colour_manual(values=c("brown3", "cornflowerblue", "#B25D91FF", "darkolivegreen4", "darkblue", "goldenrod")) +
  xlab("Distance from gene (bp)") +
  ylab("Number of LTR Retrotransposons")

pdf("test2.pdf", width=20, height=13)
grid.arrange(p2, distance_point, ncol=2)
dev.off()



#############################################################
#.    function to filter genes by their distance to a TE    #
#############################################################


filter_by_dist<-function(species.gene_te, selection_statement){
  to_return<-species.gene_te %>%
    filter(rlang::eval_tidy(rlang::parse_expr(selection_statement)))
  return(to_return)
}
sand_0<-filter_by_dist(sandwicensis.gene_te, "gene_dist == 0")
vie_0<-filter_by_dist(vieillardii.gene_te, "gene_dist == 0")
panc_0<-filter_by_dist(pancheri.gene_te, "gene_dist < 1000")
rev_0<-filter_by_dist(revolutissima.gene_te, "gene_dist == 0")
yah_0<-filter_by_dist(yahouensis.gene_te, "gene_dist == 0")
imp_0<-filter_by_dist(impolita.gene_te, "gene_dist == 0")

sand_0_GO<-get_enriched_terms(sand_0$gene, mp_sandwicensis, return_sample_GOData=TRUE)
vie_0_GO<-get_enriched_terms(vie_0$gene, mp_vieillardii, return_sample_GOData=TRUE)
panc_0_GO<-get_enriched_terms(panc_0$gene, mp_pancheri, return_sample_GOData=TRUE)
rev_0_GO<-get_enriched_terms(rev_0$gene, mp_revolutissima, return_sample_GOData=TRUE)
yah_0_GO<-get_enriched_terms(yah_0$gene, mp_yahouensis, return_sample_GOData=TRUE)
imp_0_GO<-get_enriched_terms(imp_0$gene, mp_impolita, return_sample_GOData=TRUE)

sand_0_GO$result %>% head(10)
vie_0_GO$result %>% head(10)
panc_0_GO$result %>% head(10)
rev_0_GO$result %>% head(10)
yah_0_GO$result %>% head(10)
imp_0_GO$result %>% head(10)




listInput<-list(sandwicensis = sand_0_GO$result %>% filter(classicFisher < 0.05) %>% pull(Term),
                pancheri = panc_0_GO$result %>% filter(classicFisher < 0.05) %>% pull(Term),
                revolutissima = rev_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term),
                yahouensis = yah_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term),
                impolita = imp_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term),
                vieillardii=vie_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term))

upset(fromList(listInput), order.by = "freq", nsets = 6, shade.alpha = 0.25)

pdf("upset_goterms.pdf", height=10, width=17)
ComplexUpset::upset(fromList(listInput), 
                    c("sandwicensis", "pancheri", "revolutissima", "yahouensis", "impolita", "vieillardii"),
                    base_annotations=list('Intersection size'=intersection_size(counts=FALSE)),
                    themes=upset_default_themes(text=element_text(size=30)),
                    queries=list(
                      upset_query(set='sandwicensis', fill='dodgerblue'),
                      upset_query(set='pancheri', fill='orchid3'),
                      upset_query(set='revolutissima', fill='limegreen'),
                      upset_query(set='yahouensis', fill='orchid3'),
                      upset_query(set='impolita', fill='orchid3'),
                      upset_query(set='vieillardii', fill='limegreen')))
dev.off()




Reduce(intersect, list(rev_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term),
                       yah_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term),
                       imp_0_GO$result  %>% filter(classicFisher < 0.05) %>% pull(Term)))


#sand_0_GO0000373<-get_de_genes_in_term(sand_0$transcript, "GO:0000373", sand_0_GO$goData)
#vie_0_GO0000373<-get_de_genes_in_term(vie_0$transcript, "GO:0000373", vie_0_GO$goData)
#panc_0_GO0000373<-get_de_genes_in_term(panc_0$transcript, "GO:0000373", panc_0_GO$goData)
#rev_0_GO0000373<-get_de_genes_in_term(rev_0$transcript, "GO:0000373", rev_0_GO$goData)
#yah_0_GO0000373<-get_de_genes_in_term(yah_0$transcript, "GO:0000373", yah_0_GO$goData)
#imp_0_GO0000373<-get_de_genes_in_term(imp_0$transcript, "GO:0000373", imp_0_GO$goData)

#sand_0 %>% filter(transcript %in% sand_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#vie_0 %>% filter(transcript %in% vie_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#panc_0 %>% filter(transcript %in% panc_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#rev_0 %>% filter(transcript %in% rev_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#yah_0 %>% filter(transcript %in% yah_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#imp_0 %>% filter(transcript %in% imp_0_GO0000373) %>% dplyr::select(transcript, term, ins) %>% pull(ins) %>% as.numeric() %>% summary()
#
#sand_0_GO0000373<-get_de_genes_in_term(sand_0$transcript, "GO:0000373", sand_0_GO$goData) %>% paste("sandwicensis", ., sep="_")
#vie_0_GO0000373<-get_de_genes_in_term(vie_0$transcript, "GO:0000373", vie_0_GO$goData) %>% paste("vieillardii", ., sep="_")
#panc_0_GO0000373<-get_de_genes_in_term(panc_0$transcript, "GO:0000373", panc_0_GO$goData) %>% paste("pancheri", ., sep="_")
#rev_0_GO0000373<-get_de_genes_in_term(rev_0$transcript, "GO:0000373", rev_0_GO$goData) %>% paste("revolutissima", ., sep="_")
#yah_0_GO0000373<-get_de_genes_in_term(yah_0$transcript, "GO:0000373", yah_0_GO$goData) %>% paste("yahouensis", ., sep="_")
#imp_0_GO0000373<-get_de_genes_in_term(imp_0$transcript, "GO:0000373", imp_0_GO$goData) %>% paste("impolita", ., sep="_")
#
#write.table(c(sand_0_GO0000373, vie_0_GO0000373, panc_0_GO0000373, 
#              rev_0_GO0000373, yah_0_GO0000373, imp_0_GO0000373), 
#            file="testing.txt", 
#            quote=FALSE, 
#            row.names = FALSE, 
#            col.names = FALSE)
#

###################################################################################################
#.     get genes with a TE in them and an equal number of rand picked genes without TE in them.   #
#.     tpo run each on orthofinder and ask if there are more duplicated genes in TE pile.         #
###################################################################################################


# get genes with a TE at least 1KB away
sand_mt1K<-filter_by_dist(sandwicensis_GO, sandwicensis.gene_te, "dist > 5000")
vie_mt1K<-filter_by_dist(vieillardii_GO, vieillardii.gene_te, "dist > 5000")
panc_mt1K<-filter_by_dist(pancheri_GO, pancheri.gene_te, "dist > 5000")
rev_mt1K<-filter_by_dist(revolutissima_GO, revolutissima.gene_te, "dist > 5000")
yah_mt1K<-filter_by_dist(yahouensis_GO, yahouensis.gene_te, "dist > 5000")
imp_mt1K<-filter_by_dist(impolita_GO, impolita.gene_te, "dist > 5000")

# randomly sample the same number of genes >1kb to next TE as those with a TE 0bp away
sand_mt1K %>% sample_n(sand_0$transcript %>% unique %>% length()) %>% pull(transcript)
vie_mt1K %>% sample_n(vie_0$transcript %>% unique %>% length()) %>% pull(transcript)
panc_mt1K %>% sample_n(panc_0$transcript %>% unique %>% length()) %>% pull(transcript)
rev_mt1K %>% sample_n(rev_0$transcript %>% unique %>% length()) %>% pull(transcript)
yah_mt1K %>% sample_n(yah_0$transcript %>% unique %>% length()) %>% pull(transcript)
imp_mt1K %>% sample_n(imp_0$transcript %>% unique %>% length()) %>% pull(transcript)






orthogroups<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_family_analysis/fastas/OrthoFinder/Results_Sep09/Orthogroups/Orthogroups.tsv") %>% 
  set_colnames(c("orthogroup", "oleifera", "impolita", "pancheri", "revolutissima", "sandwicensis", "vieillardii", "yahouensis"))
orthocounts<-read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_family_analysis/fastas/OrthoFinder/Results_Sep09/Orthogroups/Orthogroups.GeneCount.tsv") %>% 
  set_colnames(c("orthogroup", "oleifera", "impolita", "pancheri", "revolutissima", "sandwicensis", "vieillardii", "yahouensis", "total"))




sand_0_orthogroups<-sapply(sand_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, sandwicensis)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
vie_0_orthogroups<-sapply(vie_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, vieillardii)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
panc_0_orthogroups<-sapply(panc_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, pancheri)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
rev_0_orthogroups<-sapply(rev_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, revolutissima)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
yah_0_orthogroups<-sapply(yah_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, yahouensis)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
imp_0_orthogroups<-sapply(imp_0 %>% pull(transcript) %>% unique(), function(x) orthogroups %>% dplyr::filter(grepl(x, impolita)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()

sig_cafe<-read.table("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_family_analysis/cafe/results/Base_family_results_significant.txt")
all_0_orthogroups<-c(sand_0_orthogroups, vie_0_orthogroups, panc_0_orthogroups, rev_0_orthogroups, yah_0_orthogroups, imp_0_orthogroups)
intersect(sig_cafe$V1, all_0_orthogroups)

sand_mt1K_orthogroups<-sapply(sand_mt1K %>% sample_n(sand_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, sandwicensis)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
vie_mt1K_orthogroups<-sapply(vie_mt1K %>% sample_n(vie_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, vieillardii)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
panc_mt1K_orthogroups<-sapply(panc_mt1K %>% sample_n(panc_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, pancheri)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
rev_mt1K_orthogroups<-sapply(rev_mt1K %>% sample_n(rev_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, revolutissima)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
yah_mt1K_orthogroups<-sapply(yah_mt1K %>% sample_n(yah_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, yahouensis)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()
imp_mt1K_orthogroups<-sapply(imp_mt1K %>% sample_n(imp_0$transcript %>% unique %>% length()) %>% pull(transcript), function(x) orthogroups %>% dplyr::filter(grepl(x, impolita)) %>% pull(orthogroup), simplify = TRUE) %>% unlist() %>% as.vector()

all_mt1K_orthogroups<-c(sand_mt1K_orthogroups, vie_mt1K_orthogroups, panc_mt1K_orthogroups, 
  rev_mt1K_orthogroups, yah_mt1K_orthogroups, imp_mt1K_orthogroups)

intersect(sig_cafe$V1, all_mt1K_orthogroups)

sig_cafe$V1

intersect(sig_cafe$V1, sand_0_orthogroups) %>% length()
intersect(sig_cafe$V1, vie_0_orthogroups) %>% length()
intersect(sig_cafe$V1, panc_0_orthogroups) %>% length()
intersect(sig_cafe$V1, rev_0_orthogroups) %>% length()
intersect(sig_cafe$V1, yah_0_orthogroups) %>% length()
intersect(sig_cafe$V1, imp_0_orthogroups) %>% length()


intersect(sig_cafe$V1, sand_mt1K_orthogroups) %>% length()
intersect(sig_cafe$V1, vie_mt1K_orthogroups) %>% length()
intersect(sig_cafe$V1, panc_mt1K_orthogroups) %>% length()
intersect(sig_cafe$V1, rev_mt1K_orthogroups) %>% length()
intersect(sig_cafe$V1, yah_mt1K_orthogroups) %>% length()
intersect(sig_cafe$V1, imp_mt1K_orthogroups) %>% length()


test<-rbind(data.frame(count=orthocounts %>% filter(orthogroup %in% sand_0_orthogroups) %>% pull(sandwicensis),
           species="sandwicensis",
           tedist="0"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% sand_mt1K_orthogroups) %>% pull(sandwicensis),
                 species="sandwicensis",
                 tedist="mt1k"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% vie_0_orthogroups) %>% pull(vieillardii),
                 species="vieillardii",
                 tedist="0"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% vie_mt1K_orthogroups) %>% pull(vieillardii),
                 species="vieillardii",
                 tedist="mt1k"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% panc_0_orthogroups) %>% pull(pancheri),
                 species="pancheri",
                 tedist="0"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% panc_mt1K_orthogroups) %>% pull(pancheri),
                 species="pancheri",
                 tedist="mt1k"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% rev_0_orthogroups) %>% pull(revolutissima),
                 species="revolutissima",
                 tedist="0"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% rev_mt1K_orthogroups) %>% pull(revolutissima),
                 species="revolutissima",
                 tedist="mt1k"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% yah_0_orthogroups) %>% pull(yahouensis),
                 species="yahouensis",
                 tedist="0"),
      data.frame(count=orthocounts %>% filter(orthogroup %in% yah_mt1K_orthogroups) %>% pull(yahouensis),
                 species="yahouensis",
                 tedist="mt1k"))

test %>%
  filter(count < 20 & count > 0) %>%
  ggplot(aes(x=species, y=(count), fill=tedist)) + 
  #geom_boxplot() +
  #geom_jitter(color="black", size=0.4, alpha=0.9)
  geom_violin()

test %>% filter(species == "sandwicensis") %>% group_by(tedist) %>% summarise(count=mean(count))
test %>% filter(species == "pancheri") %>% group_by(tedist) %>% summarise(count=mean(count))
test %>% filter(species == "revolutissima") %>% group_by(tedist) %>% summarise(count=mean(count))



test %>% group_by(species, tedist) %>% summarise(count=median(count))





sand_mt1K %>% sample_n(sand_0$transcript %>% unique %>% length()) %>% pull(transcript)
vie_mt1K %>% sample_n(vie_0$transcript %>% unique %>% length()) %>% pull(transcript)
panc_mt1K %>% sample_n(panc_0$transcript %>% unique %>% length()) %>% pull(transcript)
rev_mt1K %>% sample_n(rev_0$transcript %>% unique %>% length()) %>% pull(transcript)
yah_mt1K %>% sample_n(yah_0$transcript %>% unique %>% length()) %>% pull(transcript)
imp_mt1K %>% sample_n(imp_0$transcript %>% unique %>% length()) %>% pull(transcript)


write.table(sand_mt1K, "sand_mt1K",
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


write.table(c(sand_0_GO0000373, vie_0_GO0000373, panc_0_GO0000373, 
              rev_0_GO0000373, yah_0_GO0000373, imp_0_GO0000373), 
            file="testing.txt", 
            quote=FALSE, 
            row.names = FALSE, 
            col.names = FALSE)









sand_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)
vie_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)
panc_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)
rev_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)
yah_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)
imp_0_GO$result %>% filter(as.numeric(classicFisher) < 0.005)


test_all<-rbind(sand_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "sandwicensis"),
            vie_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "viellardii"),
            imp_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "impolita"),
            rev_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "revolutissima"),
            yah_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "yahouensis"),
            panc_0$result %>% filter(as.numeric(classicFisher) < 0.5) %>% dplyr::select(GO.ID, Term, Significant) %>% mutate(species = "pancheri"))


############################################################################################################
#         get table of GO term and number of genes with a TE in them annotated to that term per species.   #
#############################################################################################################

gene_te_0_table<-rbind(sand_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "sandwicensis"),
vie_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "viellardii"),
imp_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "impolita"),
rev_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "revolutissima"),
yah_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "yahouensis"),
panc_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "pancheri")) %>% 
  reshape(idvar = "Term", timevar = "species", direction = "wide") %>% 
  replace_na(list(Significant.sandwicensis = 0,
                  Significant.viellardii = 0, 
                  Significant.impolita = 0, 
                  Significant.revolutissima = 0, 
                  Significant.yahouensis = 0, 
                  Significant.pancheri = 0))

gene_te_0_table_nonsig<-rbind(sand_0$result %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Significant) %>% mutate(species = "sandwicensis"),
                       vie_0$result %>% filter(as.numeric(classicFisher) > 0.1) %>% dplyr::select(Term, Significant) %>% mutate(species = "viellardii"),
                       imp_0$result %>% filter(as.numeric(classicFisher) > 0.1) %>% dplyr::select(Term, Significant) %>% mutate(species = "impolita"),
                       rev_0$result %>% filter(as.numeric(classicFisher) > 0.1) %>% dplyr::select(Term, Significant) %>% mutate(species = "revolutissima"),
                       yah_0$result %>% filter(as.numeric(classicFisher) > 0.1) %>% dplyr::select(Term, Significant) %>% mutate(species = "yahouensis"),
                       panc_0$result %>% filter(as.numeric(classicFisher) > 0.1) %>% dplyr::select(Term, Significant) %>% mutate(species = "pancheri")) %>% 
  reshape(idvar = "Term", timevar = "species", direction = "wide") %>% 
  replace_na(list(Significant.sandwicensis = 0,
                  Significant.viellardii = 0, 
                  Significant.impolita = 0, 
                  Significant.revolutissima = 0, 
                  Significant.yahouensis = 0, 
                  Significant.pancheri = 0))


# get terms with at least N genes annotated to them per species
gene_te_0_table[rowSums(gene_te_0_table>1) >=5,] %>% rownames_to_column(var="remove") %>% dplyr::select(-"remove") %>% column_to_rownames(var="Term") %>% pheatmap::pheatmap(cluster_rows=F, cluster_cols=F) 
 
##########################################################################
#.  HEATMAP minimum 10 distance between ultramafic and volcnic species.  #
##########################################################################

gene_te_0_table %>% rownames_to_column(var="remove") %>% 
  dplyr::select(-"remove") %>% 
  column_to_rownames(var="Term") %>% 
  set_colnames(c("sandwicensis", "viellardii", "impolita", "revolutissima", "yahouensis", "pancheri")) %>%
  filter(revolutissima - impolita > 10 & 
         pancheri - yahouensis > 10) %>% 
  pheatmap::pheatmap(cluster_rows=F, cluster_cols=F, fontsize=20, angle_col=315)


#############################################################
#.  HEATMAP mean normalised hatmap mean normalised by row.  #
#############################################################

test$row_mean<-rowMeans(test)
test %>% mutate(Significant.sandwicensis=Significant.sandwicensis/row_mean,
                Significant.viellardii =Significant.viellardii /row_mean,
                Significant.impolita=Significant.impolita/row_mean,
                Significant.revolutissima=Significant.revolutissima/row_mean,
                Significant.yahouensis=Significant.yahouensis/row_mean,
                Significant.pancheri=Significant.pancheri/row_mean) %>%
  dplyr::select(-"row_mean") %>%
  set_colnames(c("sandwicensis", "viellardii", "impolita", "revolutissima", "yahouensis", "pancheri")) %>%
  pheatmap::pheatmap(cluster_rows=F, cluster_cols=F)



#############################################
#.  get genes in term for DNA metabolism .  #
#############################################


dna_met_revo<-get_de_genes_in_term(inner_join(revolutissima_GO, revolutissima.gene_te, by="gene") %>% 
                                     arrange(gene) %>% dplyr::select(-none) %>% 
                                     filter(abs(dist) == 0 & proc == "biological_process") %>%
                                     pull(transcript), "GO:0006259", rev_0$goData)

dna_met_impo<-get_de_genes_in_term(inner_join(impolita_GO, impolita.gene_te, by="gene") %>% 
                                     arrange(gene) %>% dplyr::select(-none) %>% 
                                     filter(abs(dist) == 0 & proc == "biological_process") %>%
                                     pull(transcript), "GO:0006259", imp_0$goData)

dna_met_panc<-get_de_genes_in_term(inner_join(pancheri_GO, pancheri.gene_te, by="gene") %>% 
                                     arrange(gene) %>% dplyr::select(-none) %>% 
                                     filter(abs(dist) == 0 & proc == "biological_process") %>%
                                     pull(transcript), "GO:0006259", panc_0$goData)

dna_met_yaho<-get_de_genes_in_term(inner_join(yahouensis_GO, yahouensis.gene_te, by="gene") %>% 
                                     arrange(gene) %>% dplyr::select(-none) %>% 
                                     filter(abs(dist) == 0 & proc == "biological_process") %>%
                                     pull(transcript), "GO:0006259", yah_0$goData)

dna_met_vie<-get_de_genes_in_term(inner_join(vieillardii_GO, vieillardii.gene_te, by="gene") %>% 
                                     arrange(gene) %>% dplyr::select(-none) %>% 
                                     filter(abs(dist) == 0 & proc == "biological_process") %>%
                                     pull(transcript), 
                                  "GO:0006259", 
                                  vie_0$goData)



inner_join(impolita_GO, impolita.gene_te, by="gene") %>% filter(transcript %in% dna_met_impo & proc == "biological_process") %>% data.frame()
inner_join(revolutissima_GO, revolutissima.gene_te, by="gene") %>% filter(transcript %in% dna_met_revo& proc == "biological_process") %>% data.frame()
inner_join(yahouensis_GO, yahouensis.gene_te, by="gene") %>% filter(transcript %in% dna_met_yaho& proc == "biological_process") %>% data.frame()
inner_join(pancheri_GO, pancheri.gene_te, by="gene") %>% filter(transcript %in% dna_met_panc& proc == "biological_process") %>% data.frame()

toplot<-rbind(inner_join(impolita_GO, impolita.gene_te, by="gene") %>% filter(transcript %in% dna_met_impo & proc == "biological_process") %>% dplyr::select(ins) %>% mutate(species="impolita"),
inner_join(revolutissima_GO, revolutissima.gene_te, by="gene") %>% filter(transcript %in% dna_met_revo& proc == "biological_process") %>% dplyr::select(ins) %>% mutate(species="revolutissima"),
inner_join(yahouensis_GO, yahouensis.gene_te, by="gene") %>% filter(transcript %in% dna_met_yaho& proc == "biological_process") %>% dplyr::select(ins) %>% mutate(species="yahounensis"),
inner_join(pancheri_GO, pancheri.gene_te, by="gene") %>% filter(transcript %in% dna_met_panc& proc == "biological_process") %>% dplyr::select(ins) %>% mutate(species="pancheri"))

toplot %>% ggplot(aes(x=as.numeric(ins), fill=species)) + geom_histogram() + facet_wrap(~species)

toplot$species <-factor(toplot$species, levels=c("yahounensis", "pancheri", "impolita", "revolutissima"))

png("te_gene_enrichment_date", width=800, height=600)
toplot %>% 
  ggplot(aes(x=as.numeric(ins), fill=species)) + 
  geom_density(aes(y = after_stat(count)), alpha = 0.25) + 
  facet_wrap(~species) +
  theme(legend.text = element_text(size=22),
        legend.position = "none",
        strip.text = element_text(size = 20, face="italic"),
        axis.text.x = element_text(size=13),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=20)) +
  ylab("Number of LTR Retrotrannsposons") +
  xlab("How recently inserted into the genome") +
  scale_fill_manual(values=c("orchid3", "limegreen", "orchid3", "limegreen", "limegreen", "cornflowerblue"))
dev.off()







read_delim("to_r/pancheri_GO_database_result.txt")
read_delim("to_r/pancheri.gene_te_dists_annotation")


ann<-read_delim("to_r/GO_database_result.txt", col_names = c("transcript", "go", "term", "proc")) %>% mutate(gene=str_split_i(transcript, "\\.", 1))
dis<-read_delim("to_r/revolutissima.gene_te_dists", col_names = c("gene", "dist", "te", "ins", "none"))

inner_join(ann, dis, by="gene") %>% 
  arrange(gene) %>% dplyr::select(-none) %>% 
  filter(abs(dist) > 10000 & proc == "biological_process") %>%
  group_by(term) %>%
  summarise(count=n()) %>%
  arrange(desc(count))

revtest<-inner_join(ann, dis, by="gene") %>% 
  arrange(gene) %>% dplyr::select(-none) %>% 
  filter(abs(dist) < 500 & proc == "biological_process") %>%
  pull(transcript) %>%
  get_enriched_terms(., mp_revolutissima, return_sample_GOData=TRUE)

revtest_0<-inner_join(ann, dis, by="gene") %>% 
  arrange(gene) %>% dplyr::select(-none) %>% 
  filter(abs(dist) == 0 & proc == "biological_process") %>%
  pull(transcript) %>%
  get_enriched_terms(., mp_revolutissima, return_sample_GOData=TRUE)

revtest_mt10k<-inner_join(ann, dis, by="gene") %>% 
  arrange(gene) %>% dplyr::select(-none) %>% 
  filter(abs(dist) > 10000 & proc == "biological_process") %>%
  pull(transcript) %>%
  get_enriched_terms(., mp_revolutissima, return_sample_GOData=TRUE)






inner_join(ann, dis, by="gene") %>% 
  arrange(gene) %>% dplyr::select(-none) %>% 
  filter(abs(dist) == 0 & proc == "biological_process" & go %in% c("GO:0000373")) %>%
  pull(transcript)

rt_0_go<-revtest_0$result %>% filter(as.numeric(classicFisher) < 0.00005) %>% pull(GO.ID)
revtest_0$result %>% filter(as.numeric(classicFisher) < 0.0005)
revtest_mt10k$result %>% filter(as.numeric(classicFisher) < 0.0005)



#############################################################
#    read in annotation of genes within 1KB of a TE         #
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
     


test_mean_centred<- test %>% mutate(pancheri=pancheri/sum(pancheri, revolutissima, vieillardii, yahouensis, impolita)/5,
                                    revolutissima=revolutissima/sum(pancheri, revolutissima, vieillardii, yahouensis, impolita)/5,
                                    vieillardii=vieillardii/sum(pancheri, revolutissima, vieillardii, yahouensis, impolita)/5,
                                    yahouensis=yahouensis/sum(pancheri, revolutissima, vieillardii, yahouensis, impolita)/5,
                                    impolita=impolita/sum(pancheri, revolutissima, vieillardii, yahouensis, impolita)/5)

gf<-test_mean_centred %>% filter(stdev > 0.002) %>% pull(GO)

test %>% filter(GO %in% gf)

test_mean_centred$stdev<- test_mean_centred %>% dplyr::select(pancheri, revolutissima, vieillardii, yahouensis, impolita) %>% apply(1, sd)



test_mean_centred$stdev %>% hist()

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



intersect(impolita.gypsyLTR_go$result$GO.ID, revolutissima.gypsyLTR_go$result$GO.ID)
#setdiff(impolita.gypsyLTR_go$result$GO.ID, revolutissima.gypsyLTR_go$result$GO.ID)


#sapply(testing4, function(x) get_de_genes_in_term(impolita.gypsyLTR$X1, x, impolita.gypsyLTR_go$goData) %>% data.frame())



#genes2test<-get_de_genes_in_term(impolita.gypsyLTR$X1, "GO:0016192", impolita.gypsyLTR_go$goData) %>% data.frame() %>% pull() %>% str_split_i("\\.", 1)
genes2test_impolita<-get_de_genes_in_term(impolita.gypsyLTR$X1, "GO:0016192", impolita.gypsyLTR_go$goData) %>% data.frame() %>% pull() %>% str_split_i("\\.", 1)
genes2test_revolutissima<-get_de_genes_in_term(revolutissima.gypsyLTR$X1, "GO:0016192", revolutissima.gypsyLTR_go$goData) %>% data.frame() %>% pull() %>% str_split_i("\\.", 1)


impolita.gene_te_dists<-read.table("impolita.gene_te_dists") %>% set_colnames(c("geneid", "basepairs", "classification", "insertion"))
revolutissima.gene_te_dists<-read.table("revolutissima.gene_te_dists") %>% set_colnames(c("geneid", "basepairs", "classification", "insertion"))


impolita.gene_te_dists %>% filter(abs(basepairs) == 0 & classification %in% c("LTR/Copia") & geneid %in% genes2test_impolita) %>% pull(geneid)
revolutissima.gene_te_dists %>% filter(abs(basepairs) == 0 & classification %in% c("LTR/Copia") & geneid %in% genes2test_revolutissima) %>% pull(geneid)

genes2test_impolita<-"g1123"
genes2test_revolutissima <-"g4"

#impolita.gene_te_dists %>% filter(geneid == "g1123")
#revolutissima.gene_te_dists %>% filter(geneid == "g4")  

pdf("impolita_dists.pdf", height=6, width=4)
impolita.gene_te_dists %>% 
  mutate(colour_by=ifelse(geneid %in% genes2test_impolita, "yes", "no")) %>% 
  filter(classification %in% c("LTR/Copia", "LTR/Gypsy", "LTR/unknown")) %>% 
  arrange((colour_by)) %>% 
  #filter(abs(basepairs) < 10000 & classification %in% c("LTR/Copia", "LTR/Gypsy")) %>% 
  drop_na() %>%
  ggplot(aes(x=basepairs, y=insertion, colour=colour_by)) +
  geom_point(size=1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~classification) + 
  scale_colour_manual(values = c("grey89", "red"))
dev.off()



pdf("revolutissima_dists.pdf", height=6, width=4)
revolutissima.gene_te_dists %>% 
  mutate(colour_by=ifelse(geneid %in% genes2test_revolutissima, "yes", "no")) %>% 
  filter(classification %in% c("LTR/Copia", "LTR/Gypsy", "LTR/unknown")) %>% 
  arrange((colour_by)) %>% 
  #filter(abs(basepairs) < 10000 & classification %in% c("LTR/Copia", "LTR/Gypsy")) %>% 
  drop_na() %>%
  ggplot(aes(x=basepairs, y=insertion, colour=colour_by)) +
  geom_point(size=1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(~classification) + 
  scale_colour_manual(values = c("grey89", "red"))
dev.off()












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






