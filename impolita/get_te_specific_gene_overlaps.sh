
#########################################################
#              set species name and window width        #
#########################################################

species="impolita"
window="1000"

#########################################
#             symlink relevant files    #
#########################################

#ln -s /lisc/scratch/botany/katie/annotation/braker3/$species/braker/braker.gtf .
#ln -s /lisc/scratch/botany/katie/annotation/edta_repeatmodeller_output/$species.fasta.mod.EDTA.final/$species.fasta.mod.EDTA.intact.gff3 .


###########################################################
#              extract gtf records for gene only          #
###########################################################

cat braker.gtf | awk '$3 == "gene"' > braker_gene.gtf

#####################################################################################
#   modify braker gtf to reflect the modifications made by EDTA in the gff3 file    #
#####################################################################################
sed -i 's/Scaffolds_//g' braker_gene.gtf

###########################################################
#   extract main TE families into separate GFF3 files     #
###########################################################

grep "helitron" $species.fasta.mod.EDTA.intact.gff3 > $species.fasta.mod.EDTA.intact.helitron.gff3
grep "CACTA_TIR_transposon" $species.fasta.mod.EDTA.intact.gff3 > $species.fasta.mod.EDTA.intact.cactaTIR.gff3
grep "PIF_Harbinger_TIR_transposon" $species.fasta.mod.EDTA.intact.gff3 > $species.fasta.mod.EDTA.intact.harbingerTIR.gff3
grep "Tc1_Mariner_TIR_transposon" $species.fasta.mod.EDTA.intact.gff3  > $species.fasta.mod.EDTA.intact.marinerTIR.gff3
grep "Mutator_TIR_transposon" $species.fasta.mod.EDTA.intact.gff3  > $species.fasta.mod.EDTA.intact.mutatorTIR.gff3
grep "Copia_LTR_retrotransposon" $species.fasta.mod.EDTA.intact.gff3 > $species.fasta.mod.EDTA.intact.copiaLTR.gff3
grep "Gypsy_LTR_retrotransposon" $species.fasta.mod.EDTA.intact.gff3 > $species.fasta.mod.EDTA.intact.gypsyLTR.gff3

###################################################################
#     use bedtools window to get TEs in 10KB window around genes  #
###################################################################

bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.helitron.gff3 > $species.gene_helitron_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.cactaTIR.gff3 > $species.gene_cactaTIR_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.harbingerTIR.gff3 > $species.gene_harbingerTIR_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.marinerTIR.gff3 > $species.gene_marinerTIR_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.mutatorTIR.gff3 > $species.gene_mutatorTIR_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.copiaLTR.gff3 > $species.gene_copiaLTR_window$window.gff3
bedtools window -w $window -a braker_gene.gtf -b $species.fasta.mod.EDTA.intact.gypsyLTR.gff3 > $species.gene_gypsyLTR_window$window.gff3

###################################################################
#                   extract genes with TE within 5KB              #
###################################################################

cut -f 9 $species.gene_helitron_window$window.gff3 > $species.gene_helitron_window$window.genes
cut -f 9 $species.gene_cactaTIR_window$window.gff3 > $species.gene_cactaTIR_window$window.genes
cut -f 9 $species.gene_harbingerTIR_window$window.gff3 > $species.gene_harbingerTIR_window$window.genes
cut -f 9 $species.gene_marinerTIR_window$window.gff3 > $species.gene_marinerTIR_window$window.genes
cut -f 9 $species.gene_mutatorTIR_window$window.gff3 > $species.gene_mutatorTIR_window$window.genes
cut -f 9 $species.gene_copiaLTR_window$window.gff3 > $species.gene_copiaLTR_window$window.genes
cut -f 9 $species.gene_gypsyLTR_window$window.gff3 > $species.gene_gypsyLTR_window$window.genes


###################################################################
#                  get annotation of these genes                  #
###################################################################

annotation_file="/lisc/scratch/botany/katie/annotation/braker3/$species/gfap_annotation/GO_database_result.txt"

while read i; do grep $i $annotation_file; done < $species.gene_helitron_window$window.genes | grep "biological_process" > $species.gene_helitron_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_cactaTIR_window$window.genes | grep "biological_process" > $species.gene_cactaTIR_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_harbingerTIR_window$window.genes | grep "biological_process" > $species.gene_harbingerTIR_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_marinerTIR_window$window.genes | grep "biological_process" > $species.gene_marinerTIR_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_mutatorTIR_window$window.genes | grep "biological_process" > $species.gene_mutatorTIR_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_copiaLTR_window$window.genes | grep "biological_process" > $species.gene_copiaLTR_window$window.annotation
while read i; do grep $i $annotation_file; done < $species.gene_gypsyLTR_window$window.genes | grep "biological_process" > $species.gene_gypsyLTR_window$window.annotation

########################################
#             cleanup                  #
########################################


rm *window*genes *window*gff3
