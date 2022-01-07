library(ade4)
library(genoPlotR)

## Saving data
pdfPath <- ". ./pdfs"

##Mauve version
#bbone <- read_mauve_backbone(bbone_file)
#names <- c("SAP40", "RN6390")
#names(bbone$dna_segs) <- names
#pdf("barto_seg1.pdf")
#plot_gene_map(bbone$dna_segs, bbone$comparisons)
#dev.off()


##Blast version
#In here change genome/comparison names and files
SAP40 <- read_dna_seg_from_genbank('Titus.gbk')
RN6390 <- read_dna_seg_from_genbank('Mak.gb')
comp <- read_comparison_from_blast('WAPA2F13114-Alignment.txt')
pdf("Titusvsmak.pdf",h=4, w=7)  #edit h for more distance between the genes
plot_gene_map(dna_segs=list(SAP40, RN6390),comparisons=list(comp),
              main="Titus vs Mak",
              gene_type="arrows",
              dna_seg_scale=TRUE, scale=FALSE,
              global_color_scheme = c("per_id", "increasing", "blue_red", 1))
dev.off()



