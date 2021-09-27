#set your working directory to the dada_project_local directory with your files from Charlie
setwd("/path/to/local/dada_project_local")
taxa <- read.csv("ASV_taxa.csv", header=T, row.names=1)
head(taxa)
seqtab.nochim <- read.csv("seqtab-nochim.csv", header=T, row.names=1)
head(seqtab.nochim)
metadata <- read.csv("info.csv", header=T, row.names=1)
head(metadata)

#install phyloseq if you have do not have this package
# You can install the latest version of phyloseq from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install("phyloseq")
library(phyloseq)
library(ggplot2) #if you do not have ggplot2 type install.packages("ggplot2")

###############################################PHYLOSEQ_START###############################################################
rawotus <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=T), 
                    sample_data(metadata), 
                    tax_table(as.matrix(taxa)))
rawotus@sam_data
head(rawotus@otu_table)
head(rawotus@tax_table)

wh0 <-  genefilter_sample(rawotus, filterfun_sample(function(x) x > 2), A=2)
ps <- prune_taxa(wh0, rawotus)

#make a barplot of the top 200 most abundant ASVs, color coded by the highest taxonomic level (Supergroup)
top200 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:200]
ps.top200 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top200 <- prune_taxa(top200, ps.top200)
plot_bar(ps.top200, x="Name", fill="Supergroup") + 
  theme(text = element_text(size = 18), legend.position = "right")

#dig deeper and compare just Alveolata proportions
ps_alv <- subset_taxa(ps, Supergroup == "Alveolata")
top20 <- names(sort(taxa_sums(ps_alv), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_alv, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Name", fill="Genus") + 
  theme(text = element_text(size = 18), legend.position = "right")

#NMDS - not very meaningful for 3 samples but including code for folks who may want to modify for their dataset
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
p <- plot_ordination(ps.prop, ord.nmds.bray, color="Site", shape="Collection.Month", title="Bray NMDS")+
  geom_point(size = 7) +
  geom_text(mapping = aes(label = Name), size = 4,vjust = 2.5)#+
p

#heatmap to compare ASVs between samples
top25 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:25]
ps.top25 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top25 <- prune_taxa(top25, ps.top25)
plot_heatmap(ps.top25, taxa.label="Genus", sample.label="Name")
