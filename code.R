library(dada2); packageVersion("dada2")
set.seed(100)

#saving plots. width and height variables storing mm length of A4 landscape page to be used throughout
width = 297
height = 210

path <- "/data1/ankit/16s/run1" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

#trimming was done to remove 2 nucleotide in R1 and 11 in R2 reads. Further, primer sequences were removed by triming 17 in R1 and 21 in R2 reads.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE, maxEE = 2, truncLen = c(249,240), trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

#learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#denoise data
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

#merge pairs and extract sequences with merged length in range 400-430
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:430]
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

#remove bimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)

#track reads through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats.txt", sep = "\t", quote = FALSE)

saveRDS(seqtab.nochim, "16sall.rds")
save.image("16s.RData")


# Assign taxonomy using different database
taxSilva <- assignTaxonomy(seqtab.nochim, "/data1/ankit/16s/DADA2-db/silva_nr_v132_train_set.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
saveRDS(taxSilva, "/data1/ankit/16s/taxSilva.rds") 
taxRdp <- assignTaxonomy(seqtab.nochim, "/data1/ankit/16s/DADA2-db/rdp_train_set_16.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
saveRDS(taxRdp, "/data1/ankit/16s/taxRdp.rds") 
taxRefseqRdp <- assignTaxonomy(seqtab.nochim, "/data1/ankit/16s/DADA2-db/RefSeq-RDP16S_v3_May2018.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
saveRDS(taxRefseqRdp, "/data1/ankit/16s/taxRefseqRdp.rds") 
taxGTDB86 <- assignTaxonomy(seqtab.nochim, "/data1/ankit/16s/DADA2-db/GTDB_bac-arc_ssu_r86.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
saveRDS(taxGTDB86, "/data1/ankit/16s/taxGTDB86.rds") 
taxGTDB89 <- assignTaxonomy(seqtab.nochim, "/data1/ankit/16s/DADA2-db/GTDB_bac120_arc122_ssu_r89.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
saveRDS(taxGTDB89, "/data1/ankit/16s/taxGTDB89.rds") 

#ADD SPECIES infomration using same database
speciesSilva <- addSpecies(taxSilva, "/data1/ankit/16s/DADA2-db/silva_species_assignment_v132.fa.gz", tryRC = TRUE, n = 10000)
saveRDS(speciesSilva, "/data1/ankit/16s/speciesSilva.rds")
speciesRdp <- addSpecies(taxRdp, "/data1/ankit/16s/DADA2-db/rdp_species_assignment_16.fa.gz", tryRC = TRUE, n = 10000)
saveRDS(speciesRdp, "/data1/ankit/16s/speciesRdp.rds")
speciesRefseqRdp <- addSpecies(taxRefseqRdp, "/data1/ankit/16s/DADA2-db/RefSeq-RDP_dada2_assignment_species.fa.gz", tryRC = TRUE, n = 10000)
saveRDS(speciesRefseqRdp, "/data1/ankit/16s/speciesRefseqRdp.rds")
speciesGTDB86 <- addSpecies(taxGTDB86, "/data1/ankit/16s/DADA2-db/GTDB_dada2_assignment_species.fa.gz", tryRC = TRUE, n = 10000)
saveRDS(speciesGTDB86, "/data1/ankit/16s/speciesGTDB86.rds")
speciesGTDB89 <- addSpecies(taxGTDB89, "/data1/ankit/16s/DADA2-db/GTDB_dada2_assignment_species.fa", tryRC = TRUE, n = 10000)
saveRDS(speciesGTDB89, "/data1/ankit/16s/speciesGTDB89.rds")

write.table(taxGTDB86, file = "taxGTDB86.tsv", sep = "\t")
write.table(taxGTDB89, file = "taxGTDB89.tsv", sep = "\t")
write.table(taxRdp, file = "taxRdp.tsv", sep = "\t")
write.table(taxRefseqRdp, file = "taxRefseqRdp.tsv", sep = "\t")
write.table(taxSilva, file = "taxSilva.tsv", sep = "\t")
write.table(speciesGTDB86, file = "speciesGTDB86.tsv", sep = "\t")
write.table(speciesGTDB89, file = "speciesGTDB89.tsv", sep = "\t")
write.table(speciesRdp, file = "speciesRdp.tsv", sep = "\t")
write.table(speciesRefseqRdp, file = "speciesRefseqRdp.tsv", sep = "\t")
write.table(speciesSilva, file = "speciesSilva.tsv", sep = "\t")

#summarise assignments for all database at each taxonomic levels and plot them
Silva <- apply(speciesSilva, 2, function(x) length(which(!is.na(x))))
Rdp <- apply(speciesRdp, 2, function(x) length(which(!is.na(x))))
RefseqRdp <- apply(speciesRefseqRdp, 2, function(x) length(which(!is.na(x))))
GTDB86 <- apply(speciesGTDB86, 2, function(x) length(which(!is.na(x))))
GTDB89 <- apply(speciesGTDB89, 2, function(x) length(which(!is.na(x))))
read.counts <- as.data.frame(rbind(Silva, Rdp, RefseqRdp, GTDB86, GTDB89))
read.counts$Database <- row.names(read.counts)
read.counts <- reshape2::melt(data = read.counts, id = "Database")
dbplot <- ggplot(read.counts,aes(x=Database,y=value,fill=variable)) + geom_bar(stat="identity",position = "identity", alpha = 0.8) + scale_y_log10() + labs(fill = "Taxonomy assignment level", y = "Reads assigned") + annotation_logticks(sides = "lr") + theme(legend.position = "bottom") + theme_classic() + scale_fill_brewer(palette = "Dark2")
ggpubr::ggsave(filename = "database-assignment.pdf", plot = dbplot, device = "pdf", width = height, height = width, units = "mm")
#

###metadata comparison####
#read metadata information file containing information of groups and environmental parameters. Also compare all parameters across groups.
sample.data <- read.delim("sample-data.txt")
rownames(sample.data) = sample.data$Names
sample.data$Stage <- factor(x = sample.data$Stage, levels = c("PreSowing", "Seedling", "PreNodulating", "Nodulating", "Flowering", "Matured", "PostHarvest"), ordered = TRUE)
sample.data$Type <- factor(x = sample.data$Type, levels = c("PreSowing", "Bulk", "Rhizosphere", "PostHarvest"), ordered = TRUE)

all.kruskal <- ggpubr::compare_means(formula = c(pH, EC, OC, P2O5, K2O, S, Mn, Cu, Fe, Zn) ~ Type, data = sample.data, method = "kruskal.test", p.adjust.method = "BH")
BR.wilcox <- ggpubr::compare_means(formula = c(pH, EC, OC, P2O5, K2O, S, Mn, Cu, Fe, Zn) ~ Type, data = sample.data[sample.data$Type %in% c("Bulk", "Rhizosphere"), ], method = "wilcox.test", p.adjust.method = "BH", paired = TRUE)
rhizosphere.kruskal <- ggpubr::compare_means(formula = c(pH, EC, OC, P2O5, K2O, S, Mn, Cu, Fe, Zn) ~ Stage, data = sample.data[sample.data$Type == "Rhizosphere",], method = "kruskal.test", p.adjust.method = "BH")
bulk.kruskal <- ggpubr::compare_means(formula = c(pH, EC, OC, P2O5, K2O, S, Mn, Cu, Fe, Zn) ~ Stage, data = sample.data[sample.data$Type == "Bulk",], method = "kruskal.test", p.adjust.method = "BH")

write.table(x = all.kruskal, file = "Type-comparison-metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = BR.wilcox, file = "BR-comparison-metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = rhizosphere.kruskal, file = "Rhizosphere-stage-comparison-metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = bulk.kruskal, file = "Bulk-stage-comparison-metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#


library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(data.table)
library(randomcoloR)
library(tidyr)
library(scales)
library(vegan)
library(ggConvexHull)
library(grid)
library(RColorBrewer)
library(ggnewscale)
library(dplyr)
library(plyr)
library(pairwiseAdonis)
library(DESeq2)
library(ggrepel)
library(corrplot)
library(ape)
library(gridExtra)
library(stringr)
library(UpSetR)
library(viridis)
library(ComplexHeatmap)


#create phyloseq object with metadata information and taxonomy from GTDB
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sample.data), tax_table(speciesGTDB89))

#renaming sequence as ASV names to something more comfortable. "ASV" as prefix followed by numbers.
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
saveRDS(object = ps, file = "ps.RDS")

#to check distribution of reads per sample and counts per ASVs. This will be used to determine further filters.
sort(sample_sums(ps))
table(taxa_sums(ps))
ls <- seq(0, 8500, by = 5000)
#c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000)
hist(sample_sums(ps), main = "Histogram: Read Counts", xlab = "Total Reads", border = "blue", col = c("blue", "blue", rep("green", 15)), breaks = ls, freq = TRUE, ylab = "Frequency of samples")


#removing ASVs with toal count less than 6
ps1 <- prune_taxa(taxa_sums(ps) > 5, ps)
ps1
#removing samples with read counts less than 10000. This will be the final object analysed throughout.
ps2 <- prune_samples(sample_sums(ps1) > 10000, ps1)
ps2
saveRDS(object = ps2, file = "ps2.RDS")

#create a object with relative abundance and subset with Rhizosphere samples only. This is done for simplicity when required.
ps2.rel <- microbiome::transform(ps2, "compositional")
Rhizosphereps2 <- subset_samples(ps2, Type == "Rhizosphere")
Rhizosphereps2.rel <- microbiome::transform(Rhizosphereps2, "compositional")



#########ALPHA DIVERSITY
#create a metadata data frame to add diversity values.
metadata <- data.frame(sample_data(ps2))
#calculating diversity
div <- microbiome::alpha(x = ps2, index = "all")
metadata$ShannonDiversity = div$diversity_shannon
metadata$ObservedASV = div$observed
#storing list to be used in plots
comparisonsList = list( c("Presowing", "Rhizosphere1"), c("Rhizosphere1","Rhizosphere2"), c("Rhizosphere2", "Rhizosphere3"), c("Rhizosphere3", "Rhizosphere4"), c("Rhizosphere4", "Rhizosphere5"), c("Rhizosphere5","Rhizosphere6"), c("Rhizosphere6","Rhizosphere7"), c("Rhizosphere7", "Rhizosphere8"), c("Rhizosphere8","Postharvest"), c("Presowing","Bulk1"), c("Bulk1","Bulk2"), c("Bulk2", "Bulk3"), c("Bulk3", "Bulk4"), c("Bulk4", "Bulk5"), c("Bulk5", "Bulk6"), c("Bulk6", "Bulk7"), c("Bulk7","Bulk8"), c("Bulk8", "Postharvest"), c("Bulk1", "Rhizosphere1"), c("Bulk2", "Rhizosphere2"), c("Bulk3", "Rhizosphere3"), c("Bulk4", "Rhizosphere4"), c("Bulk5", "Rhizosphere5"), c("Bulk6", "Rhizosphere6"), c("Bulk7", "Rhizosphere7"), c("Bulk8", "Rhizosphere8")) 
lable_order <- c("Presowing", "Bulk1", "Rhizosphere1", "Bulk2", "Rhizosphere2", "Bulk3", "Rhizosphere3", "Bulk4", "Rhizosphere4", "Bulk5", "Rhizosphere5", "Bulk6", "Rhizosphere6", "Bulk7", "Rhizosphere7", "Bulk8", "Rhizosphere8", "Postharvest" )
ObservedOrder = c(3100,3100,3100,3100,3100,3100,3100,3100,3100, 3400,3400,3400,3400,3400,3400,3400,3400,3400,3700,3700,3700,3700,3700,3700,3700,3700)
ShannonOrder = c(7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7, 8.2,8.2,8.2,8.2,8.2,8.2,8.2,8.2,8.2, 8.7,8.7,8.7,8.7,8.7,8.7,8.7,8.7)
#plot Observed ASVs and Shannon diversity measures
alphaObserved <- ggplot(data = metadata, aes(x = factor(Description, levels = lable_order), y = ObservedASV, color = factor(Type, levels = c("PreSowing", "Bulk", "Rhizosphere", "PostHarvest"))))  + geom_boxplot() + geom_point(size = 3 )   + stat_compare_means(comparisons = comparisonsList, label = "p.format", label.y = ObservedOrder) + theme(legend.position = "top") + ylab("Observed ASVs") + ggtitle("Alpha diversity measures") + theme(plot.title = element_text(hjust = 0.5)) + labs(color = "Sample type")  + theme(axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(plot.margin = unit(c(5,5,1,5), "pt")) + geom_text(mapping = aes(x = 10, y = 4100,  fontface = "plain", size = 3), color = "black", label = paste0("Kruskal-Wallis (plotted groups) p-value: ", compare_means(formula = ObservedASV ~ Description, data = metadata, method = "kruskal.test")$p.format , "          Kruskal-Wallis (R1 to R8) p-value: ", compare_means(formula = ObservedASV ~ Description, data = metadata[metadata$Type == "Rhizosphere",], method = "kruskal.test")$p.format , "          Kruskal-Wallis (B1 to B8) p-value: ", compare_means(formula = ObservedASV ~ Description, data = metadata[metadata$Type == "Bulk",], method = "kruskal.test")$p.format ), show.legend = FALSE, inherit.aes = FALSE) + geom_text(mapping = aes(x = 10, y = 4300, fontface = "plain", size = 3 ), color = "black", label = paste0("Kruskal-Wallis (among types: presowing, bulk, rhizosphere and postharvest soil) p-value: ", compare_means(formula = ObservedASV ~ Type, data = metadata, method = "kruskal.test")$p.format), show.legend = FALSE)
alphaShannon <- ggplot(data = metadata, aes(x = factor(Description, levels = lable_order), y = ShannonDiversity, color = factor(Type, levels = c("PreSowing", "Bulk", "Rhizosphere", "PostHarvest")))) + geom_boxplot() + geom_point(size = 3 )    + stat_compare_means(comparisons = comparisonsList, label = "p.format", label.y = ShannonOrder) + theme(legend.position = "none") + xlab("Sample groups") + ylab("Shannon diversity") + theme(strip.background.x = element_blank(), strip.text.x = element_blank() ) + theme(plot.margin = unit(c(1,5,5,21), "pt")) + geom_text(mapping = aes(x = 10, y = 9.2, fontface = "plain", size = 3 ), color = "black", label = paste0("Kruskal-Wallis (plotted groups) p-value: ", compare_means(formula = ShannonDiversity ~ Description, data = metadata, method = "kruskal.test")$p.format , "          Kruskal-Wallis (R1 to R8) p-value: ", compare_means(formula = ShannonDiversity ~ Description, data = metadata[metadata$Type == "Rhizosphere",], method = "kruskal.test")$p.format , "          Kruskal-Wallis (B1 to B8) p-value: ", compare_means(formula = ShannonDiversity ~ Description, data = metadata[metadata$Type == "Bulk",], method = "kruskal.test")$p.format )) + geom_text(mapping = aes(x = 10, y = 9.5, fontface = "plain", size = 3 ), color = "black", label = paste0("Kruskal-Wallis (among types: presowing, bulk, rhizosphere and postharvest soil) p-value: ", compare_means(formula = ShannonDiversity ~ Type, data = metadata, method = "kruskal.test")$p.format  ))
palpha <- ggarrange(alphaObserved, alphaShannon, ncol = 1, align = "hv")
ggsave(filename = "alpha-plot.pdf", plot = palpha, device = "pdf", width = width*1.25, height = height*1.25, units = "mm")

########plotting taxonomy
#####PHYLUM
# define the levels to glom 
ps2.Phylum.rel <- tax_glom(physeq = ps2.rel, taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified taxa
taxa_names(ps2.Phylum.rel) <- tax_table(ps2.Phylum.rel)[, 2]
taxa_names(ps2.Phylum.rel)[is.na(taxa_names(ps2.Phylum.rel))] <- "Unclassified Phylum"
#create dataframe
ps2.Phylum.rel.df <- data.table(psmelt(ps2.Phylum.rel))
#group df by taxa and calculate median rel. abundance
ps2.Phylum.rel.df[, median := median(Abundance, na.rm = TRUE), by = "OTU"]
#calculate the phyla-wise group significance
phylum.type <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df, method = "kruskal.test", p.adjust.method = "BH", group.by = "OTU")
phylum.type.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df[ps2.Phylum.rel.df$Description != "Rhizosphere1", ], method = "kruskal.test", p.adjust.method = "BH", group.by = "OTU")
phylum.type.wilcox <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df, method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")
phylum.type.wilcox.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df[ps2.Phylum.rel.df$Description != "Rhizosphere1", ], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")
phylum.BR <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df[ps2.Phylum.rel.df$Type %in% c("Bulk", "Rhizosphere"), ], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")
phylum.BR.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Phylum.rel.df[ps2.Phylum.rel.df$Type %in% c("Bulk", "Rhizosphere") & ps2.Phylum.rel.df$Description != "Rhizosphere1", ], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")

# Change name of taxa less than 1%
ps2.Phylum.rel.df[(median <= 0.01), OTU := "Less Abundant Phyla"]
# Creating df with summarized lesser abundant taxa abundance
ps2.Phylum.rel.df <- ps2.Phylum.rel.df[, sum(Abundance), by = list(OTU,Sample,Type,Description,CollectionTimepoint,Stage)]
colnames(ps2.Phylum.rel.df)[7] <- "Abundance"
#get color codes, store in object and remember to not run that code again or different colors everytime
colcodes.phylum <- distinctColorPalette(length(unique(ps2.Phylum.rel.df$OTU))+13)
#plotting stacked bar chart
Phylum.bar <- ggplot(data=ps2.Phylum.rel.df, aes(x=Sample, y=Abundance*100, fill=OTU)) + geom_bar(position="stack", stat="identity", color = "black", size = 0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Phyla")   + scale_fill_manual(values = colcodes.phylum) + labs(fill = "Phylum") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = sort(ps2.Phylum.rel.df[,unique(Type), by = Sample]$V1))) + coord_cartesian(expand = FALSE ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))
ggsave(filename = "Phylum-bar-plot.pdf", plot = Phylum.bar, device = "pdf", width = width, height = height, units = "mm")
##

#####Genus
# # define the levels to glom 
ps2.Genus.rel <- tax_glom(physeq = ps2.rel, taxrank = "Genus", NArm = FALSE)
#melt to dataframe and convert factor to character
ps2.Genus.rel.df <- data.table(psmelt(ps2.Genus.rel))
ps2.Genus.rel.df$Kingdom <- as.character(ps2.Genus.rel.df$Kingdom)
ps2.Genus.rel.df$Phylum <- as.character(ps2.Genus.rel.df$Phylum)
ps2.Genus.rel.df$Class <- as.character(ps2.Genus.rel.df$Class)
ps2.Genus.rel.df$Order <- as.character(ps2.Genus.rel.df$Order)
ps2.Genus.rel.df$Family <- as.character(ps2.Genus.rel.df$Family)
ps2.Genus.rel.df$Genus <- as.character(ps2.Genus.rel.df$Genus)
#NA entries in Genera were renamed to "highest annotated level"_X
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Phylum <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Kingdom, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Class <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Phylum, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Order <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Class, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Family <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Order, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Genus <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Family, "_X")

genus.type <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df, method = "kruskal.test", p.adjust.method = "BH", group.by = "Genus")
genus.type.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[ps2.Genus.rel.df$Description != "Rhizosphere1", ], method = "kruskal.test", p.adjust.method = "BH", group.by = "Genus")
genus.type.wilcox <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df, method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")
genus.type.wilcox.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[ps2.Genus.rel.df$Description != "Rhizosphere1", ], method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")
genus.BR <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[ps2.Genus.rel.df$Type %in% c("Bulk", "Rhizosphere"), ], method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")
genus.BR.R1 <- ggpubr::compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[ps2.Genus.rel.df$Type %in% c("Bulk", "Rhizosphere") & ps2.Genus.rel.df$Description != "Rhizosphere1", ], method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")
write.table(x = genus.type, file = "genus-type.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.type.R1, file = "genus-type-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.type.wilcox, file = "genus-type-wilcoxon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.type.wilcox.R1, file = "genus-type-wilcoxon-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.BR, file = "genus-BR.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = genus.BR.R1, file = "genus-BR-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

ps2Genus.df = ps2.Genus.rel.df
# #group df by taxa and calculate mean rel. abundance and maximum abundance per genus
ps2.Genus.rel.df[, mean := ((sum(Abundance, na.rm = TRUE))/88), by = "Genus"]
ps2.Genus.rel.df[, maximum := max(Abundance[Abundance > 0], na.rm = TRUE), by = "Genus"]
#The median value was not considered for filtering as the values seems to differ drastically for Rhizosphere and other groups. 
#To get a better picture, mean > 0.01(1%) was considered. Further, considering the highly varying range of each genera, additinal condition with maximum abundance of genera > 0.05(5%) was also considered. 
#Those not meeting these criteria were renamed as Lesser abundant genera
ps2.Genus.rel.df[(mean < 0.01 & maximum < 0.05), Genus := " Less Abundant Genera(n=867)"]
#Creating df with summarized lesser abundant taxa abundance
ps2.Genus.rel.df <- ps2.Genus.rel.df[, sum(Abundance), by = list(Genus,Sample,Type,Description,CollectionTimepoint,Stage)]
colnames(ps2.Genus.rel.df)[7] <- "Abundance"
# get color codes, store in object and remember to not run that code again or different colors everytime
colcodes.genus <- distinctColorPalette(length(unique(ps2.Genus.rel.df$Genus))+7)
#plotting stacked bar chart
Genus.bar <- ggplot(data=ps2.Genus.rel.df, aes(x=Sample, y=Abundance*100, fill=Genus)) + geom_bar(position="stack", stat="identity", color = "black", size=0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Genera")   + scale_fill_manual(values = colcodes.genus) + labs(fill = "Genus") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = sort(ps2.Genus.rel.df[,unique(Type), by = Sample]$V1))) + coord_cartesian(expand = FALSE, ylim = c(0, 85) ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))
ggsave(filename = "Genus-bar-plot.pdf", plot = Genus.bar, device = "pdf",dpi = 300, width = width, height = height, units = "mm")
# #plotting box plot for these taxa.
#Since Enterobacteriaceae_X was absent from all samples of PreSowing, it was causing troubles in proper alignment of dots and boxplot. An arbitary close to 0 value was assigned for the purpose. Also, y-axis labels were manually given to avoid scientific notation; the last label was required to adjust length of labels.
ps2.Genus.rel.df[ps2.Genus.rel.df$Genus == "Enterobacteriaceae_X" & ps2.Genus.rel.df$Sample == "PS1"]$Abundance <- 0.00000001
Genus.box <- ggplot(ps2.Genus.rel.df[Genus != " Less Abundant Genera(n=867)"], aes(x=Genus, y=Abundance*100, fill = factor(Type, levels = c("PreSowing", "Bulk", "Rhizosphere", "PostHarvest")))) + geom_point(position=position_dodge(width=0.75), show.legend = FALSE) + geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75, preserve = "single")) + scale_y_continuous(trans = "log10", labels = c("0", "0.01", "0.1", "1", "10", "100", "101")) + theme(legend.position = "left") + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Genus") + xlab("") + theme(plot.margin=unit(c(1,0.25,0.25,0.25), "cm")) + labs(fill = "Type") + coord_cartesian(ylim = c(0.01, 100))
# #compare means of all pairs using Wilcoxon test and melt the data frame
Genus.means <- compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[Genus != " Less Abundant Genera(n=867)"], group.by = "Genus")
Genus.means$group <- paste(Genus.means$group2, Genus.means$group1, sep = " vs ")
Genus.means <- Genus.means[,c("Genus", "p.signif", "group")] %>% spread(Genus, p.signif)
# #compare group means using Kruskal-Wallis test
Genus.group.means <- compare_means(formula = Abundance ~ Type, data = ps2.Genus.rel.df[Genus != " Less Abundant Genera(n=867)"], group.by = "Genus", method = "kruskal.test")
Genus.group.means <- Genus.group.means[order(Genus.group.means$Genus), c("Genus","p.signif", "method")] %>% spread(Genus, p.signif)
# #merge both comparisons
Genus.group.means <- rbind(Genus.group.means[,2:25], Genus.means[,2:25])
rownames(Genus.group.means) <- c("Kruskal-Wallis group test", Genus.means$group)
# #create table
Genus.text <- ggtexttable(Genus.group.means, rows = rownames(Genus.group.means), cols = rep(x = 888888, times = 24), theme = ttheme(rownames.style = rownames_style(size = 10), colnames.style = colnames_style(size = 12, color = NA, fill = NA, linecolor = NA), tbody.style = tbody_style(size = 15))) + theme(plot.margin=unit(c(-13,0.25,0.25,0.25), "cm"))
Genus.plot2 <- ggarrange(Genus.box, Genus.text, ncol = 1)
ggsave(filename = "Genus-box-plot.pdf", plot = Genus.plot2, device = "pdf", width = width*1.66, height = height*1.66, units = "mm")
 

#SPECIES & ASV
#get a species level ps object and identify top 15 species to plot separately
ps2.Species.rel <- tax_glom(physeq = ps2.rel, taxrank = "Species", NArm = FALSE)
Species.name <- names(sort(taxa_sums(subset_taxa(physeq = ps2.Species.rel, !is.na(Species))), decreasing = TRUE)[1:15])
#separate top 15 species in a new ps object and remove them from original ps object. Both of this will be plotted.
ps2.Species <- subset_taxa(physeq = ps2.Species.rel, rownames(tax_table(ps2.Species.rel)) %in% Species.name)
`%notin%` <- Negate(`%in%`)
ps2.Species.rel <- subset_taxa(physeq = ps2.Species.rel, rownames(tax_table(ps2.Species.rel)) %notin% Species.name)
#rename species name to include genus name as well
ps2.Species@tax_table@.Data[,7] = paste(ps2.Species@tax_table@.Data[,6], ps2.Species@tax_table@.Data[,7], sep = " ")

#plot the heatmap
heatmap.species <- plot_heatmap(physeq = ps2.Species, sample.order = "Description", taxa.label = "Species")
heatmap.species <- heatmap.species + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = sort(metadata$Type)), axis.text.y = element_text(size = 7)) + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), labels[i]))}) + aes(fill = Abundance * 100) + theme(plot.margin = unit(c(1,5,5,5), "pt")) + labs(fill = "Abundance") + scale_fill_gradient(limits=c(0.0005,75) , trans=log_trans(4), na.value = "black", breaks = c(0.001,0.01,0.1,1,10,50)) + theme(legend.position="none")
heatmap.ASV <- plot_heatmap(physeq = ps2.Species.rel, sample.order = "Description", max.label = 1000)
heatmap.ASV <- heatmap.ASV  + aes(fill = Abundance * 100) + labs(fill = "Abundance") + scale_fill_gradient(limits=c(0.0005,75), trans=log_trans(4), na.value = "black", breaks = c(0.001,0.01,0.1,1,10,50)) + theme(legend.position="left") + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(5,5,1,25), "pt"))
heatmapPlot <- ggarrange(heatmap.ASV, heatmap.species, nrow = 2, ncol = 1, heights = c(4,1), align = "hv" )
ggsave(filename = "heatmap-ASV-Species.pdf", plot = heatmapPlot, device = "pdf", width = height*2, height = width*2, units = "mm")



########BETA DIVERSITY

#calculating bray-curtis distance and NMDS ordination
ps2.bray.dist = phyloseq::distance(ps2.rel, method="bray")
ps2.bray.ordination = ordinate(ps2.rel, method="NMDS", distance=ps2.bray.dist)
#perform permanova on calculated distance across Type of samples
ps2.adonis <- vegan::adonis(ps2.bray.dist ~ Type, data=metadata)
#taken from microbiomeseq package (http://www.github.com/umerijaz/microbiomeSeq)
ps2.adn_pvalue <- ps2.adonis[[1]][["Pr(>F)"]][1]
ps2.adn_rsquared <- round(ps2.adonis[[1]][["R2"]][1],3)
#use the bquote function to format adonis results to be annotated on the ordination plot.
ps2.signi_label <- paste(cut(ps2.adn_pvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ".")))
ps2.adn_res_format <- bquote(atop(atop("PERMANOVA", R^2==~.(ps2.adn_rsquared)), atop("p-value="~.(ps2.adn_pvalue)~.(ps2.signi_label), phantom())))
#extract NMDS stress for plotting
ps2.stress.label <- paste("NMDS Stress = ", round(ps2.bray.ordination$stress, 3))
#plot ordination. Shapes are numbers 0-9 for respective collections.
shapes <- 48:57
betaplot.ps2 <- plot_ordination(ps2.rel, ps2.bray.ordination, color = "Type", shape = "CollectionTimepoint") + theme(aspect.ratio = 1) + geom_point(size = 5) + scale_shape_manual(values = shapes) + geom_convexhull(inherit.aes = TRUE, alpha = 0.1, show.legend = FALSE ,aes(fill = "Description")) + ggtitle("NMDS plot using Bray-Curtis distance for relative abundance of all samples") + labs(color = "Type")
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.ps2 <- betaplot.ps2 + annotation_custom(grob = textGrob(label = ps2.stress.label, hjust = 0, gp = gpar(fontsize = 20, fontface = "bold")), xmin = -0.45, xmax = -0.45, ymin = 0.55, ymax = 0.55)
#add PERMANOVA results on the graph....adjust the position based on the actual graph. 
betaplot.ps2 <- betaplot.ps2 + annotation_custom(grob = textGrob(label = ps2.adn_res_format, hjust = 0, gp = gpar(fontsize = 20, fontface = "bold")), xmin = -0.3, xmax = -0.3, ymin = 0.4, ymax = 0.4)
ggsave(filename = "beta-plot-all.pdf", plot = betaplot.ps2, device = "pdf", width = width, height = height, units = "mm")
#plotting stressplot for calculated NMDS ordinates
stressplot(ps2.bray.ordination)


############################Looking through Rhizosphere samples###############
metadata.rhizo <- metadata[metadata$Type == "Rhizosphere",]
########plotting taxonomy
#####PHYLUM
# define the levels to glom 
Rhizo.Phylum.rel <- tax_glom(physeq = Rhizosphereps2.rel, taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified taxa
taxa_names(Rhizo.Phylum.rel) <- tax_table(Rhizo.Phylum.rel)[, 2]
taxa_names(Rhizo.Phylum.rel)[is.na(taxa_names(Rhizo.Phylum.rel))] <- "Unclassified Phylum"
#create dataframe
Rhizo.Phylum.rel.df <- data.table(psmelt(Rhizo.Phylum.rel))
#group df by taxa and calculate mean rel. abundance
Rhizo.Phylum.rel.df[, mean := mean(Abundance, na.rm = TRUE), by = "OTU"]

#group-wise significance
Rhizo.stage <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df[Rhizo.Phylum.rel.df$mean > 0.0001,], method = "kruskal.test", p.adjust.method = "BH", group.by = "OTU")
Rhizo.stage.R1 <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df[Rhizo.Phylum.rel.df$mean > 0.0001 & Rhizo.Phylum.rel.df$Description != "Rhizosphere1",], method = "kruskal.test", p.adjust.method = "BH", group.by = "OTU")
Rhizo.stage.wilcox <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df[Rhizo.Phylum.rel.df$mean > 0.0001,], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")
Rhizo.stage.wilcox.R1 <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df[Rhizo.Phylum.rel.df$mean > 0.0001 & Rhizo.Phylum.rel.df$Description != "Rhizosphere1",], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")

write.table(x = Rhizo.stage, file = "Rhizosphere-stage.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.stage.R1, file = "Rhizosphere-stage-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.stage.wilcox, file = "Rhizosphere-stage-wilcoxon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.stage.wilcox.R1, file = "Rhizosphere-stage-wilcoxon-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Change name of taxa less than 1%, storing the names in new column
Rhizo.Phylum.rel.df[(mean <= 0.01), OTU := "Less Abundant Phyla"]
# Creating df with summarized lesser abundant taxa abundance
Rhizo.Phylum.rel.df <- Rhizo.Phylum.rel.df[, sum(Abundance), by = list(OTU,Sample,Type,Description,CollectionTimepoint,Stage)]
colnames(Rhizo.Phylum.rel.df)[7] <- "Abundance"
#compare means of all pairs using Wilcoxon test and melt the data frame
Rhizo.Phylum.means <-  compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df, group.by = "OTU")
Rhizo.Phylum.means$group <- paste(Rhizo.Phylum.means$group2, Rhizo.Phylum.means$group1, sep = " vs ")
Rhizo.Phylum.means <- Rhizo.Phylum.means[,c("OTU", "p.signif", "group")] %>% spread(OTU, p.signif)
#compare group means using Kruskal-Wallis test
Rhizo.Phylum.group.means <- compare_means(formula = Abundance ~ Stage, data = Rhizo.Phylum.rel.df, group.by = "OTU", method = "kruskal.test")
Rhizo.Phylum.group.means <- Rhizo.Phylum.group.means[order(Rhizo.Phylum.group.means$OTU), c("OTU","p.signif", "method")] %>% spread(OTU, p.signif)
#merge both comparisons
Rhizo.Phylum.group.means <- rbind(Rhizo.Phylum.group.means[,2:13], Rhizo.Phylum.means[,2:13])
rownames(Rhizo.Phylum.group.means) <- c("Kruskal-Wallis group test", Rhizo.Phylum.means$group)
#plot box plot
Rhizo.Phylum.box <- ggplot(Rhizo.Phylum.rel.df, aes(x=OTU, y=Abundance*100, fill = factor(Stage, levels=c("Seedling","PreNodulating","Nodulating","Flowering","Matured")))) + geom_point(position=position_dodge(width=0.75), show.legend = FALSE, size=1) + geom_boxplot(alpha = 0.7) + scale_y_log10() + theme(legend.position = "left") + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Phylum") + xlab("") + theme(plot.margin=unit(c(1,0.25,0.25,0.25), "cm")) + labs(fill = "Stage")
#create table. Column names were added and were hidden to make the column width uniform and adjust with other plot
Rhizo.Phylum.text <- ggtexttable(Rhizo.Phylum.group.means, rows = rownames(Rhizo.Phylum.group.means), cols = rep(x = 8888, times = 12), theme = ttheme(rownames.style = rownames_style(size = 11), colnames.style = colnames_style(size = 20, fill = NA, color = NA, linecolor = NA), tbody.style = tbody_style(size = 15))) + theme(plot.margin=unit(c(-3,0.25,0.25,0.25), "cm"))
Rhizo.Phylum.plot2 <- ggarrange(Rhizo.Phylum.box, Rhizo.Phylum.text, nrow=2)
Rhizo.Phylum.plot2
ggsave(filename = "Rhizo.Phylum-box-plot.pdf", plot = Rhizo.Phylum.plot2, device = "pdf", width = width, height = height, units = "mm", dpi = 300)
#

#####Genus
# # define the levels to glom 
Rhizo.ps2.Genus.rel <- tax_glom(physeq = Rhizosphereps2.rel, taxrank = "Genus", NArm = FALSE)
#melt to dataframe and convert factor to character
Rhizo.ps2.Genus.rel.df <- data.table(psmelt(Rhizo.ps2.Genus.rel))
Rhizo.ps2.Genus.rel.df$Kingdom <- as.character(Rhizo.ps2.Genus.rel.df$Kingdom)
Rhizo.ps2.Genus.rel.df$Phylum <- as.character(Rhizo.ps2.Genus.rel.df$Phylum)
Rhizo.ps2.Genus.rel.df$Class <- as.character(Rhizo.ps2.Genus.rel.df$Class)
Rhizo.ps2.Genus.rel.df$Order <- as.character(Rhizo.ps2.Genus.rel.df$Order)
Rhizo.ps2.Genus.rel.df$Family <- as.character(Rhizo.ps2.Genus.rel.df$Family)
Rhizo.ps2.Genus.rel.df$Genus <- as.character(Rhizo.ps2.Genus.rel.df$Genus)
#NA entries in Genera were renamed to "highest annotated level"_X
Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Phylum),]$Phylum <- paste0(Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Phylum),]$Kingdom, "_X")
Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Class),]$Class <- paste0(Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Class),]$Phylum, "_X")
Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Order),]$Order <- paste0(Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Order),]$Class, "_X")
Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Family),]$Family <- paste0(Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Family),]$Order, "_X")
Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Genus),]$Genus <- paste0(Rhizo.ps2.Genus.rel.df[is.na(Rhizo.ps2.Genus.rel.df$Genus),]$Family, "_X")
#group df by taxa and calculate mean rel. abundance
Rhizo.ps2.Genus.rel.df[, mean := (mean(Abundance, na.rm = TRUE)), by = "Genus"]
#The median value was not considered for filtering as the values seems to differ drastically for Rhizosphere and other groups. 
#Comparison of all genera
#group-wise significance
Rhizo.genus.stage <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Rhizo.ps2.Genus.rel.df$mean > 0.00001,], method = "kruskal.test", p.adjust.method = "BH", group.by = "Genus")
Rhizo.genus.stage.R1 <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Rhizo.ps2.Genus.rel.df$mean > 0.00001 & Rhizo.ps2.Genus.rel.df$Description != "Rhizosphere1",], method = "kruskal.test", p.adjust.method = "BH", group.by = "Genus")
Rhizo.genus.stage.wilcox <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Rhizo.ps2.Genus.rel.df$mean > 0.00001,], method = "wilcox.test", p.adjust.method = "BH", group.by = "Genus")
Rhizo.genus.stage.wilcox.R1 <- ggpubr::compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Rhizo.ps2.Genus.rel.df$mean > 0.00001 & Rhizo.ps2.Genus.rel.df$Description != "Rhizosphere1",], method = "wilcox.test", p.adjust.method = "BH", group.by = "OTU")

write.table(x = Rhizo.genus.stage, file = "Rhizosphere-genus-stage.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.genus.stage.R1, file = "Rhizosphere-genus-stage-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.genus.stage.wilcox, file = "Rhizosphere-genus-stage-wilcoxon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = Rhizo.genus.stage.wilcox.R1, file = "Rhizosphere-genus-stage-wilcoxon-R1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#To get a better picture, mean > 0.01(1%) was considered for plotting.
#Those not meeting these criteria were renamed as Lesser abundant genera
Rhizo.ps2.Genus.rel.df[(mean < 0.01), Genus := "Less Abundant Genera"]
#Creating df with summarized lesser abundant taxa abundance
Rhizo.ps2.Genus.rel.df <- Rhizo.ps2.Genus.rel.df[, sum(Abundance), by = list(Genus,Sample,Type,Description,CollectionTimepoint,Stage)]
colnames(Rhizo.ps2.Genus.rel.df)[7] <- "Abundance"
# #plot box plot
Rhizo.Genus.box <- ggplot(Rhizo.ps2.Genus.rel.df[Genus != "Less Abundant Genera"], aes(x=Genus, y=Abundance*100, fill = factor(Stage, levels=c("Seedling","PreNodulating","Nodulating","Flowering","Matured")))) + geom_point(position=position_dodge(width=0.75), show.legend = FALSE, size = 1) + geom_boxplot(alpha = 0.7) + scale_y_log10() + theme(legend.position = "left", axis.text.x = element_text(size = 7)) + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Genus") + xlab("") + theme(plot.margin=unit(c(1,0.25,0.25,0.5), "cm")) + labs(fill = "Stage")
# #compare means of all pairs using Wilcoxon test and melt the data frame
Rhizo.Genus.means <- compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Genus != "Less Abundant Genera"], group.by = "Genus")
Rhizo.Genus.means$group <- paste(Rhizo.Genus.means$group2, Rhizo.Genus.means$group1, sep = " vs ")
Rhizo.Genus.means <- Rhizo.Genus.means[,c("Genus", "p.signif", "group")] %>% spread(Genus, p.signif)
# #compare group means using Kruskal-Wallis test
Rhizo.Genus.group.means <- compare_means(formula = Abundance ~ Stage, data = Rhizo.ps2.Genus.rel.df[Genus != "Less Abundant Genera"], group.by = "Genus", method = "kruskal.test")
Rhizo.Genus.group.means <- Rhizo.Genus.group.means[order(Rhizo.Genus.group.means$Genus), c("Genus","p.signif", "method")] %>% spread(Genus, p.signif)
# #merge both comparisons
Rhizo.Genus.group.means <- rbind(Rhizo.Genus.group.means[,2:19], Rhizo.Genus.means[,2:19])
rownames(Rhizo.Genus.group.means) <- c("Kruskal-Wallis group test", Rhizo.Genus.means$group)
# #create table
Rhizo.Genus.text <- ggtexttable(Rhizo.Genus.group.means, rows = rownames(Rhizo.Genus.group.means), cols = rep(x = 8888, times = 18), theme = ttheme(rownames.style = rownames_style(size = 11), colnames.style = colnames_style(size = 16, fill = NA, color = NA, linecolor = NA), tbody.style = tbody_style(size = 11))) + theme(plot.margin=unit(c(-6,0.25,0.25,0.25), "cm"))
Rhizo.Genus.plot2 <- ggarrange(Rhizo.Genus.box, Rhizo.Genus.text, nrow = 2)
ggsave(filename = "Rhizo.Genus-box-plot.pdf", plot = Rhizo.Genus.plot2, device = "pdf", width = width*1.2, height = height*1.2, units = "mm", dpi = 300)
####


########BETA DIVERSITY
#calculating bray-curtis distance and NMDS ordination
Rhizo.ps2.bray.dist = phyloseq::distance(Rhizosphereps2.rel, method="bray")
Rhizo.ps2.bray.ordination = ordinate(Rhizosphereps2.rel, method="NMDS", distance=Rhizo.ps2.bray.dist)
#fitting all the environmental variables(physical properties and nutrient composition) with calculated NMDS and identifying significantly associated
ef.Rhizo.ps2 <- envfit(Rhizo.ps2.bray.ordination, metadata.rhizo[,6:17], permu = 999)
ef.Rhizo.ps2.df <- as.data.frame(scores(ef.Rhizo.ps2, display = "vectors"))
ef.Rhizo.ps2.df <- cbind(ef.Rhizo.ps2.df, EnvironFactor = rownames(ef.Rhizo.ps2.df), pvalue = ef.Rhizo.ps2$vectors$pvals, color = "black")
ef.Rhizo.ps2.df$color <- as.character(ef.Rhizo.ps2.df$color)
ef.Rhizo.ps2.df[ef.Rhizo.ps2.df$pvalue < 0.05,]$color <- "red"
#perform permanova on calculated distance across Type of samples
Rhizo.ps2.adonis <- vegan::adonis(Rhizo.ps2.bray.dist ~ Stage, data=metadata.rhizo)
#taken from microbiomeseq package (http://www.github.com/umerijaz/microbiomeSeq)
Rhizo.ps2.adn_pvalue<-Rhizo.ps2.adonis[[1]][["Pr(>F)"]][1]
Rhizo.ps2.adn_rsquared<-round(Rhizo.ps2.adonis[[1]][["R2"]][1],3)
#use the bquote function to format adonis results to be annotated on the ordination plot.
Rhizo.ps2.signi_label <- paste(cut(Rhizo.ps2.adn_pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ".")))
Rhizo.ps2.adn_res_format <- bquote(atop(atop("PERMANOVA",R^2==~.(Rhizo.ps2.adn_rsquared)), atop("p-value="~.(Rhizo.ps2.adn_pvalue)~.(Rhizo.ps2.signi_label), phantom())))
#extract NMDS stress for plotting
Rhizo.ps2.stress.label <- paste("NMDS Stress = ", round(Rhizo.ps2.bray.ordination$stress, 3))
Rhizo.shapes <- 49:56
#plot ordination
Rhizo.betaplot.ps2 <- plot_ordination(Rhizosphereps2.rel, Rhizo.ps2.bray.ordination, color="Stage", shape = "CollectionTimepoint") + theme(aspect.ratio=1) + geom_point(size=5) + scale_shape_manual(values=Rhizo.shapes) + geom_convexhull(inherit.aes = TRUE, alpha=0.1, show.legend = FALSE ,aes(fill = "Description")) + ggtitle("NMDS plot using Bray-Curtis distance for relative abundance of Rhizosphere samples")
#add NMDS stress on the graph....adjust the position based on the actual graph. 
Rhizo.betaplot.ps2 <- Rhizo.betaplot.ps2 + annotation_custom(grob = textGrob(label = Rhizo.ps2.stress.label, hjust = 0, gp = gpar(fontsize = 20, fontface = "bold")), xmin = -0.65, xmax = -0.65, ymin = 0.5, ymax = 0.5)
#add PERMANOVA results on the graph....adjust the position based on the actual graph. 
Rhizo.betaplot.ps2 <- Rhizo.betaplot.ps2 + annotation_custom(grob = textGrob(label = Rhizo.ps2.adn_res_format, hjust = 0, gp = gpar(fontsize = 20, fontface = "bold")), xmin = -0.6, xmax = -0.6, ymin = 0.35, ymax = 0.35)
#add the results of environment fit. Plot arrows representing each variable. Red color indicates significant association
Rhizo.betaplot.ps2 <- Rhizo.betaplot.ps2 + geom_segment(data = ef.Rhizo.ps2.df,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow = arrow(length = unit(0.25, "cm")), colour = ef.Rhizo.ps2.df$color, inherit.aes = FALSE) +  geom_text(data = ef.Rhizo.ps2.df, aes(x = NMDS1*1.05, y = NMDS2*1.05, label = EnvironFactor), size = 3, inherit.aes = FALSE)
ggsave(filename = "beta-plot-Rhizo.pdf", plot = Rhizo.betaplot.ps2, device = "pdf", width = width*1.1, height = height*1.1, units = "mm")
#plotting stressplot for calculated NMDS ordinates
stressplot(Rhizo.ps2.bray.ordination)
#ordisurf
Rhizo.ps2.ordisurf.df=data.frame(x=Rhizo.ps2.bray.ordination$point[,1],y=Rhizo.ps2.bray.ordination$point[,2],Groups=metadata.rhizo[,"Stage"], Shapes = metadata.rhizo[,"CollectionTimepoint"], Description = metadata.rhizo[,"Description"])
Rhizo.ps2.ordisurf.plot <- list()
for (i in 1:length(variables))
{
  var <- variables[i]
  ordi<- vegan::ordisurf(Rhizo.ps2.bray.ordination,metadata.rhizo[,eval(var)],plot = FALSE, bs="ds")
  ordi.grid <- ordi$grid #extracts the ordisurf object
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  ordi.mite.na <- data.frame(na.omit(ordi.mite))
  ordiplot <- ggplot() + stat_contour(data = ordi.mite.na, mapping = aes(x = x, y = y, z = z, colour = ..level..),positon="identity", size = 1) + scale_colour_gradientn(colours = myPalette(100)) + ggplot2::labs(colour = paste(var)) + xlab("NMDS1") + ylab("NMDS2") + theme(aspect.ratio=1) + new_scale_color() + geom_point(data=Rhizo.ps2.ordisurf.df, aes(x,y,shape = Shapes, color = Groups),size=3, inherit.aes = FALSE) + scale_shape_manual(values=Rhizo.shapes, guide=FALSE) + geom_convexhull(inherit.aes = FALSE, data = Rhizo.ps2.ordisurf.df, mapping = aes(x = x, y = y, fill = Description), alpha=0.1, show.legend = FALSE) + ggtitle(label = paste("Ordisurf for envirnomental parameter ",var, sep = "")) + guides(colour = guide_legend(order = 1))
  Rhizo.ps2.ordisurf.plot[[eval(var)]] <- ordiplot
}
Rhizo.ps2.ordisurf.plot.all <- ggarrange(nrow = 4, ncol = 3, plotlist = Rhizo.ps2.ordisurf.plot)
ggsave(filename = "beta-plot-ordisurf-rhizo.pdf", plot = Rhizo.ps2.ordisurf.plot.all, device = "pdf", width = height * 2, height = width * 2, units = "mm")

#pairwise adonis for bray-curtis distance
Rhizo.pairwise.adonis.bray <- pairwise.adonis(x = Rhizo.ps2.bray.dist, factors = metadata.rhizo$CollectionTimepoint)
#modify dataframe to be upper matrix
Rhizo.pairwise.adonis.bray$pairs <- as.character(Rhizo.pairwise.adonis.bray$pairs)
Rhizo.pairwise.adonis.bray <- separate(data = Rhizo.pairwise.adonis.bray[,c(1,6)], col = pairs, into = c("g1", "g2"), sep = " vs ")
Rhizo.pairwise.adonis.bray <- spread(data = Rhizo.pairwise.adonis.bray, g2, p.value)
rownames(Rhizo.pairwise.adonis.bray) <- Rhizo.pairwise.adonis.bray[,1]
Rhizo.pairwise.adonis.bray <- Rhizo.pairwise.adonis.bray[,-1]
write.table(x = Rhizo.pairwise.adonis.bray, file = "Rhizo.pairwise.adonis.bray.txt", quote = FALSE, sep = "\t", na = "", row.names = TRUE, col.names = TRUE)
##


#####core microbiome
#core microbiome for genera
#genus level glommed phyloseq object already exists. However, need to modify the taxa names, so will make a new phyloseq object
ps2Rhizogenus <- Rhizo.ps2.Genus.rel
Rhizo.taxa.genus <- tax_table(ps2Rhizogenus)
Rhizo.taxa.genus <- as.data.frame(Rhizo.taxa.genus@.Data)
colname <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
Rhizo.taxa.genus[colname] <- sapply(Rhizo.taxa.genus[colname],as.character)
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Phylum),]$Phylum <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Phylum),]$Kingdom, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Class),]$Class <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Class),]$Phylum, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Order),]$Order <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Order),]$Class, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Family),]$Family <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Family),]$Order, "_X")
Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Genus),]$Genus <- paste0(Rhizo.taxa.genus[is.na(Rhizo.taxa.genus$Genus),]$Family, "_X")
taxa_names(ps2Rhizogenus) <- Rhizo.taxa.genus[, 6]
Rhizo.genus.core <- plot_core(x = ps2Rhizogenus, min.prevalence = 0.4, prevalences = seq(.05, 1, .05), detections = round(10^seq(log10(0.001), log10(1), length = 50), digits = 4), plot.type = 'heatmap', colours = rev(brewer.pal(5, "Spectral"))) +  labs(x = "Detection Threshold", y = "Genus")  + ggtitle("Genera distribution at minimum prevalence 0.4 across rhizosphere samples") + theme(axis.text.x = element_text(size = 5))
ggsave(filename = "Core-microbiome-rhizosphere-genus.pdf", plot = Rhizo.genus.core, device = "pdf", width = height *1.5, height = width*1.5, units = "mm")
#Extract the core microbiome genera for correlations
rhizo.core <- core(ps2Rhizogenus, detection = 0.001, prevalence = .4)

#correlation among all core genera
Rhizo.ps2Genus.df <- psmelt(rhizo.core)
Rhizo.genus.corr <- Rhizo.ps2Genus.df[,c("Sample", "OTU", "Abundance")]
Rhizo.genus.corr <- spread(data = Rhizo.genus.corr, key = "OTU", value = "Abundance")
rownames(Rhizo.genus.corr) <- Rhizo.genus.corr$Sample
Rhizo.genus.corr <- Rhizo.genus.corr[,-1]
res <- cor(Rhizo.genus.corr)
res2<- Hmisc::rcorr(as.matrix(Rhizo.genus.corr))
pdf(file = "Rhizosphere-genus-0.001-correlation.pdf", width = height/15, height = width/15)
corrplot(res2$r, p.mat = res2$P, sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 0.5, order = "hclust", method = "square", addrect = 6)
dev.off()
##



#####ASV significance based on DESeq2 based log2FoldChange
#between types of samples
ds = phyloseq_to_deseq2(ps2, ~ Type)
ds = DESeq(ds)
alpha = 0.01
levels.comp <- c("PreSowing", "Rhizosphere", "Bulk", "PostHarvest")
type.deseq <- list()
for ( i in 1:(length(levels.comp)-1))
{
  for (j in (i+1):length(levels.comp))
  {
    comp1 <- levels.comp[i]
    comp2 <- levels.comp[j]
    
    res = results(ds, contrast=c("Type", comp1, comp2), alpha=alpha)
    res = res[order(res$padj, na.last=NA), ]
    res_sig = res[(res$padj < alpha), ]
    res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps2)[rownames(res_sig), ], "matrix"))
    label.t <- paste0(comp1, " vs ", comp2, "  (Number of ASVs = ", nrow(res_sig), ")")
    list.name <- paste(comp1,comp2, sep = "_")
    type.deseq[[list.name]] <- ggplot(res_sig, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_hline(yintercept = c(10,-10)) + geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme(legend.position = "none") + ggtitle(label = label.t) + geom_text_repel(aes(label=ifelse(log2FoldChange >= 10 | log2FoldChange <= -10,ifelse(test = (is.na(Genus)), as.character(paste(Family,"_X", sep = "")), as.character(Genus) ),'')),hjust=0,vjust=0)
    
   }
}
type.deseq.plot = ggarrange(plotlist = type.deseq[c(1,4,2,5,3,6)], ncol = 2, nrow = 3, align = "hv")
ggsave(filename = "Deseq-type.pdf", plot = type.deseq.plot, device = "pdf", width = height*2.5, height = width*2.5, units = "mm")

#between Rhizosphere and bulk samples of respective collections
ds.description = phyloseq_to_deseq2(ps2, ~ Description)
ds.description = DESeq(ds.description)
description.deseq <- list()
description.deseq.df <- list()
for (j in 1:8)
  {
    comp1 <- paste("Rhizosphere",j, sep = "")
    comp2 <- paste("Bulk",j, sep = "")
    
    res = results(ds.description, contrast=c("Description", comp1, comp2), alpha=alpha)
    res = res[order(res$padj, na.last=NA), ]
    res_sig = res[(res$padj < alpha), ]
    res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps2)[rownames(res_sig), ], "matrix"))
#    label.t <- paste0(comp1, " vs ", comp2, "  (Number of ASVs = ", nrow(res_sig), ")")
#    list.name <- paste(comp1,comp2, sep = "_")
#    description.deseq[[list.name]] <- ggplot(res_sig, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_hline(yintercept = c(10,-10)) + geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme(legend.position = "none") + ggtitle(label = label.t) + geom_text_repel(aes(label=ifelse(log2FoldChange >= 10 | log2FoldChange <= -10,ifelse(test = (is.na(Genus)), as.character(paste(Family,"_X", sep = "")), as.character(Genus) ),'')),hjust=0,vjust=0)
    res_sig.df <- ifelse(test = (is.na(res_sig$Family)), as.character(paste(ifelse(test = (is.na(res_sig$Order)), as.character(paste(ifelse(test = (is.na(res_sig$Class)), as.character(paste(ifelse(test = (is.na(res_sig$Phylum)), as.character(paste(res_sig$Kingdom,"_X", sep = "")), as.character(res_sig$Phylum) ),"_X", sep = "")), as.character(res_sig$Class) ),"_X", sep = "")), as.character(res_sig$Order) ),"_X", sep = "")), as.character(res_sig$Family))
    list.df.name <- paste0("Collection",j)
    description.deseq.df[[list.df.name]] <- as.character(res_sig.df)
}
description.deseq.plot = ggarrange(plotlist = description.deseq, ncol = 2, nrow = 4, align = "hv")
ggsave(filename = "Deseq-description.pdf", plot = description.deseq.plot, device = "pdf", width = height*3, height = width*3, units = "mm")

#between crop development stages for rhizosphere samples
ds.stage = phyloseq_to_deseq2(subset_samples(ps2, Type != "Bulk"), ~ Stage)
ds.stage = DESeq(ds.stage)
alpha = 0.01
levels.comp.stage <- c("PreSowing", "Seedling", "PreNodulating", "Nodulating", "Flowering", "Matured", "Postharvest")
stage.deseq <- list()
for ( i in 1:(length(levels.comp.stage)-1))
{
  for (j in (i+1):length(levels.comp.stage))
  {
    comp1 <- levels.comp.stage[i]
    comp2 <- levels.comp.stage[j]
    
    res = results(ds.stage, contrast=c("Stage", comp1, comp2), alpha=alpha)
    res = res[order(res$padj, na.last=NA), ]
    res_sig = res[(res$padj < alpha), ]
    res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps2)[rownames(res_sig), ], "matrix"))
    label.t <- paste0(comp1, " vs ", comp2, "  (Number of ASVs = ", nrow(res_sig), ")")
    list.name <- paste(comp1,comp2, sep = "_")
    stage.deseq[[list.name]] <- ggplot(res_sig, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_hline(yintercept = c(10,-10)) + geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme(legend.position = "none") + ggtitle(label = label.t) + geom_text_repel(aes(label=ifelse(log2FoldChange >= 10 | log2FoldChange <= -10, ifelse(test = (is.na(Genus)), as.character(paste(Family,"_X", sep = "")), as.character(Genus) ),'')), max.overlaps = 100, force_pull = 0, force = 0.5)
     
  }
}

stage.deseq.plot = ggarrange(plotlist = stage.deseq[c(1,7,12,16,19,21)], ncol = 2, nrow = 3, align = "hv")
ggsave(filename = "Deseq-stage.pdf", plot = stage.deseq.plot, device = "pdf", width = height*2.5, height = width*2.5, units = "mm")

