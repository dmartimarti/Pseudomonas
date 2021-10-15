# RNA seq analysis from worms fed with Pseudomonas aeruginosa 14 strains
# and with metformin treatment

# In this script we will analyse the RNA seq with the following conditions:
# 	- WT +- metformin (NGM media)
# 	- WT +- metformin (LB/NGM media)
# 	- bioF mutant +- metformin (LB/NGM media)
# 	- gacA mutant w/o metformin (LB/NGM media)

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts


### libraries ####
library(tximport)
library(tidyverse)
library(DESeq2)
# notice that DESeq2 library masks 'rename' function from dplyr 
# library(ensembldb) # use only if you are going to deal with db
library(here)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(viridis)
library(cowplot)


# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

# load kegg tables from wormenrichr
kegg = read.delim("KEGG_2019.txt", header = FALSE) 
kegg = kegg[,-2]


rownames(kegg) = kegg[,1] ; kegg = kegg[,-1]

# prepare a list with file names
files = file.path(dir,"quants", samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist

# create an object with all the reference genes and transcripts
# THIS DOES NOT SEEM TO BE WORKING
# txdb = TxDb.Celegans.UCSC.ce11.refGene::TxDb.Celegans.UCSC.ce11.refGene
# k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
# tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")


# let's make our database from ensembldb 
ah = AnnotationHub::AnnotationHub(localHub = FALSE)
ahDb = AnnotationHub::query(ah, pattern = c("Caenorhabditis elegans", "EnsDb", 98))
ahEdb = ahDb[[1]]
# generate the database 
tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>% 
  tbl_df() %>%
  dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
  unnest

# join transcription info with gene ids and entrezids
info.join = tx2gene.complete %>% 
  tbl_df() %>%
  dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
  left_join(info)

write.csv(info, here('summary','gene_ids_mapping.csv'))

# subset to have tx_id in first column, and gene_id in second
tx2gene = tx2gene.complete[,c(1,7)]

# import quantification data 
txi = tximport(files, type = "salmon", tx2gene = tx2gene)




### starting analysis with DESeq2
# create DESeq data type to be analysed
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, design = ~ Sample)

# # prefilter, but that might not be necessary
# keep = rowSums(counts(ddsTxi)) >= 10
# ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "WT_0")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)


### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% 
  cbind(gene_list,.) %>% tbl_df()

gene_counts = gene_counts %>% 
  gather(Name, counts, W0L1:G0L4) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(tbl_df(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) %>%
  mutate(Sample = factor(Sample, levels = c('WT_0', 'WT_50', 'WTN_0', 'WTN_50',
                                            'B_0', 'B_50', 'G_0'))) # refactor levels to show them in this order in the plots





### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
norm_gene_counts = counts(dds, normalized = TRUE)
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
g_counts = counts(dds, normalized = FALSE)

write.table(g_counts, here('summary','gene_counts.txt'), quote = FALSE, sep = '\t') 
write.table(norm_gene_counts, here('summary','norm_gene_counts.txt'), quote = FALSE, sep = '\t') 




# KEGG databases
# get databases for genes and pathways from KEGG
kegg.links.entrez = limma::getGeneKEGGLinks('cel', convert = TRUE) 
kegg.links.ids = limma::getGeneKEGGLinks('cel')
path.ids = limma::getKEGGPathwayNames('cel', remove.qualifier = TRUE)
kegg.links = cbind(kegg.links.entrez, kegg.links.ids[,1])
colnames(kegg.links) = c('entrezid', 'PathwayID', 'KEGG_genes')

kegg.links = kegg.links %>% 
  as_tibble %>% 
  mutate(entrezid = as.integer(entrezid)) %>% 
  left_join(path.ids) %>%
  mutate(PathwayID = str_replace(PathwayID, 'path:cel', ''))




# transofrm data

vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)


# plot differences between different transformation data methods
df = bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("mean", "sd")  

ggplot(df, aes(x = mean, y = sd)) + 
  geom_hex(bins = 100) +
  # coord_fixed() + 
  facet_grid( . ~ transformation) +
  theme_light()


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'transformation_comparison.pdf'),
             height = 4.5, width = 12, useDingbats = FALSE)


###
# Sample distances

sampleDists = dist(t(assay(rld)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
names = colnames(sampleDists) %>%
  # str_split('_', simplify = T) %>%
  data.frame %>% tbl_df() %>%
  unite(sample, X1, X2, sep = " - ") %>%
  dplyr::select(sample) %>%
  t %>% as.vector

colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
        col = colors)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Euclidean_distances_samples.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


### PCA ############################# 
# PCA data

plotPCA(rld, intgroup = c("Bacteria", "Metformin", "Media"))

pcaData = plotPCA(rld, intgroup = c("Bacteria", "Metformin", "Media"), returnData = TRUE)
pcaData

names(pcaData) = c('PC1', 'PC2', 'Groups', 'Bacteria', 'Metformin','Media', 'name')

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}



# get info for the ellipses
# ell = pcaData %>% group_by(Bacteria, Metformin, Media) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame
ell = pcaData %>% group_by(Groups) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame
# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
pcaData %>% 
  ggplot(aes(x = PC1, y = PC2, color = Groups)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, linetype = Groups), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, fill = Groups), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_cowplot(14) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)




# plot!

ell = ell %>% 
  mutate(Metformin = factor(Metformin, levels = c(0,50)))

pcaData %>% 
  as_tibble() %>% 
  mutate(Metformin = factor(Metformin)) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Bacteria, group = interaction(Bacteria, Metformin, Media))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Bacteria, Metformin, Media)), 
            size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Bacteria, Metformin, Media),
                               fill = Metformin , linetype=Metformin), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black')) 








