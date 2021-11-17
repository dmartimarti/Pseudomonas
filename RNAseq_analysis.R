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
library(glue)

theme_set(theme_cowplot())


# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

samples = read.delim("sampleInfo.txt") 

samples = samples %>% 
  mutate(Batch = Replicate,
         Batch = as.factor(Batch),
         Sample = as.factor(Sample))

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





# DESeq2 analysis ---------------------------------------------------------



### starting analysis with DESeq2
# create DESeq data type to be analysed
# let's introduce a batch effect just in case
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, 
                                  design = ~ Batch + Sample)

# # prefilter, but that might not be necessary
keep = rowSums(counts(ddsTxi)) >= 10
ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "WT_0")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)
  

### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% 
  cbind(gene_list,.) %>% as_tibble()

gene_counts = gene_counts %>% 
  gather(Name, counts, W0L1:G0L4) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) %>%
  mutate(Sample = factor(Sample, levels = c('WT_0', 'WT_50', 'WTN_0', 'WTN_50',
                                            'B_0', 'B_50', 'G_0'))) # refactor levels to show them in this order in the plots


write_csv(gene_counts, here('summary', 'gene_counts.csv'))


## correct for batch effect 
## this way we can show the plots in a more nicer way

gene_counts_norm = limma::removeBatchEffect(counts(dds, normalized = TRUE), dds$batch) 
gene_list = rownames(gene_counts_norm)
gene_counts_norm = gene_counts_norm %>% 
  cbind(gene_list,.) %>% as_tibble()

gene_counts_norm = gene_counts_norm %>% 
  gather(Name, counts, W0L1:G0L4) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) %>%
  mutate(Sample = factor(Sample, levels = c('WT_0', 'WT_50', 'WTN_0', 'WTN_50',
                                            'B_0', 'B_50', 'G_0'))) # refactor levels to show them in this order in the plots


write_csv(gene_counts_norm, here('summary', 'gene_counts_batch_normalized.csv'))




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



pcaData %>%
  mutate(replicate = str_sub(name, -1),
         replicate = factor(replicate, levels = c(1,2,3,4))) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Groups)) + 
  geom_point(aes(shape = replicate), size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, linetype = Groups), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, fill = Groups), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_cowplot(14) 
  
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld_replicates.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


#### limma batch effect ####

dds$batch = factor(rep(c("1", "2", "3", "4")))
vsd = varianceStabilizingTransformation(dds)
plotPCA(vsd, "batch")

assay(vsd) = limma::removeBatchEffect(assay(vsd), vsd$batch)    
plotPCA(vsd, intgroup = c("Sample"))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_vsd_batch_effect.pdf'),
             height = 8, width = 9, useDingbats = FALSE)

# custom PCA plot

pcaData = plotPCA(vsd, intgroup = c("Bacteria", "Metformin", "Media"), returnData = TRUE)
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
             file = here('summary', 'PCA_main_vld_batch_effect.pdf'),
             height = 8, width = 9, useDingbats = FALSE)



# Stats MAIN --------------------------------------------------------------

# get results and tidy it
res = results(dds) 


# results with different shape of contrasts, tidy

# Wt50 - Wt0
res.WT = results(dds,   contrast = c("Sample", "WT_50" , "WT_0"))  
res.WT = lfcShrink(dds, contrast = c("Sample", "WT_50" , "WT_0"), res = res.WT, type = 'ashr')

# WTN50 - WTN0
res.WTN = results(dds,  contrast = c("Sample",  "WTN_50", "WTN_0")) 
res.WTN = lfcShrink(dds, contrast = c("Sample",  "WTN_50", "WTN_0"), res = res.WTN, type = 'ashr')

# bioF: B_50 - B_0
res.B = results(dds,  contrast = c("Sample",  "B_50", "B_0"))   
res.B = lfcShrink(dds, contrast = c("Sample", "B_50", "B_0"), res = res.B, type = 'ashr')

# B_0 - WT_0
res.B_WT = results(dds, contrast = c("Sample",   "B_0", "WT_0")) 
res.B_WT = lfcShrink(dds, contrast = c("Sample", "B_0", "WT_0"), res = res.B_WT, type = 'ashr')

# G_0 - WT_0
res.G_WT = results(dds, contrast = c("Sample",   "G_0", "WT_0")) 
res.G_WT = lfcShrink(dds, contrast = c("Sample", "G_0", "WT_0"), res = res.G_WT, type = 'ashr')

# B_50 - G_0
res.Bmet_G = results(dds, contrast = c("Sample",   "B_50", "G_0")) 
res.Bmet_G = lfcShrink(dds, contrast = c("Sample", "B_50", "G_0"), res = res.Bmet_G, type = 'ashr')

# B_50 - WT_0
res.Bmet_WT = results(dds, contrast = c("Sample",   "B_50", "WT_0")) 
res.Bmet_WT = lfcShrink(dds, contrast = c("Sample", "B_50", "WT_0"), res = res.Bmet_WT, type = 'ashr')


# WTN50 - WT0
res.WTNmet_WT = results(dds,  contrast = c("Sample",  "WTN_50", "WT_0")) 
res.WTNmet_WT = lfcShrink(dds, contrast = c("Sample",  "WTN_50", "WT_0"), res = res.WTNmet_WT, type = 'ashr')



# Wt0 - WtN0 (LB effect w/o metformin)
res.NGM = results(dds,   contrast = c("Sample", "WT_0" , "WTN_0"))  
res.NGM = lfcShrink(dds, contrast = c("Sample", "WT_0" , "WTN_0"), res = res.NGM, type = 'ashr')

# WT50 - WTN50 (LB effect with metformin)
res.NGM_met = results(dds,  contrast = c("Sample",  "WT_50", "WTN_50")) 
res.NGM_met = lfcShrink(dds, contrast = c("Sample",  "WT_50", "WTN_50"), res = res.NGM_met, type = 'ashr')



#### tidying the results #####
# Wt50 - Wt0
res.WT.tidy = as_tibble(res.WT, rownames = 'gene_id') %>% mutate(
      p_adj_stars = gtools::stars.pval(padj),
      Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
      Contrast = 'WT',
      Media = 'LB/NGM',
      Target = 'WT_50',
      Reference = 'WT_0',
      Contrast_description = 'Comparison of WT (+- Metformin) in LB/NGM') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

# WTN50 - WTN0
res.WTN.tidy = as_tibble(res.WTN, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'WTN',
  Media = 'NGM',
  Target = 'WTN_50',
  Reference = 'WTN_0',
  Contrast_description = 'Comparison of WT (+- Metformin) in NGM') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

# bioF: B_50 - B_0
res.B.tidy = as_tibble(res.B, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'bioF',
  Media = 'LB/NGM',
  Target = 'B_50',
  Reference = 'B_0',
  Contrast_description = 'Comparison of bioF (+- Metformin) in LB/NGM') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# B_0 - WT_0
res.B_WT.tidy = as_tibble(res.B_WT, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'B_WT',
  Media = 'LB/NGM',
  Target = 'B_0',
  Reference = 'WT_0',
  Contrast_description = 'Comparison of bioF_0 vs WT_0 in LB/NGM') %>%
  # mutate(entrezid = unlist(entrezid)) %>% 
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# G_0 - WT_0
res.G_WT.tidy = as_tibble(res.G_WT, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'G_WT',
  Media = 'LB/NGM',
  Target = 'G_0',
  Reference = 'WT_0',
  Contrast_description = 'Comparison of gacA_0 vs WT_0 in LB/NGM') %>% 
  # mutate(entrezid = unlist(entrezid)) %>% 
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# B_50 - G_0
res.Bmet_G.tidy = as_tibble(res.Bmet_G, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'Bmet_G',
  Media = 'LB/NGM',
  Target = 'B_50',
  Reference = 'G_0',
  Contrast_description = 'Comparison of bioF_50 vs gacA_0 in LB/NGM') %>%
  # mutate(entrezid = unlist(entrezid)) %>% 
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# B_50 - WT_0
res.Bmet_WT.tidy = as_tibble(res.Bmet_WT, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'Bmet_WT',
  Media = 'LB/NGM',
  Target = 'B_50',
  Reference = 'WT_0',
  Contrast_description = 'Comparison of bioF_50 vs WT_0 in LB/NGM') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

# WTN50 - WT0
res.WTNmet_WT.tidy = as_tibble(res.WTNmet_WT, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'WTNmet_WT',
  Media = 'LB/NGM',
  Target = 'WTN_50',
  Reference = 'WT_0',
  Contrast_description = 'Comparison of WTN50 vs WT_0') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# Wt0 - WtN0 (LB effect w/o metformin)
res.NGM.tidy = as_tibble(res.NGM, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'NGM_0',
  Media = 'LB/NGM, NGM',
  Target = 'WT_0',
  Reference = 'WtN_0',
  Contrast_description = 'Comparison of media w/o metformin (WT0 - WTN0)') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# WT50 - WTN50 (LB effect with metformin)
res.NGM_met.tidy = as_tibble(res.NGM_met, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'NGM_50',
  Media = 'LB/NGM, NGM',
  Target = 'WT_50',
  Reference = 'WtN_50',
  Contrast_description = 'Comparison of media with metformin (WT0 - WTN0)') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())


# gather all results
results.complete = res.WT.tidy %>% rbind(res.WTN.tidy,res.WTNmet_WT.tidy, 
                                         res.B.tidy, res.B_WT.tidy,
                                         res.G_WT.tidy,
                                         res.Bmet_G.tidy, res.Bmet_WT.tidy,
                                         res.NGM.tidy, res.NGM_met.tidy)



# write results in excel files
list_of_datasets = list('WT_LB/NGM' = res.WT.tidy, 
                        'WT_NGM' = res.WTN.tidy, 
                        'WTNmet_WT' = res.WTNmet_WT.tidy,
                        'bioF_LB/NGM' = res.B.tidy,
                        'bioF_0vsWT_0' = res.B_WT.tidy,
                        'gacA_0vsWT_0' = res.G_WT.tidy,
                        'bioF_50vsgacA_0' = res.Bmet_G.tidy,
                        'WT_0vsbioF_50' = res.Bmet_WT.tidy,
                        'WTN_50 vs WT_0'= res.WTNmet_WT.tidy ,
                        'WT0 vs WTN0 (NGM effect)' = res.NGM.tidy,
                        'WT50 vs WTN50 (NGM effect)' = res.NGM_met.tidy)


write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'),
           colNames = T, rowNames = F, overwrite = TRUE)


write_csv(results.complete, here('summary', 'complete_stats.csv'))



#### MA plots ####


### MA plots for every comparison
plotMA(res.WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.WTN,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_WTNmet_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.WTNmet_WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_WTN.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.B,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_B.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.B_WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_B_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.G_WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_G_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.Bmet_G,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_Bmet_G.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.Bmet_WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_Bmet_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)






# interactions in nutrition ------------------------------------------------------------

# load input data, again, but with only interaction terms

samples_int = read.delim("sampleInfo_interactions.txt") 

samples_int = samples_int %>% mutate(
  Metformin = as.factor(Metformin),
  Media = as.factor(Media),
  Sample = as.factor(Sample),
  Batch = as.factor(Batch)
)

dir = getwd()
rownames(samples_int) = samples_int$Name

# prepare a list with file names
files = file.path(dir,"quants", samples_int$Name, "quant.sf")
names(files) = samples_int$Name
all(file.exists(files)) # check that files exist

# import quantification data 
txi_int = tximport(files, type = "salmon", tx2gene = tx2gene)

ddsTxi = DESeqDataSetFromTximport(txi_int, colData = samples_int, 
                                  design = ~ Batch + Media + Metformin + Media:Metformin)

# # prefilter, but that might not be necessary
keep = rowSums(counts(ddsTxi)) >= 10
ddsTxi = ddsTxi[keep,]

ddsTxi$Media = relevel(ddsTxi$Media, ref = "LB/NGM")
# ddsTxi$Metformin = relevel(ddsTxi$Metformin, ref = 0)
# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds_int = DESeq(ddsTxi)

resultsNames(dds_int)

int_res = results(dds_int, name="MediaNGM.Metformin50")

int_res.tidy = as_tibble(int_res, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% 
  drop_na(padj) %>% 
  left_join(info)


d = plotCounts(dds_int, gene='WBGene00017501', intgroup="Media",
               returnData = T)

ggplot(d, aes(x=Media, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) 



gns = int_res.tidy %>% arrange(padj) %>% 
  head(6) %>% pull(gene_id)


gns = int_res.tidy %>% 
  filter(padj < 0.01)%>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  filter(baseMean > 50)  %>% 
  head(15) %>% pull(gene_id)



gns = c('WBGene00003766')

gene_counts%>%
  dplyr::filter(gene_id %in% gns) %>%
  filter(Sample %in% c('WT_0','WT_50','WTN_0','WTN_50')) %>% 
  ggplot(aes(y = counts, x = Sample)) +
  geom_boxplot(aes(fill = Sample),
               outlier.colour = NULL,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(~gene_name, scales = 'free_y') +
  labs(x = 'Sample',
       y = 'Normalised counts (log scale)') +
  theme_cowplot(15) +
  panel_border() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))


#### get the individual stats ####

# ~ Batch + Media + Metformin + Media:Metformin
# ddsTxi$Media = relevel(ddsTxi$Media, ref = "LB/NGM")

resultsNames(dds_int)

# interaction
int_res = results(dds_int, name="MediaNGM.Metformin50")

int_res.tidy = as_tibble(int_res, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% 
  drop_na(padj) %>% 
  left_join(info) %>% 
  mutate(contrast = 'interaction', .before = 'Direction')


# WT_50 vs WT_0
wt_int_res = results(dds_int, contrast = c('Metformin', '50','0'))

wt_int_res.tidy = as_tibble(wt_int_res, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% 
  drop_na(padj) %>% 
  left_join(info) %>% 
  mutate(contrast = 'WT', .before = 'Direction')

# 
wtn_int_res = results(dds_int, list( c("Metformin_50_vs_0","MediaNGM.Metformin50") ))

wtn_int_res.tidy = as_tibble(wtn_int_res, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% 
  drop_na(padj) %>% 
  left_join(info) %>% 
  mutate(contrast = 'WTN', .before = 'Direction')





list_of_datasets = list(
  'WT'=wt_int_res.tidy %>% select(-entrezid),
  'WTN'=wtn_int_res.tidy %>% select(-entrezid),
  'interaction'=int_res.tidy %>% select(-entrezid)
)



write.xlsx(list_of_datasets, here('summary', 'interaction_stats.xlsx'),
           overwrite = T, colNames = T, rowNames = F,)


interaction_results = int_res.tidy %>% bind_rows(wt_int_res.tidy, wtn_int_res.tidy) %>% 
  select(-entrezid) 


interaction_results %>% 
  write_csv(here('summary', 'interaction_stats.csv'))


#### interaction plots ####
gns = c('WBGene00003766')

gns_stats = interaction_results %>% 
  filter(gene_id %in% gns) %>% 
  mutate(p_adj_stars = ifelse(p_adj_stars == ' ', 'no sig', p_adj_stars))

gns_counts = gene_counts %>%
  dplyr::filter(gene_id %in% gns) %>%
  filter(Sample %in% c('WT_0','WT_50','WTN_0','WTN_50')) %>% 
  mutate(max_val = max(counts))

max_val = gns_counts$max_val[1]

fc.wt = gns_stats %>% filter(contrast=='WT') %>% pull(log2FoldChange) %>% round(2)
fc.wtn = gns_stats %>% filter(contrast=='WTN') %>% pull(log2FoldChange)%>% round(2)
fc.int = gns_stats %>% filter(contrast=='interaction') %>% pull(log2FoldChange)%>% round(2)

pval.wt = gns_stats %>% filter(contrast=='WT') %>% pull(p_adj_stars)
pval.wtn = gns_stats %>% filter(contrast=='WTN') %>% pull(p_adj_stars)
pval.int = gns_stats %>% filter(contrast=='interaction') %>% pull(p_adj_stars)


gns_counts %>% 
  ggplot(aes(y = counts, x = Sample)) +
  geom_boxplot(aes(fill = Sample),
               outlier.colour = NULL,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 6,
    shape = 21,
    fill = "blue"
  ) +
  # WT segment
  geom_segment(aes(x = 1, y = max_val*1.05, 
                   xend = 2, yend = max_val*1.05),
               colour = 'black', size = 1) +
  annotate('text', x = 1.5, y = max_val*1.07,
           label = glue('log2FC:{fc.wt}, pval:{pval.wt}')) +
  # WTN segment
  geom_segment(aes(x = 3, y = max_val*1.05, 
                   xend = 4, yend = max_val*1.05),
               colour = 'black', size = 1) +
  annotate('text', x = 3.5, y = max_val*1.07,
           label = glue('log2FC:{fc.wtn}, pval:{pval.wtn}')) +
  # Int segment
  geom_segment(aes(x = 1.5, y = max_val*1.15, 
                   xend = 3.5, yend = max_val*1.15),
               colour = 'red', size = 1) +
  annotate('text', x = 2.5, y = max_val*1.17,
           label = glue('log2FC:{fc.int}, pval:{pval.int}')) +
  facet_wrap(~gene_name, scales = 'free_y') +
  labs(x = 'Sample',
       y = 'Normalised counts (log scale)') +
  theme_cowplot(15) +
  panel_border() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))



# DESeq2 analysis ---------------------------------------------------------



### starting analysis with DESeq2
# create DESeq data type to be analysed
# let's introduce a batch effect just in case
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, 
                                  design = ~ Batch + Sample)

# # prefilter, but that might not be necessary
keep = rowSums(counts(ddsTxi)) >= 10
ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "WT_0")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)





### Exploratory plots/analysis ############################# 

# barplot that shows the % of DE genes in the 4 main comparisons
results.complete %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(Sig = ifelse(padj < 0.05, 1, 0)) %>%
  group_by(Contrast, Sig) %>%
  summarise(N = n()) %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100, 2)) %>%
  dplyr::filter(Sig == 1) %>%
  ggplot(aes(x = fct_reorder(Contrast, Fraction, .desc = TRUE), y = Fraction)) +
  geom_bar(stat = 'identity', width = 0.5, aes(fill = Contrast)) +
  # scale_fill_brewer(palette = "Dark2") + 
  scale_y_continuous(limits = c(0, 50), 
                     breaks = c(0, 15, 30, 50),
                     expand = expansion(mult = c(0, 0.05))) +
  geom_text(aes(label = Fraction, y = (Fraction + 1))) +
  labs(y = '% of DE genes',
       x = 'Condition') +
  guides(fill = 'none') +
  theme_cowplot(14)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'DEgenes_barplot.pdf'),
             height = 8, width = 8, useDingbats = FALSE)




# barplot that shows the % of DE genes in the 4 main comparisons


results.complete %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(Sig = ifelse(padj < 0.05, 1, 0)) %>%
  group_by(Contrast_description, Contrast, Direction, Sig) %>%
  summarise(N = n()) %>%
  group_by(Contrast_description, Contrast) %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1)) %>%
  group_by(Contrast_description, Contrast,  Sig) %>%
  arrange(desc(Direction)) %>%
  mutate(label_ypos = cumsum(Fraction)) %>%
  dplyr::filter(Sig == 1) %>%
  ggplot(aes(x = fct_reorder(Contrast, Fraction, .desc = TRUE), y = Fraction, fill = Direction)) +
  geom_bar(stat = 'identity', width = 0.5) +
  geom_text(aes(y = label_ypos, label = Fraction), 
            vjust = 1.6, size = 3.5) +
  scale_fill_brewer(palette = "Dark2") + 
  # scale_x_discrete(limits = c('N2', 'skpo', 'OP50', 'OG1RF')) +	
  scale_y_continuous(limits = c(0, 50), 
                     breaks = c(0, 25, 50),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(y = '% of DE genes',
       x = 'Condition') +
  theme_cowplot(14)
  

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'DEgenes_direction_barplot.pdf'),
             height = 8, width = 8, useDingbats = FALSE)


#### # # # ## # # ## # 
#### Gene boxplots ####
#### # # # ## # # ## # 


# argK-1 
gns = c('WBGene00009706')

gene_counts%>%
  dplyr::filter(gene_id %in% gns) %>%
  ggplot(aes(y = counts, x = Sample)) +
  geom_boxplot(aes(fill = Sample),
               outlier.colour = NULL,
               outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(~gene_name, scales = 'free_y') +
  labs(x = 'Sample',
       y = 'Normalised counts') +
  theme_cowplot(15) +
  panel_border() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))


# quartz.save(file = here('summary', 'nol-6.pdf'),
#             type = 'pdf', dpi = 300, height = 7, width = 6)
# 


### small function to help me
gplot = function(genes = c('WBGene00009706'), fw_nrows = 1){
  gene_counts_norm %>% 
    dplyr::filter(gene_id %in% genes) %>%
    ggplot(aes(y = counts, x = Sample)) +
    geom_boxplot(aes(fill = Sample),
                 outlier.colour = NULL,
                 outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2)) +
    facet_wrap(~gene_name, scales = 'free_y', nrow = fw_nrows) +
    labs(x = 'Sample',
         y = 'Normalised counts') +
    theme_cowplot(15) +
    panel_border() +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
    guides(fill = 'none')
}



gplot(c('WBGene00009706', 'WBGene00000185'), fw_nrows = 2)



#### interesting genes ####

# which genes share a similar behaviour as argk-1

selected_genes = results.complete %>% 
  filter(Contrast %in% c('WT', 'WTN', 'bioF')) %>% 
  filter(padj < 0.05, Direction == 'Up', log2FoldChange > 1) %>% 
  select(gene_name, gene_id, Contrast, log2FoldChange) %>% 
  pivot_wider(names_from = Contrast, values_from = log2FoldChange) %>% 
  drop_na() 

gplot(genes = selected_genes%>% pull(gene_id), fw_nrows = 4)

ggsave(file = here('summary', 'genes_of_interest.pdf'),
            height = 10, width = 12)









#### Volcano plots ####

temp_results = results.complete %>% 
  # dplyr::filter(Contrast %in% c('WT')) %>%
  mutate(padj_log = -log10(padj),
         significant = ifelse(padj < 0.05, 'Significant', ''),
         log2FC = ifelse(abs(log2FoldChange) > 3.5, 3.5, log2FoldChange),
         logPval = ifelse(padj_log > 15, 15, padj_log)) %>% 
  drop_na(significant) 


target = 'WT_50'
ref = 'WT_0'

temp_results %>% 
  dplyr::filter(Contrast %in% c('WT')) %>%
  ggplot(aes(x = log2FC, y = logPval, 
             text = paste('Gene: ',gene_name))) +
  geom_point(aes(color = significant), alpha = 0.7) +
  scale_colour_manual(values=c('#BABABA', '#0000FF')) +
  labs(x = glue::glue('log<sub>2</sub> Fold Change<br>Contrast: {target} vs {ref}'),
       y = '-log<sub>10</sub>(P-value adj)',
       title = '',
       color = NULL) +
  theme_cowplot(15) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown()) +
  guides(color = 'none')


ggplotly(p,  width=700, height=800) %>% 
  layout(
    title = list(text="Volcano plot of comparison <br> WT vs WT 50",
                 size = 4))

temp_results %>% 
  dplyr::filter(Contrast %in% c('WT')) %>%
  plot_ly() %>% 
  add_trace(x = ~log2FC,
            y = ~logPval,
            type = 'scatter',
            mode = 'markers', 
            # text = ~gene_name,
            color = ~significant, 
            colors = c('#BABABA', '#0000FF'),
            hovertemplate = paste('<i>Gene</i>: {~gene_name}',
                                  '<br>log2FC: %{log2FC}',
                                  '<br>log2Pval: %{logPval}</b>'),
            showlegend = FALSE
            )
  



## create a function to plot this

volcano = function(
  contrast = 'WT',
  padj_thres = 15,
  log2_thres = 3.5){
  
  temp = results.complete %>% 
    dplyr::filter(Contrast %in% contrast) %>%
    mutate(padj_log = -log10(padj),
           significant = ifelse(padj < 0.05, 'Significant', ''),
           log2FC = ifelse(abs(log2FoldChange) > log2_thres, log2_thres, log2FoldChange),
           logPval = ifelse(padj_log > padj_thres, padj_thres, padj_log)) %>% 
    drop_na(significant)
  
  title = temp %>% distinct(Contrast_description) %>% pull(Contrast_description)
  
  p = temp %>% 
    ggplot(aes(x = log2FC, y = logPval, 
               text = paste('Gene: ',gene_name))) +
    geom_point(aes(color = significant), alpha = 0.7) +
    scale_colour_manual(values=c('#BABABA', '#0000FF')) +
    labs(x = 'log<sub>2</sub> Fold Change',
         y = '-log<sub>10</sub>(P-value adj)',
         title = title,
         color = NULL) +
    theme_cowplot(15) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    guides(color = 'none')
  
  
  ggplotly(p,  width=700, height=800) %>% 
    layout(
      title = list(size = 4)) 
}


volcano()

volcano(contrast = 'bioF', padj_thres = 50, log2_thres = 9)



# # # # # # # # # # # # # # # # # #
### IMPORTANT: gene directions ####
# # # # # # # # # # # # # # # # # #

# positive threshold
pos = 0
# negative threshold
neg = -0

WT.gns.up = res.WT.tidy %>% 
  filter(padj <= 0.05, log2FoldChange >= pos) %>% 
  pull(gene_id)
WT.gns.down = res.WT.tidy %>% 
  filter(padj <= 0.05, log2FoldChange <= neg) %>% 
  pull(gene_id)
write.table(c('Wormbase.ID', unique(WT.gns.up)), 
            here('summary/gene_list', 'WT_UP_genes.txt'), quote = FALSE, col.names = F, row.names = F)
write.table(c('Wormbase.ID', unique(WT.gns.down)), 'WT_DOWN_genes.txt', quote = FALSE, col.names = F, row.names = F)





# # # # # # # # # # # # # # # # # # # # #
### Scatter plots #######################
# # # # # # # # # # # # # # # # # # # # #

# to plot regression info 
ggplotRegression = function(fit){
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


## WT vs WTN 

cosa = results.complete %>%
  filter(Contrast %in% c('WT', 'WTN'))

gns = cosa %>% filter(padj <=1) %>% select(gene_name) %>% unique %>% t %>% as.character

cosa = cosa %>% 
  filter(gene_name %in% gns) %>%
  select(Contrast, gene_name, log2FoldChange) %>%
  pivot_wider(names_from = Contrast, values_from = log2FoldChange) %>%
  data.frame

# ggplotRegression(lm(skpo ~ N2, data = cosa)) + theme_classic() +
#   ylim(-10,10) + 
#   xlim(-10,10)

model = summary(lm(cosa[,3] ~ cosa[,2]))
stats = paste0('y = ', signif(model$coef[[2]], 5) , ' * x + (', signif(model$coef[[1]], 5), '), \n R2 = ', signif(model$adj.r.squared, 5), ' P-value < 0.0001')

# cosa = cosa %>%
#   mutate(N2 = ifelse(N2 > 5, 4.99, N2),
#          N2 = ifelse(N2 < -5, -4.99, N2),
#          skpo = ifelse(skpo > 5, 4.99, skpo),
#          skpo = ifelse(skpo < -5, -4.99, skpo))

ggplot(cosa, aes(x = WT, y = WTN)) +
  geom_smooth(method = lm) +
  geom_point(alpha = 0.5) +
  ylim(-5,5) +
  xlim(-5,5) +
  geom_text(x = -3, y = 4, label = stats) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  # labs(x = expression(paste('N2 on', italic(' E. faecalis'),'/N2 on', italic(' E. coli'), '\n (Log' ['2'], 'FC)')),
  #      y = expression(paste(italic('skpo-1'),' on', italic(' E. faecalis'),'/',italic('skpo-1'),' on', italic(' E. coli'), '\n (Log' ['2'], 'FC)'))) +
  theme_light()


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'scatter_WT_vs_WTN.pdf'),
             height = 6, width = 8, useDingbats = FALSE)





# STRING DB enrich --------------------------------------------------------


# prepare the gene lists to feed my script

# comparisons that I made
unique(results.complete$Contrast)

# function to extract gene sets
gene_sets = function(
  contrast = 'WT',     # contrast to extract
  fc_thres_pos = 0.5,  # log2FC positive threshold
  fc_thres_neg = -0.5, # log2FC negative threshold
  base_thres = 50      # base threshold
  ){
  positive = results.complete %>% 
    filter(Contrast == contrast, log2FoldChange > fc_thres_pos, 
           padj < 0.05, baseMean > base_thres) %>% 
    arrange(desc(log2FoldChange)) %>% head(500) %>% 
    pull(gene_name)
  
  negative = results.complete %>% 
    filter(Contrast == contrast, log2FoldChange < fc_thres_neg, 
           padj < 0.05, baseMean > base_thres) %>% 
    arrange(log2FoldChange) %>% head(500) %>%  
    pull(gene_name)
  
  out_list = list()
  
  out_list[[paste0(contrast,'_UP')]] = positive
  out_list[[paste0(contrast,'_DOWN')]] = negative
  
  return(out_list)
  }


# loop to get the gene sets from the different contrasts
study_contrasts = unique(results.complete$Contrast)
list_of_genes = list()
for (contrast in study_contrasts){
  up = gene_sets(contrast=contrast)[[1]]
  down = gene_sets(contrast=contrast)[[2]]
  
  list_of_genes[[glue::glue('{contrast}_UP')]] = up
  list_of_genes[[glue::glue('{contrast}_DOWN')]] = down
  
}


write.xlsx(list_of_genes, here('summary', 'multi_gene_selection.xlsx'),
           overwrite = TRUE)




# GSEA analysis -----------------------------------------------------------

# Enrichment Browser

library(EnrichmentBrowser)

# user function to read samples from txt
read.samples = function(samp_file = "sampleInfo.txt", quants = 'quants') {
  
  samples = read.delim(samp_file) 
  
  dir = getwd()
  rownames(samples) = samples$Name
  
  quants_dir = quants
  
  # prepare a list with file names
  files = file.path(dir,quants_dir, samples$Name, "quant.sf")
  names(files) = samples$Name
  
  return(files)
}



# user function to get the DESeq2 summarizeExperiment from sample info

get.experiment <- function(
  thres=40, # count threshold to filter out genes
  sampleInfo = 'sampleInfo_WT.txt', # sample info file with samples to load
  ref = 'WT_0' # reference of the comparison
) {
  
  cat(glue::glue('Reading file {sampleInfo} \n\n'))
  
  # need to do again DESeq2 for hct per separate
  files.wt = read.samples(samp_file = sampleInfo)
  
  cat('Importing files to the system!\n')
  
  # import quantification data
  txi.wt = tximport::tximport(files.wt, type = "salmon", tx2gene = tx2gene)
  
  samples.red = read.delim(sampleInfo) 
  ddsTxi.wt = DESeqDataSetFromTximport(txi.wt,
                                       colData = samples.red,
                                       design = ~ Batch + Sample)
  
  cat(glue::glue('Filtering rows with less than {thres} counts.\n'))
  
  # filtering by min number of sequences
  keep = rowSums(counts(ddsTxi.wt)) >= thres
  ddsTxi.wt = ddsTxi.wt[keep,]
  
  ddsTxi.wt$Sample = relevel(ddsTxi.wt$Sample, ref = ref)
  
  cat('Running DDSeq2 pipeline...\n')
  # run the Differential Expression Analysis
  dds.wt = DESeq(ddsTxi.wt)
  
  dds.wt.pure = results(dds.wt)
  
  cat('Importing gene names and normalising data.\n')
  # import for gene enrichment
  WT.SE = import(dds.wt, dds.wt.pure, from = c('DESeq2'), anno = 'cel')
  # map IDs
  WT.SE = idMap(WT.SE, org = "cel", from = "ENSEMBL", to = "ENTREZID")
  
  # normalize counts
  WT.SE = normalize(WT.SE, norm.method = "vst")
  
  cat('Finished!\n')
  return(WT.SE)
}



## Run this once!

# obtaining gene sets
kegg.gs = getGenesets(org = "cel", db = "kegg")
go.gs = getGenesets(org = "cel", db = "go", onto = "BP")




#### WT ####

WT.SE = get.experiment()

head(rownames(WT.SE))

# set based enrichment analysis
sbeaMethods()

# run GSEA enrichment
# kegg 
alpha = 0.2
wt.gsea = sbea(method = "gsea", se = WT.SE, gs = kegg.gs, alpha = alpha)
wt.ora = sbea(method = "ora", se = WT.SE, gs = kegg.gs, alpha = alpha)

wt.go.gsea = sbea(method = "gsea", se = WT.SE, gs = go.gs, alpha = 0.05)

gsRanking(wt.go.gsea)

gsRanking(wt.gsea)

gsRanking(wt.ora)

eaBrowse(wt.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/wt.gsea', 
         report.name = 'wt.gsea')

eaBrowse(wt.ora, html.only = FALSE, out.dir = 'EnrichmentBrowser/KEGG/wt.ora', 
         report.name = 'wt.ora')


eaBrowse(wt.go.gsea, html.only = FALSE, out.dir = 'EnrichmentBrowser/GO_BP/wt.go.gsea', 
         report.name = 'wt.go.gsea')


# network regulation analysis

cel.grn = compileGRN(org="cel", db="kegg")


nbeaMethods()

nbea.res = nbea(method="ggea", 
                se=WT.SE, 
                gs=kegg.gs, 
                grn=cel.grn)

gsRanking(nbea.res) %>% view




#### WT NGM ####

WTN.SE = get.experiment(sampleInfo = 'sampleInfo_WTN.txt',
                       ref = 'WTN_0')

# run GSEA enrichment
# kegg 
alpha = 0.2
wtn.gsea = sbea(method = "gsea", se = WTN.SE, gs = kegg.gs, alpha = alpha)
wtn.ora = sbea(method = "ora", se = WTN.SE, gs = kegg.gs, alpha = alpha)

wtn.go.gsea = sbea(method = "gsea", se = WTN.SE, gs = go.gs, alpha = 0.05)

gsRanking(wtn.gsea)
gsRanking(wtn.ora)
gsRanking(wtn.go.gsea)

eaBrowse(wtn.gsea, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/KEGG/wtn.gsea', 
         report.name = 'wtn.gsea')

eaBrowse(wtn.ora, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/KEGG/wtn.ora', 
         report.name = 'wtn.ora')


eaBrowse(wtn.go.gsea, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/GO_BP/wtn.go.gsea', 
         report.name = 'wtn.go.gsea')


# network regulation analysis

nbea.res = nbea(method="ggea", 
                se=WTN.SE, 
                gs=kegg.gs, 
                grn=cel.grn)

gsRanking(nbea.res) 






#### bioF  ####

bioF.SE = get.experiment(sampleInfo = 'sampleInfo_bioF.txt',
                        ref = 'B_0')

# run GSEA enrichment
# kegg 
alpha = 0.2
bioF.gsea = sbea(method = "gsea", se = bioF.SE, gs = kegg.gs, alpha = alpha)
bioF.ora = sbea(method = "ora", se = bioF.SE, gs = kegg.gs, alpha = alpha)

bioF.go.gsea = sbea(method = "gsea", se = bioF.SE, gs = go.gs, alpha = 0.05)

gsRanking(bioF.gsea)
gsRanking(bioF.ora)
gsRanking(bioF.go.gsea)

eaBrowse(bioF.gsea, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/KEGG/bioF.gsea', 
         report.name = 'bioF.gsea')

eaBrowse(bioF.ora, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/KEGG/bioF.ora', 
         report.name = 'bioF.ora')


eaBrowse(bioF.go.gsea, html.only = FALSE, 
         out.dir = 'EnrichmentBrowser/GO_BP/bioF.go.gsea', 
         report.name = 'bioF.go.gsea')


# network regulation analysis

nbea.res = nbea(method="ggea", 
                se=bioF.SE, 
                gs=kegg.gs, 
                grn=cel.grn)

gsRanking(nbea.res) 












