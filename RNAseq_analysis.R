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
  Contrast_description = 'Comparison of WT_0 vs bioF_50 in LB/NGM') %>%
  left_join(info) %>%
  mutate(entrezid = unlist(entrezid)) %>% 
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())






results.complete = res.WT.tidy %>% rbind(res.WTN.tidy, res.B.tidy, res.B_WT.tidy,
                                         res.G_WT.tidy,
                                         res.Bmet_G.tidy, res.Bmet_WT.tidy)



# write results in excel files
list_of_datasets = list('WT_LB/NGM' = res.WT.tidy, 
                        'WT_NGM' = res.WTN.tidy, 
                        'bioF_LB/NGM' = res.B.tidy,
                        'bioF_0vsWT_0' = res.B_WT.tidy,
                        'gacA_0vsWT_0' = res.G_WT.tidy,
                        'bioF_50vsgacA_0' = res.Bmet_G.tidy,
                        'WT_0vsbioF_50' = res.Bmet_WT.tidy)


write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'),
           colNames = T, rowNames = F)

# write.xlsx(res.WT.tidy, here('summary', 'stats.xlsx'),
#            colNames = T, rowNames = F) 

write_csv(results.complete, here('summary', 'complete_stats.csv'))




#### MA plots ####


### MA plots for every comparison
plotMA(res.WT,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_WT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.WTN,  ylim=c(-3,3),  alpha = 0.05)
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
  scale_fill_brewer(palette = "Dark2") + 
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





# # # # # # # # # # # # # # # # # # # # #
### Scatter plots ########################
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


