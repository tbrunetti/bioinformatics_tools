####################################################
########## LOAD LIBRARIES AND SET OPTIONS ##########
####################################################
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(patchwork)
options(bitmapType="cairo")
allowWGCNAThreads() # allow multi-threading for running WGCNA functions

####################################################
############## USER DEFINED VARIABLES ##############
####################################################
counts_data = "raw_counts_data.csv" # should be raw counts of comma-separated data with rownames as gene names and colum names as sample names
sample_metadata_data = "sample_metadata.csv" # should be sample metadata of comma-separated values with rownames as sample names and column names as metadata names
outdir = getwd()
min_counts = 10 # minimum number of counts in a gene per sample in order to be counted as true for gene filtering
min_fraction = 0.30 # minimum fraction of samples that need to have min_counts in order for gene to be retained
power_select=10 # ! User must update based on power calculation in analysis step 2 !
network_type = "signed" # signed, unsigned, or signed hybrid
threads = 10

##########################################################
##################### START ANALYSIS #####################
##########################################################
timestamp=str_replace_all(string = Sys.time(), pattern = " ", replacement = "_")
##########################################
###### 1. FILTER AND NORMALIZE DATA ######
##########################################
# read in raw count data
raw_counts <- read.delim(file = counts_data, header = T, sep = ",", row.names = 1)

# read in metadata
metadata <- read.delim(file = sample_metadata_data,  header = T, sep = ",", row.names = 1)

# rearrange row and column names to match order
metadata <- metadata[colnames(raw_counts),]

# exclude any genes or samples that are outliers as well
check_for_outliers <- WGCNA::goodSamplesGenes(t(raw_counts))
table(check_for_outliers$goodGenes)
raw_counts <- raw_counts[check_for_outliers$goodGenes == TRUE,]
# filter genes and only keep genes that have at least min_counts counts in at least min_fraction of the samples
filtered_raw_counts <- raw_counts[names(which(rowSums(raw_counts[,] >= min_counts) >= round(ncol(raw_counts)*min_fraction))),]

# Normalize data using DESeq2 variance stabilizing transform
de_obj <- DESeqDataSetFromMatrix(round(filtered_raw_counts),
                              metadata,
                              design = ~ 1) # not running DE analysis, no model necessary, just getting normalized counts via variable stabilizing transformation
de_norm_obj <- vst(de_obj)
normalized_data <- assay(de_norm_obj)

# transpose data for WGCNA
wgcna_mat = t(normalized_data)


##################################################################
######## 2. DETERMINE POWER VALUE FOR NETWORK CONSTUCTION ########
##################################################################
# powerVector is the default by WGCNA but provide a list/range of powers to test in order to find best power for analysis
# networkType: unsigned means the direction of correlation does not matter; genes are connected in a network weather it is a positive or negative correlation
#              signed means the direction of correlation does matter; whether a gene is + or - correlated is taken into account when deciding on a connection between genes
select_soft_thresh = pickSoftThreshold(data = wgcna_mat, 
                                       dataIsExpr = T, 
                                       powerVector = c(seq(1, 10, by = 1), seq(12, 50, by = 2)),
                                       RsquaredCut = 0.85, # default
                                       nBreaks = 10, # default
                                       networkType = network_type,
                                       corFnc = cor, # default
                                       corOptions = list(use ='p'), # default
                                       verbose = 5)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!! USER WILL NEED TO UPDATE power_select VARIABLE BASED ON GRAPHS BELOW !!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# plot data to find power threshold to use; similar to an elbow plot pick a power near the elbow curve
# aim to pick a threshold that maximizes Rsquared while minimizing mean K connectivity

# power for r-squared (pick something above the line, but don't overfit)
rsq <- ggplot(select_soft_thresh$fitIndices, aes(x = Power, y = SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.05) +
  geom_hline(yintercept = 0.8, color = 'orange', linetype = "dashed") +
  labs(x = 'Power', y = paste('Scale free toplogy model fit,', network_type, 'R-squared', sep = ' ')) +
  ggtitle("Power determination for network construction") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# power for mean connectivity (smaller is better, but don't overfit)
conn <- ggplot(select_soft_thresh$fitIndices, aes(x = Power, y = mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_x = 0.7) +
  labs(x = 'Power', y ='mean connectivity') +
  theme_minimal() + ggtitle("")

# plot graph
rsq/conn


##################################################################
############### 3. GENERATE NETWORK BASED ON POWER ###############
##################################################################
temp_cor <- cor # save current cor function namespace to temp var since will use WGCNA cor function
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

# Call the network topology analysis function
network <- blockwiseModules(datExpr = wgcna_mat, power = power_select, networkType = network_type, 
                            corType = "pearson",  # default
                            TOMType = network_type,
                            mergeCutHeight = 0.25, #threshold to use to merge similar modules, this is default
                            numericLabels = F, # modules will be color names instead of numbers
                            deepSplit = 2, # default
                            detectCutHeight = 0.995, # default
                            minModuleSize = 30, # default: min(20, ncol(datExpr)/2 )
                            maxBlockSize = 22000, # default: 5000 genes, based on gigbytes of memory available on computer (standard laptop ~8-10k genes, 16G, can handle ~20k genes, 32G can handle 30k genes)
                            randomSeed = 42,
                            nThreads = threads,
                            useCorOptionsThroughout = T,
                            verbose = 5)
cor <- temp_cor # Return cor function to original namespace

save(list = c("select_soft_thresh", "wgcna_mat", "filtered_raw_counts", "normalized_data", "metadata", "network"), file = paste(outdir, "WGCNA_analysis_", timestamp, ".Rdat", sep = ""))

##################################################################
############## 4. RETRIEVE MODULE EIGENGENES (MEs) ###############
##################################################################
# Retrieve Module Eigengenes (MEs)
module_eigengenes <- network$MEs
# number of genes in each MEs
table(network$colors)

# plot dendrograms and module colors before and after merging
# merging means some modules were very similar and if they hit the merge threshold listed above
# then they were merged into one module
# use the merged modules for further analysis, grey modules are genes that don't fit in any of the modules
# so they are assigned to the grey category
plotDendroAndColors(
  network$dendrograms[[1]],
  colors = cbind(network$unmergedColors, network$colors),
  c("modules before merge", "modules after merge"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# get the genes associated with the module colors based on the heatmap
genes_in_modules <- as.data.frame(network$colors)
names(genes_in_modules)[1] <- "MEcolor"
genes_in_modules$genes <- rownames(genes_in_modules)
# split genes into MEs for each file
genes_in_MEs <- split(genes_in_modules$genes, genes_in_modules$MEcolor)
for (me in names(genes_in_MEs)){
  writeLines(text = genes_in_MEs[[me]], con = paste(outdir, "WGCNA_", me, "_eigengene_gene_members.txt", sep = ""))
}


######################################################################
#### 5. ASSOCIATE MODULE EIGENGENE ASSIGNMENTS TO GROUPS/TRAITS  #####
######################################################################

# WARNING!  THIS IS NOT AUTOMATED! Must be manually set and worked through
# per experiment.  DO NOT RUN without changing code for the dataset you are 
# working on

# identifies modules that are associated to traits
# categorize traits by binarization to find associations
# create a separate column for each level

# remove spaces in metadata
metadata$timepoint <- str_replace_all(metadata$timepoint, pattern = " ", replace = "_")
# change reference level to Iso_6_hr
metadata$timepoint <- factor(metadata$timepoint, levels = c("Iso_6_hr", "0", "1_week", "12_hr", "24_hr", "3_hr", "30_min", "48_hr", "6_hr", "72_hr", "Iso_24_hr"))
# binarize all associations
binarized_assocs <- binarizeCategoricalColumns(metadata$timepoint, includePairwise = T, includeLevelVsAll = T, minCount = 1)

##################################################
############## vs all associations ###############
##################################################

# will also correlate each group versus with all remaining samples
binarized_assocs_vs_all <- binarized_assocs %>% select(ends_with(".vs.all")) # update to pull out the "all" associations

# add sample names to rownames; should be the same order as the metadata
rownames(binarized_assocs_vs_all) <- rownames(metadata)

# get number of genes and samples for the all calculations
totalSamples <- nrow(binarized_assocs)
totalGenes <- ncol(wgcna_mat)

# correlate the eigengenes to traits
# generates correlations and pvalues but note, not required for the visualization; this is just if you want the
# data calculated as a matrix; could have also modeled with linear model or maybe a logistic regression instead of a 0/1 correlation but it does get you similar resultseigengene_trait_assocs_ref_vs_all<- cor(module_eigengenes, binarized_assocs_vs_all, method = "pearson", use = "everything")
eigengene_trait_assocs_pvals_ref_vs_all <- corPvalueStudent(eigengene_trait_assocs_ref_vs_all, totalSamples) # the important thing here is the association significance not necessarily the strength of the correlation

# visualize
names(binarized_assocs_vs_all) <- str_remove(string = names(binarized_assocs_vs_all), pattern = "^data.") # change name of samples on x-axis to drop the data. string 
vis_assoc_data_all <- merge(module_eigengenes, binarized_assocs_vs_all, by = 'row.names') # need to merge the eigengene info from step 4 with the binarized associations
# drop the rownames column since that is not part of the heatmap calculation but make sure rownames are sample names
rownames(vis_assoc_data_all) <- vis_assoc_data_all$Row.names 
vis_assoc_data_all <- vis_assoc_data_all[,c(names(module_eigengenes), names(binarized_assocs_vs_all))] #select only the columns you want to visualize
# make the heatmap plot of the pvalues and correlations; note this is calculated on the gly, no need to calculate above unless you want the matrix calculated but this will auto calculate=d
CorLevelPlot(vis_assoc_data_all,
             x = names(binarized_assocs_vs_all), # col names of the traits
             y = names(module_eigengenes), # col names of the eigengene modules
             col = c("blue1", "skyblue", "white", "pink", "red"), 
             main = "Module eigengene associations with binarized phenotype trait comparisons", cexMain = 2.5,
             titleY = "Module Eigengene (ME)", rotTitleY = 90, cexTitleY = 1.5, fontLabY = 1,
             titleX = "binarized traits", cexTitleX = 1.5, fontLabX = 1, 
)

##################################################
############ vs Iso_6_hr associations ############
##################################################

# the user wants to compare each timepoint (0.5h - 1 week) to the isoflurane 6hr controls.
binarized_assocs_vs_iso6hr <- binarized_assocs %>% select(ends_with(".vs.Iso_6_hr")) # update to pull out the Iso_6_hr associations

# add sample names to rownames; should be the same order as the metadata
rownames(binarized_assocs_vs_iso6hr) <- rownames(metadata)

# correlate the eigengenes to traits
# generates correlations and pvalues but note, not required for the visualization; this is just if you want the
# data calculated as a matrix; could have also modeled with linear model or maybe a logistic regression instead of a 0/1 correlation but it does get you similar results
# Uses pairwise complete observations, meaning each pair of values will use the complete data available for that pair. 
# It may result in different numbers of observations for each pair. 
eigengene_trait_assocs_ref_vs_iso6hr<- cor(module_eigengenes, binarized_assocs_vs_iso6hr, method = "pearson", use = "pairwise.complete.obs") #
eigengene_trait_assocs_ref_vs_iso6hr_pvals = ""
initialized = F
for (comps in names(binarized_assocs_vs_iso6hr)){ # the for loop is required here because the total number of samples is not the same unlike the all comparison
  if (initialized == T){
    print(comps)
    tmp_df <- as.data.frame(eigengene_trait_assocs_ref_vs_iso6hr[,comps], drop = F)
    names(tmp_df) <- comps
    print(tmp_df)
    totalSamplesInCorr = colSums(is.na(binarized_assocs_vs_iso6hr[comps]) == FALSE)
    print(totalSamplesInCorr)
    tmp_pval <- corPvalueStudent(as.matrix(tmp_df), totalSamplesInCorr)
    eigengene_trait_assocs_ref_vs_iso6hr_pvals <- cbind(eigengene_trait_assocs_ref_vs_iso6hr_pvals, tmp_pval)
  }else{
    print(comps)
    tmp_df <- as.data.frame(eigengene_trait_assocs_ref_vs_iso6hr[,comps], drop = F)
    names(tmp_df) <- comps
    print(tmp_df)
    totalSamplesInCorr = colSums(is.na(binarized_assocs_vs_iso6hr[comps]) == FALSE)
    print(totalSamplesInCorr)
    eigengene_trait_assocs_ref_vs_iso6hr_pvals <- corPvalueStudent(as.matrix(tmp_df), totalSamplesInCorr)
    initialized = T
  }
}
eigengene_trait_assocs_pvals_ref_vs_iso6hr <- corPvalueStudent(eigengene_trait_assocs_ref_vs_iso6hr, totalSamples) # the important thing here is the association significance not necessarily the strength of the correlation



# visualize
names(binarized_assocs_vs_iso6hr) <- str_remove(string = names(binarized_assocs_vs_iso6hr), pattern = "^data.") # change name of samples on x-axis to drop the data. string 
rownames(binarized_assocs_vs_iso6hr) <- rownames(metadata) 
vis_assoc_data_iso6hr <- merge(module_eigengenes, binarized_assocs_vs_iso6hr, by = 'row.names') # need to merge the eigengene info from step 4 with the binarized associations
# drop the rownames column since that is not part of the heatmap calculation but make sure rownames are sample names
rownames(vis_assoc_data_iso6hr) <- vis_assoc_data_iso6hr$Row.names
vis_assoc_data_iso6hr <- vis_assoc_data_iso6hr[,c(names(module_eigengenes), names(binarized_assocs_vs_iso6hr))] #select only the columns you want to visualize
# make the heatmap plot of the pvalues and correlations; note this is calculated on the gly, no need to calculate above unless you want the matrix calculated but this will auto calculate=d
CorLevelPlot(vis_assoc_data_iso6hr,
             x = names(binarized_assocs_vs_iso6hr), # col names of the traits
             y = names(module_eigengenes),# col names of the eigengene modules
             col = c("blue1", "skyblue", "white", "pink", "red"), 
             corUSE = "pairwise.complete.obs", corFUN = "pearson",
             main = "Module eigengene associations with binarized phenotype trait comparisons", cexMain = 2.5,
             titleY = "Module Eigengene (ME)", rotTitleY = 90, cexTitleY = 1.5, fontLabY = 1,
             titleX = "binarized traits", cexTitleX = 1.5, fontLabX = 0.5, rotLabX = 45
)



######################################################################
################### 6. OPTIONAL - FIND HUB GENES  ####################
######################################################################

# name of module color to explore
module_color_of_interest = "yellow"

# get all the genes in that module
genes_in_modules %>% filter(`network$colors` == !!module_color_of_interest) %>% rownames() # for tidyverse !! treats the string as a variable substitution


## Find hub genes in the modules (genes with high module membership)
module_membership_strength <- cor(module_eigengenes, t(normalized_data), method = "pearson")
totalSamples = nrow(t(normalized_data))
module_membership_strength_pvals <- corPvalueStudent(module_membership_strength, totalSamples)

write.csv(x = t(module_membership_strength), file = "signed_pearson_correlation_weights_for_each_gene_contributing_to_each_eigengene_module.csv", row.names = T, col.names = T)
write.csv(x = t(module_membership_strength_pvals), file = "signed_pearson_correlation_pvalues_for_each_gene_contributing_to_each_eigengene_module.csv", row.names = T, col.names = T)

# get correlation and pvalues for comparisons associated with genes
# this calculates the pvalues; note that totalSamples will need to be udpated for the association; 
# in particular when not using the "all" association as NAs should not count towards the sample count when phenotypes are binarized
for (comp in colnames(binarized_assocs_vs_iso6hr)){
  print(paste("Running", comp))
  samples = rownames(binarized_assocs_vs_iso6hr[which(is.na(binarized_assocs_vs_iso6hr[,comp]) == FALSE),])
  counts = t(normalized_data)[samples,]
  gene_sig_corr <- cor(counts, binarized_assocs_vs_iso6hr[samples, comp], method = "pearson", use = "complete.obs")
  colnames(gene_sig_corr) <- "pearson_correlation_coeff"
  gene_sig_corr_pvals <- corPvalueStudent(gene_sig_corr, length(samples)) 
  colnames(gene_sig_corr_pvals) <- "pvalue"
  assocs <- merge(gene_sig_corr, gene_sig_corr_pvals, by = "row.names")
  colnames(assocs)[1] <- "gene"
  write.csv(x = assocs, file = paste(comp, "_gene_correlation_pval_matrix.csv", sep = ""), row.names = F, col.names = T)
}


# the top 25 genes that are significantly correlated with severe_vs_all; top 25 hub genes for a module of interest
gene_sig_corr_pvals %>% as.data.frame() %>% arrange(V1) %>% head(25)
