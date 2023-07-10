# Progenitor Analysis
library(pacman)
p_load(
  ChIPseeker,
  readxl,
  TxDb.Mmusculus.UCSC.mm39.refGene,
  GenomicRanges,
  DESeq2,
  clusterProfiler,
  EnhancedVolcano,
  org.Mm.eg.db,
  edgeR,
  biomaRt,
  sva,
  ggpubr,
  cowplot,
  stringr
)
# Call in meta data for each sample
meta_data_all <-
  read_excel("~/Dropbox (ASU)/Mac (2)/Documents/ATAC/SampleLog/All_Samples.xlsx")

# Ensure meta data has all the proper samples in order of sample ID
meta_data1 <- meta_data_all[which(
  meta_data_all$Genotype == 'WT' &
    (meta_data_all$Diet == "Control" |
       meta_data_all$Diet == "HFD")
), ]
meta_data2 <- meta_data_all[which(
  meta_data_all$Genotype == 'WT' &
    (meta_data_all$Diet == "HFD1wk" |
       meta_data_all$Diet == "HFD4wk")
), ]
meta_data3 <- meta_data_all[which(
  meta_data_all$Genotype == "DA" |
    meta_data_all$Genotype == "Cpt1a" |
    meta_data_all$Genotype == "Apc" |
    meta_data_all$Diet == "Rereversal" |
    meta_data_all$Diet == "fasted"
), ]
meta_data <- rbind(meta_data1, meta_data2, meta_data3)
# Import multibamcov matrix, if not yet generated refer to lines [61-62]
countMatrix <-
  read.delim(
    "~/Dropbox (ASU)/Mac (2)/Downloads/diet_isc_tac/diet_isc_tac_countMatrix",
    header = T
  )
countMatrix1 <-
  read.delim("~/Downloads/diet_isc_reversal/diet_isc_reversal_countMatrix",
             header = T)
countMatrix2 <- 
  read.delim("~/Dropbox (ASU)/Mac (3)/Documents/ATAC/remaining_countMatrix",
             header = T)
countMatrix <- cbind(countMatrix, countMatrix1[, -c(1:4)] ,countMatrix2[, -c(1:4)])
# Subset regions based on autosomes only
# countMatrix <-
#   countMatrix[which(countMatrix$chr != "chrY" &
#                       countMatrix$chr != "chrX"),]

rownames(countMatrix) <- countMatrix$name # set rownames()

# Define batches in data
batch <- as.factor(meta_data$batch)

# Adjust based on batches
adjusted <- ComBat_seq(as.matrix(countMatrix[,-c(1:4)]), batch = batch, group = NULL)

# Define covariates
diet <- as.factor(meta_data$Diet)
sex <- as.factor(meta_data$Sex)
genotype <- as.factor(meta_data$Genotype)
GFP <- as.factor(meta_data$GFP)

# Combine covariates in a data frame
covar_mat<-cbind(diet,sex,genotype,GFP)

# Generate an adjusted count matrix
adjusted_counts <-
  ComBat_seq(as.matrix(countMatrix[,-c(1:4)]),
             batch = batch,
             group = NULL,
             covar_mod = covar_mat)

# Use DESeq2 to determine differentially accessible regions
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = adjusted,
  colData = meta_data,
  design =  ~ batch + Sex + GFP + Genotype + Diet,
  tidy = F
)

# Perform differential expression analysis with the batch effect removed
dds <- DESeq2::DESeq(dds,
             test = "Wald")

dds$group <-
  factor(paste0(dds$Diet, dds$GFP, dds$Genotype)) # Build all combinations of groups
design(dds) <- ~ group
dds <- DESeq2::DESeq(dds)

# Use DESeq2 to determine differentially accessible regions
dds0 <- DESeq2::DESeqDataSetFromMatrix(
  countData = adjusted_counts,
  colData = meta_data,
  design =  ~ batch + Sex + GFP + Genotype + Diet,
  tidy = F
)

# Perform differential expression analysis with the batch effect removed
dds0 <- DESeq2::DESeq(dds0,
             test = "Wald")

dds0$group <-
  factor(paste0(dds0$Diet, dds0$GFP, dds0$Genotype)) # Build all combinations of groups
design(dds0) <- ~ group
dds0 <- DESeq2::DESeq(dds0)

levels(dds$group)
# Generate annotated peak files that can be applied to all results
annotatedPeak <- annotatePeak(
  peak = GRanges(countMatrix[, c(1:4)]),
  tssRegion = c(-5000, 5000),
  TxDb = TxDb.Mmusculus.UCSC.mm39.refGene,
  annoDb = "org.Mm.eg.db"
)
annotatedPeak_df <- as.data.frame(annotatedPeak)
annotatedPeak0 <- annotatePeak(
  peak = GRanges(countMatrix[, c(1:4)]),
  tssRegion = c(-5000, 5000),
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db"
)
annotatedPeak0_df <- as.data.frame(annotatedPeak0)

# Save DDS objects for ease of use in the future (can be found in "~/Documents/ATAC/crunchtime/")
# saveRDS(object = dds,file = "allATAC_dds.rds")
# saveRDS(object = dds0,file = "allATAC_dds0.rds")


# ContrastGen::HFD versus CONTROL ISC --------------------------------------------------
ISC_hfdVcon <-
  DESeq2::results(dds, contrast = c("group", "HFDhiWT", "ControlhiWT"))
summary(ISC_hfdVcon)
ISC_hfdVcon_df <- as.data.frame(ISC_hfdVcon)

ISC_hfdVcon0 <-
  DESeq2::results(dds0, contrast = c("group", "HFDhiWT", "ControlhiWT"))
summary(ISC_hfdVcon0, alpha = .1)
ISC_hfdVcon0_df <- as.data.frame(ISC_hfdVcon0)

annotatedPeak_ISC_hfdVcon_df <-
  cbind(annotatedPeak_df, ISC_hfdVcon_df)
annotatedPeak_ISC_hfdVcon0_df <-
  cbind(annotatedPeak_df, ISC_hfdVcon0)
annotatedPeak_ISC_hfdVcon00_df <-
  cbind(annotatedPeak0_df, ISC_hfdVcon0)
# ContrastGen::HFD DA- versus CONTROL WT ISC--------------------------------------------------
hfdDAVconWT <-
  results(dds, contrast = c("group", "HFDhiDA", "ControlhiWT"))
summary(hfdDAVconWT)
hfdDAVconWT_df <- as.data.frame(hfdDAVconWT)

hfdDAVconWT0 <-
  results(dds0, contrast = c("group", "HFDhiDA", "ControlhiWT"))
summary(hfdDAVconWT0)
hfdDAVconWT0_df <- as.data.frame(hfdDAVconWT0)

annotatedPeak_hfdDAVconWT_df <-
  cbind(annotatedPeak_df, hfdDAVconWT_df)
annotatedPeak_hfdDAVconWT0_df <-
  cbind(annotatedPeak_df, hfdDAVconWT0_df)
annotatedPeak_hfdDAVconWT00_df <-
  cbind(annotatedPeak0_df, hfdDAVconWT0_df)
# ContrastGen::CON DA- versus CONTROL WT ISC--------------------------------------------------
conDAVconWT <-
  results(dds, contrast = c("group", "ControlhiDA", "ControlhiWT"))
summary(conDAVconWT)
conDAVconWT_df <- as.data.frame(conDAVconWT)

conDAVconWT0 <-
  results(dds0, contrast = c("group", "ControlhiDA", "ControlhiWT"))
summary(conDAVconWT0)
conDAVconWT0_df <- as.data.frame(conDAVconWT0)

annotatedPeak_conDAVconWT_df <-
  cbind(annotatedPeak_df, conDAVconWT_df)
annotatedPeak_conDAVconWT0_df <-
  cbind(annotatedPeak_df, conDAVconWT0_df)
annotatedPeak_conDAVconWT_df <-
  cbind(annotatedPeak0_df, conDAVconWT0_df)
# ContrastGen::HFD DA- versus CONTROL DA ISC--------------------------------------------------
hfdDAVconDA <-
  results(dds, contrast = c("group", "HFDhiDA", "ControlhiDA"))
summary(hfdDAVconDA)
hfdDAVconDA_df <- as.data.frame(hfdDAVconDA)

hfdDAVconDA0 <-
  results(dds0, contrast = c("group", "HFDhiDA", "ControlhiDA"))
summary(hfdDAVconDA0)
hfdDAVconDA0_df <- as.data.frame(hfdDAVconDA0)

annotatedPeak_hfdDAVconDA_df <-
  cbind(annotatedPeak_df, hfdDAVconDA_df)
annotatedPeak_hfdDAVconDA0_df <-
  cbind(annotatedPeak_df, hfdDAVconDA0_df)
annotatedPeak_hfdDAVconDA00_df <-
  cbind(annotatedPeak0_df, hfdDAVconDA0_df)
# ContrastGen::HFD Cpt1a- versus CONTROL WT ISC--------------------------------------------------
hfdCpt1aVconWT <-
  results(dds, contrast = c("group", "HFDhiCpt1a", "ControlhiWT"))
summary(hfdCpt1aVconWT)
hfdCpt1aVconWT_df <- as.data.frame(hfdCpt1aVconWT)

hfdCpt1aVconWT0 <-
  results(dds0, contrast = c("group", "HFDhiCpt1a", "ControlhiWT"))
summary(hfdCpt1aVconWT0)
hfdCpt1aVconWT0_df <- as.data.frame(hfdCpt1aVconWT0)

annotatedPeak_hfdCpt1aVconWT_df <-
  cbind(annotatedPeak_df, hfdCpt1aVconWT_df)
annotatedPeak_hfdCpt1aVconWT0_df <-
  cbind(annotatedPeak_df, hfdCpt1aVconWT0_df)
annotatedPeak_hfdCpt1aVconWT00_df <-
  cbind(annotatedPeak0_df, hfdCpt1aVconWT0_df)
# ContrastGen::HFD Cpt1a- versus CONTROL Cpt1a-  ISC--------------------------------------------------
hfdCpt1aVconCpt1a <-
  results(dds, contrast = c("group", "HFDhiCpt1a", "ControlhiCpt1a"))
summary(hfdCpt1aVconCpt1a)
hfdCpt1aVconCpt1a_df <- as.data.frame(hfdCpt1aVconCpt1a)

hfdCpt1aVconCpt1a0 <-
  results(dds0, contrast = c("group", "HFDhiCpt1a", "ControlhiCpt1a"))
summary(hfdCpt1aVconCpt1a0)
hfdCpt1aVconCpt1a0_df <- as.data.frame(hfdCpt1aVconCpt1a0)

annotatedPeak_hfdCpt1aVconCpt1a_df <-
  cbind(annotatedPeak_df, hfdCpt1aVconCpt1a_df)
annotatedPeak_hfdCpt1aVconCpt1a0_df <-
  cbind(annotatedPeak_df, hfdCpt1aVconCpt1a0_df)
annotatedPeak_hfdCpt1aVconCpt1a00_df <-
  cbind(annotatedPeak0_df, hfdCpt1aVconCpt1a0_df)
# ContrastGen::CONTROL Cpt1a- versus CONTROL WT ISC --------------------------------------------------
conCpt1aVconWT <-
  results(dds, contrast = c("group", "ControlhiCpt1a", "ControlhiWT"))
summary(conCpt1aVconWT)
conCpt1aVconWT_df <- as.data.frame(conCpt1aVconWT)

conCpt1aVconWT0 <-
  results(dds0, contrast = c("group", "ControlhiCpt1a", "ControlhiWT"))
summary(conCpt1aVconWT0)
conCpt1aVconWT0_df <- as.data.frame(conCpt1aVconWT0)

annotatedPeak_conCpt1aVconWT_df <-
  cbind(annotatedPeak_df, conCpt1aVconWT_df)
annotatedPeak_conCpt1aVconWT0_df <-
  cbind(annotatedPeak_df, conCpt1aVconWT0_df)
annotatedPeak_conCpt1aVconWT00_df <-
  cbind(annotatedPeak0_df, conCpt1aVconWT0_df)
# ContrastGen::HFD WT versus HFD Cpt1a- ISC --------------------------------------------------
hfdWTVhfdCpt1a <-
  results(dds, contrast = c("group", "HFDhiWT", "HFDhiCpt1a"))
summary(hfdWTVhfdCpt1a)
hfdWTVhfdCpt1a_df <- as.data.frame(hfdWTVhfdCpt1a)

hfdWTVhfdCpt1a0 <-
  results(dds0, contrast = c("group", "HFDhiWT", "HFDhiCpt1a"))
summary(hfdWTVhfdCpt1a0)
hfdWTVhfdCpt1a0_df <- as.data.frame(hfdWTVhfdCpt1a0)

annotatedPeak_hfdWTVhfdCpt1a_df <-
  cbind(annotatedPeak_df, hfdWTVhfdCpt1a_df)
annotatedPeak_hfdWTVhfdCpt1a0_df <-
  cbind(annotatedPeak_df, hfdWTVhfdCpt1a0_df)
annotatedPeak_hfdWTVhfdCpt1a00_df <-
  cbind(annotatedPeak0_df, hfdWTVhfdCpt1a0_df)
# ContrastGen::HFD WT versus HFD DA- ISC --------------------------------------------------
hfdWTVhfdDA <-
  results(dds, contrast = c("group", "HFDhiWT", "HFDhiDA"))
summary(hfdWTVhfdDA)
hfdWTVhfdDA_df <- as.data.frame(hfdWTVhfdDA)

hfdWTVhfdDA0 <-
  results(dds0, contrast = c("group", "HFDhiWT", "HFDhiDA"))
summary(hfdWTVhfdDA0)
hfdWTVhfdDA0_df <- as.data.frame(hfdWTVhfdDA0)

annotatedPeak_hfdWTVhfdDA_df <-
  cbind(annotatedPeak_df, hfdWTVhfdDA_df)
annotatedPeak_hfdWTVhfdDA0_df <-
  cbind(annotatedPeak_df, hfdWTVhfdDA0_df)
annotatedPeak_hfdWTVhfdDA00_df <-
  cbind(annotatedPeak0_df, hfdWTVhfdDA0_df)
# ContrastGen::HFD APC versus CONTROL WT ISC  --------------------------------------------------
hfdAPCVconWT <-
  results(dds, contrast = c("group", "HFDhiApc", "ControlhiWT"))
summary(hfdAPCVconWT)
hfdAPCVconWT_df <- as.data.frame(hfdAPCVconWT)

hfdAPCVconWT0 <-
  results(dds0, contrast = c("group", "HFDhiApc", "ControlhiWT"))
summary(hfdAPCVconWT0)
hfdAPCVconWT0_df <- as.data.frame(hfdAPCVconWT0)

annotatedPeak_hfdAPCVconWT_df <-
  cbind(annotatedPeak_df, hfdAPCVconWT_df)
annotatedPeak_hfdAPCVconWT0_df <-
  cbind(annotatedPeak_df, hfdAPCVconWT0_df)
annotatedPeak_hfdAPCVconWT00_df <-
  cbind(annotatedPeak0_df, hfdAPCVconWT0_df)

# ContrastGen::CONTROL APC versus CONTROL WT ISC--------------------------------------------------
conAPCVconWT <-
  results(dds, contrast = c("group", "ControlhiApc", "ControlhiWT"))
summary(conAPCVconWT)
conAPCVconWT_df <- as.data.frame(conAPCVconWT)

conAPCVconWT0 <-
  results(dds0, contrast = c("group", "ControlhiApc", "ControlhiWT"))
summary(conAPCVconWT0)
conAPCVconWT0_df <- as.data.frame(conAPCVconWT0)

annotatedPeak_conAPCVconWT_df <-
  cbind(annotatedPeak_df, conAPCVconWT_df)
annotatedPeak_conAPCVconWT0_df <-
  cbind(annotatedPeak_df, conAPCVconWT0_df)
annotatedPeak_conAPCVconWT0_df <-
  cbind(annotatedPeak0_df, conAPCVconWT0_df)

# ContrastGen::HFD APC versus CONTROL APC ISC --------------------------------------------------
hfdAPCVconAPC <-
  results(dds, contrast = c("group", "HFDhiApc", "ControlhiApc"))
summary(hfdAPCVconAPC)
hfdAPCVconAPC_df <- as.data.frame(hfdAPCVconAPC)

hfdAPCVconAPC0 <-
  results(dds0, contrast = c("group", "HFDhiApc", "ControlhiApc"))
summary(hfdAPCVconAPC0)
hfdAPCVconAPC0_df <- as.data.frame(hfdAPCVconAPC0)

annotatedPeak_hfdAPCVconAPC_df <-
  cbind(annotatedPeak_df, hfdAPCVconAPC_df)
annotatedPeak_hfdAPCVconAPC0_df <-
  cbind(annotatedPeak_df, hfdAPCVconAPC0_df)
annotatedPeak_hfdAPCVconAPC00_df <-
  cbind(annotatedPeak0_df, hfdAPCVconAPC0_df)

# ContrastGen::FASTED versus CONTROL WT ------------------------------------------------
fastedVconWT <-
  results(dds, contrast = c("group", "fastedhiWT", "ControlhiWT"))
summary(fastedVconWT)
fastedVconWT_df <- as.data.frame(fastedVconWT)

fastedVconWT0 <-
  results(dds0, contrast = c("group", "fastedhiWT", "ControlhiWT"))
summary(fastedVconWT0)
fastedVconWT0_df <- as.data.frame(fastedVconWT0)

annotatedPeak_fastedVconWT_df <-
  cbind(annotatedPeak_df, fastedVconWT_df)
annotatedPeak_fastedVconWT0_df <-
  cbind(annotatedPeak_df, fastedVconWT0_df)
annotatedPeak_fastedVconWT00_df <-
  cbind(annotatedPeak0_df, fastedVconWT0_df)

# ContrastGen::FASTED versus HFD WT ------------------------------------------------
fastedVhfdWT <-
  results(dds, contrast = c("group", "fastedhiWT", "HFDhiWT"))
summary(fastedVhfdWT)
fastedVhfdWT_df <- as.data.frame(fastedVhfdWT)

fastedVhfdWT0 <-
  results(dds0, contrast = c("group", "fastedhiWT", "HFDhiWT"))
summary(fastedVhfdWT0)
fastedVhfdWT0_df <- as.data.frame(fastedVhfdWT0)

annotatedPeak_fastedVhfdWT_df <-
  cbind(annotatedPeak_df, fastedVhfdWT_df)
annotatedPeak_fastedVhfdWT0_df <-
  cbind(annotatedPeak_df, fastedVhfdWT0_df)
annotatedPeak_fastedVhfdWT00_df <-
  cbind(annotatedPeak0_df, fastedVhfdWT0_df)

# ContrastGen::HFD ISC versus HFD TAC ------------------------------------------------
hfdISCVhfdTAC <-
  results(dds, contrast = c("group", "HFDhiWT", "HFDlowWT"))
summary(hfdISCVhfdTAC)
hfdISCVhfdTAC_df <- as.data.frame(hfdISCVhfdTAC)

hfdISCVhfdTAC0 <-
  results(dds0, contrast = c("group", "HFDhiWT", "HFDlowWT"))
summary(hfdISCVhfdTAC0)
hfdISCVhfdTAC0_df <- as.data.frame(hfdISCVhfdTAC0)

annotatedPeak_hfdISCVhfdTAC_df <-
  cbind(annotatedPeak_df, hfdISCVhfdTAC_df)
annotatedPeak_hfdISCVhfdTAC0_df <-
  cbind(annotatedPeak_df, hfdISCVhfdTAC0_df)
annotatedPeak_hfdISCVhfdTAC00_df <-
  cbind(annotatedPeak0_df, hfdISCVhfdTAC0_df)

# ContrastGen::CONTROL ISC versus CONTROL TAC ------------------------------------------------
conISCVconTAC <-
  results(dds, contrast = c("group", "ControlhiWT", "ControllowWT"))
summary(conISCVconTAC)
conISCVconTAC_df <- as.data.frame(conISCVconTAC)

conISCVconTAC0 <-
  results(dds0, contrast = c("group", "ControlhiWT", "ControllowWT"))
summary(conISCVconTAC0)
conISCVconTAC0_df <- as.data.frame(conISCVconTAC0)

annotatedPeak_conISCVconTAC_df <-
  cbind(annotatedPeak_df, conISCVconTAC_df)
annotatedPeak_conISCVconTAC0_df <-
  cbind(annotatedPeak_df, conISCVconTAC0_df)
annotatedPeak_conISCVconTAC00_df <-
  cbind(annotatedPeak0_df, conISCVconTAC0_df)

# ContrastGen::HFD TAC versus CONTROL TAC ------------------------------------------------
hfdTACVconTAC <-
  results(dds, contrast = c("group", "HFDlowWT", "ControllowWT"))
summary(hfdTACVconTAC)
hfdTACVconTAC_df <- as.data.frame(hfdTACVconTAC)

hfdTACVconTAC0 <-
  results(dds0, contrast = c("group", "HFDlowWT", "ControllowWT"))
summary(hfdTACVconTAC0)
hfdTACVconTAC0_df <- as.data.frame(hfdTACVconTAC0)

annotatedPeak_hfdTACVconTAC_df <-
  cbind(annotatedPeak_df, hfdTACVconTAC_df)
annotatedPeak_hfdTACVconTAC0_df <-
  cbind(annotatedPeak_df, hfdTACVconTAC0_df)
annotatedPeak_hfdTACVconTAC00_df <-
  cbind(annotatedPeak0_df, hfdTACVconTAC0_df)


# ContrastGen::HFD -1WK versus CONTROL WT ----------------------------------------------
hfd1wkVconWT <-
  results(dds, contrast = c("group", "HFD1wkhiWT", "ControlhiWT"))
summary(hfd1wkVconWT)
hfd1wkVconWT_df <- as.data.frame(hfd1wkVconWT)

hfd1wkVconWT0 <-
  results(dds0, contrast = c("group", "HFD1wkhiWT", "ControlhiWT"))
summary(hfd1wkVconWT0)
hfd1wkVconWT0_df <- as.data.frame(hfd1wkVconWT0)

annotatedPeak_hfd1wkVconWT_df <-
  cbind(annotatedPeak_df, hfd1wkVconWT_df)
annotatedPeak_hfd1wkVconWT0_df <-
  cbind(annotatedPeak_df, hfd1wkVconWT0_df)
annotatedPeak_hfd1wkVconWT00_df <-
  cbind(annotatedPeak0_df, hfd1wkVconWT0_df)

# ContrastGen::HFD -4WK versus CONTROL WT ----------------------------------------------
hfd4wkVconWT <-
  results(dds, contrast = c("group", "HFD4wkhiWT", "ControlhiWT"))
summary(hfd4wkVconWT)
hfd4wkVconWT_df <- as.data.frame(hfd4wkVconWT)

hfd4wkVconWT0 <-
  results(dds0, contrast = c("group", "HFD4wkhiWT", "ControlhiWT"))
summary(hfd4wkVconWT0)
hfd4wkVconWT0_df <- as.data.frame(hfd4wkVconWT0)

annotatedPeak_hfd4wkVconWT_df <-
  cbind(annotatedPeak_df, hfd4wkVconWT_df)
annotatedPeak_hfd4wkVconWT0_df <-
  cbind(annotatedPeak_df, hfd4wkVconWT0_df)
annotatedPeak_hfd4wkVconWT00_df <-
  cbind(annotatedPeak0_df, hfd4wkVconWT0_df)


# ContrastGen::RE-REVERSAL versus CONTROL WT ------------------------------------------------
rereversalVconWT <-
  results(dds, contrast = c("group", "RereversalhiWT", "ControlhiWT"))
summary(rereversalVconWT)
rereversalVconWT_df <- as.data.frame(rereversalVconWT)

rereversalVconWT0 <-
  results(dds0, contrast = c("group", "RereversalhiWT", "ControlhiWT"))
summary(rereversalVconWT0)
rereversalVconWT0_df <- as.data.frame(rereversalVconWT0)

annotatedPeak_rereversalVconWT_df <-
  cbind(annotatedPeak_df, rereversalVconWT_df)
annotatedPeak_rereversalVconWT0_df <-
  cbind(annotatedPeak_df, rereversalVconWT0_df)
annotatedPeak_rereversalVconWT00_df <-
  cbind(annotatedPeak0_df, rereversalVconWT0_df)


# PlotCounts::ISC_HFDvCON_Analysis ----------------------------------------------
# Save plotcounts to a data frame object
d <- plotCounts(dds0, gene=10311, intgroup="group", returnData=TRUE)

# Define the desired order of the groups
group_order <- c("ControlhiWT", "HFDhiWT") # Progenitor analysis group

color_palette <- c("#E69F00", "#56B4E9") # Progenitor colors (CON, HFD, CON-TA, HFD-TA)

# Filter the dataset to include only the desired groups
filtered_d <- d %>% filter(group %in% group_order)

# Reorder the levels of the 'group' variable
filtered_d$group <- factor(filtered_d$group, levels = group_order)

# Plotting the normalized counts at given peaks
ggplot(filtered_d, aes(x = group, y = count, color = group)) + 
  geom_boxplot(linewidth = 2,notch = F) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
  #geom_text(label = rownames(filtered_d)) +
  
  # Progenitor comparisons
  geom_signif(comparisons = list(c("ControlhiWT", "HFDhiWT")),
              map_signif_level = TRUE, textsize = 6, size = 1.25, color = "black", test = wilcox.test) +
  
  # Progenitor Labels
  scale_x_discrete(labels=c("ControlhiWT" = "CON", "HFDhiWT" = "HFD")) +
  
  labs(title = "Pdk4 (Promoter)", x = "", y="Normalized Counts per Sample") +
  scale_color_manual(values = color_palette) +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 4/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PlotCount_ISC_HFDvCON_Pdk4.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# PlotCounts::ProgenitorAnalysis ----------------------------------------------
# Save plotcounts to a data frame object
d <- plotCounts(dds0, gene=53593, intgroup="group", returnData=TRUE)

# Define the desired order of the groups
group_order <- c("ControlhiWT", "HFDhiWT", "ControllowWT", "HFDlowWT") # Progenitor analysis group

color_palette <- c("#E69F00", "#56B4E9", "#DB5A00", "#54D98C") # Progenitor colors (CON, HFD, CON-TA, HFD-TA)

# Filter the dataset to include only the desired groups
filtered_d <- d %>% filter(group %in% group_order)

# Reorder the levels of the 'group' variable
filtered_d$group <- factor(filtered_d$group, levels = group_order)

# Plotting the normalized counts at given peaks
ggplot(filtered_d, aes(x = group, y = count, color = group)) + 
  geom_boxplot(linewidth = 2,notch = F) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
  
  # Progenitor comparisons
  geom_signif(comparisons = list(c("ControlhiWT", "HFDhiWT"), c("ControllowWT", "HFDlowWT")),
              map_signif_level = TRUE, textsize = 6, size = 1.25, color = "black", test = t.test) +

    # Progenitor Labels
  scale_x_discrete(labels=c("ControlhiWT" = "CON\n\nISC", "HFDhiWT" = "HFD\n\nISC",
                            "ControllowWT" = "CON\n\nTAC", "HFDlowWT" = "HFD\n\nTAC")) +

  labs(title = "Degree of accessibility", subtitle = "Pdk4", x = "", y="Normalized Counts per Sample") +
  scale_color_manual(values = color_palette) +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 2.25/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))

# PlotCounts::KnockOutAnalysis ----------------------------------------------
# Save plotcounts to a data frame object
d <- plotCounts(dds0, gene=53593, intgroup="group", returnData=TRUE)

# Define the desired order of the groups
group_order <- c("ControlhiWT", "HFDhiWT", "HFDhiDA", "HFDhiCpt1a") # Progenitor analysis group

color_palette <- c("#E69F00", "#56B4E9", "salmon", "yellow4") # Progenitor colors (CON, HFD, CON-TA, HFD-TA)

# Filter the dataset to include only the desired groups
filtered_d <- d %>% filter(group %in% group_order)

# Reorder the levels of the 'group' variable
filtered_d$group <- factor(filtered_d$group, levels = group_order)

# Plotting the normalized counts at given peaks
ggplot(filtered_d, aes(x = group, y = count, color = group)) + 
  geom_boxplot(linewidth = 2,notch = F) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
  
  # Progenitor comparisons
  geom_signif(comparisons = list(c("ControlhiWT", "HFDhiWT"), c("ControlhiWT", "HFDhiDA"),c("ControlhiWT", "HFDhiCpt1a")),
              map_signif_level = TRUE, 
              textsize = 6, 
              size = 1.25, 
              color = "black", 
              test = t.test, y_position = c(325,335,345)) +
  
  # Progenitor Labels
  scale_x_discrete(labels=c("ControlhiWT" = "CON\n\nISC", "HFDhiWT" = "HFD\n\nISC",
                            "HFDhiDA" = "HFD\nDA", "HFDhiCpt1a" = "HFD\nCpt1a")) +
  
  labs(title = "Degree of accessibility", subtitle = "Pdk4", x = "", y="Normalized Counts per Sample") +
  scale_color_manual(values = color_palette) +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 2.25/1, #legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))

# PlotCounts::RereversalAnalysis ----------------------------------------------
# Save plotcounts to a data frame object
d <- plotCounts(dds0, gene=68463, intgroup="group", returnData=TRUE)

# Define the desired order of the groups
group_order <- c("ControlhiWT", "HFDhiWT", "HFD1wkhiWT", "HFD4wkhiWT", "RereversalhiWT") # Progenitor analysis group

color_palette <- c("#E69F00", "#56B4E9", "purple3", "gold2", "cyan") # Reversal colors (CON, HFD, HFD-1wk, HFD-4wk, Re-reversal)

# Filter the dataset to include only the desired groups
filtered_d <- d %>% filter(group %in% group_order)

# Reorder the levels of the 'group' variable
filtered_d$group <- factor(filtered_d$group, levels = group_order)

# Plotting the normalized counts at given peaks
ggplot(filtered_d, aes(x = group, y = count, color = group)) + 
  geom_boxplot(linewidth = 2,notch = F) +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +

  # Rereversal comparisons
  geom_signif(comparisons = list(c("ControlhiWT", "HFDhiWT"), c("ControlhiWT", "RereversalhiWT")),
              map_signif_level = TRUE, textsize = 6, size = 1.25, color = "black", test = t.test, y_position = c(65,70)) +

  # Rereversal Labels
  scale_x_discrete(labels=c("ControlhiWT" = "CON",
                            "HFDhiWT" = "HFD",
                            "HFD1wkhiWT" = "HFD\n-1wk",
                            "HFD4wkhiWT" = "HFD\n-4wk",
                            "RereversalhiWT" = "HFD\n-4wk\n+1wk")) +
  
  labs(title = "Degree of accessibility", subtitle = "Pdk4", x = "", y="Normalized Counts per Sample") +
  scale_color_manual(values = color_palette) +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1.5/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))


# PCA::Making HFD versus Control ISC --------------------------------------------------------------
group <- paste0(meta_data$Diet,meta_data$GFP,meta_data$Genotype)
meta_data$group <- group

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$Diet[c(1,3,5,7,9,11,13,15,17:22)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = `meta_data$Diet[c(1, 3, 5, 7, 9, 11, 13, 15, 17:22)]`)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = `meta_data$Diet[c(1, 3, 5, 7, 9, 11, 13, 15, 17:22)]`)) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (45.69%)", y = "PC2 (12.88%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFDvCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# PCA::Making HFD versus Control ISC & TAC --------------------------------------------------------------
# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1:22)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1:22)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = `meta_data$group[c(1:22)]`)) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = `meta_data$group[c(1:22)]`)) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "#DB5A00", "#56B4E9", "#54D98C")) +
  scale_fill_manual(values = c("#E69F00", "#DB5A00", "#56B4E9", "#54D98C")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD & Control ISC & TAC", x = "PC1 (32.18%)", y = "PC2 (18.30%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFDvCON_ISCandTAC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


# PCA::Making HFD, HFD1wk, HFD4wk versus Control ISC --------------------------------------------------------------
# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:33)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1,3,5,7,9,11,13,15,17:33)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:33)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:33)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "purple3", "gold2", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "purple3", "gold2", "#56B4E9")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD, HFD1wk, HFD4wk & Control ISC", x = "PC1 (34.36%)", y = "PC2 (9.127%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD1wk_HFD4wkvCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)



# PCA::Making HFD, HFD1wk, HFD4wk, Rereversal versus Control ISC --------------------------------------------------------------
# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:33,42:46)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1,3,5,7,9,11,13,15,17:33,42:46)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:33,42:46)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:33,42:46)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "purple3", "gold2", "#56B4E9","cyan")) +
  scale_fill_manual(values = c("#E69F00", "purple3", "gold2", "#56B4E9","cyan")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD, HFD1wk, HFD4wk, Rereversal & Control ISC", x = "PC1 (31.96%)", y = "PC2 (8.594%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD1wk_HFD4wk_rerev_vCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# PCA::Making HFD, HFD1wk, HFD4wk versus Control ISC & TAC --------------------------------------------------------------
# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1:33)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1:33)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1:33)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1:33)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "#DB5A00", "purple3", "gold2","#56B4E9", "#54D98C")) +
  scale_fill_manual(values = c("#E69F00", "#DB5A00", "purple3", "gold2","#56B4E9", "#54D98C")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD, HFD1wk, HFD4wk & Control ISC", x = "PC1 (26.44%)", y = "PC2 (15.85%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD,HFD1wk,HFD4wkvCON_ISC&TACS.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# PCA::Making HFD, HFD-DA versus Control ISC --------------------------------------------------------------

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22,37:39,41,48,50)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$Diet[c(1,3,5,7,9,11,13,15,17:22,37:39,41,48,50)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,37:39,41,48,50)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,37:39,41,48,50)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "salmon", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "salmon", "#56B4E9")) +  # Set fill colors manually
  #xlim(-22.5, 22.5) + ylim(-10, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (39.77%)", y = "PC2 (11.19%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD-DA_vCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)



# PCA::Making HFD, HFD-Cpt1a versus Control ISC --------------------------------------------------------------

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22,52,54,56,58)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$Diet[c(1,3,5,7,9,11,13,15,17:22,52,54,56,58)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,52,54,56,58)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,52,54,56,58)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "thistle3", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "thistle3", "#56B4E9")) +  # Set fill colors manually
  #xlim(-22.5, 22.5) + ylim(-10, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (39.77%)", y = "PC2 (11.19%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD-DA_vCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)




# PCA::Making HFD, Cpt1as versus Control ISC --------------------------------------------------------------

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22,51:58)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$Diet[c(1,3,5,7,9,11,13,15,17:22,51:58)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,51:58)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,51:58)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "thistle3", "#56B4E9", "thistle3")) +
  scale_fill_manual(values = c("#E69F00", "thistle3", "#56B4E9")) +  # Set fill colors manually
  #xlim(-22.5, 22.5) + ylim(-10, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (39.77%)", y = "PC2 (11.19%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD-DA_vCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)





# PCA::Making HFD, HFD-APC, Control-APC versus Control ISC --------------------------------------------------------------

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$pvalue < 0.05),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22,59:66)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1,3,5,7,9,11,13,15,17:22,59:66)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,59:66)])) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = meta_data$group[c(1,3,5,7,9,11,13,15,17:22,59:66)])) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("red4", "#E69F00", "navy", "#56B4E9")) +
  scale_fill_manual(values = c("red4", "#E69F00", "navy", "#56B4E9")) +  # Set fill colors manually
  #xlim(-22.5, 22.5) + ylim(-10, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (36.21%)", y = "PC2 (17.50%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFD_HFD-DA_vCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

components$`meta_data$Diet[c(1, 3, 5, 7, 9, 11, 13, 15, 17:22, 59:66)]`


# VolcanoPlot::HFD DA vs CONTROL DA ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdDAVconDA0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD/Control)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-4.5,4.5) +
  ylim(0,5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_DAvCON_DA_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD WT ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_ISC_hfdVcon00_df,
  lab = annotatedPeak_ISC_hfdVcon00_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(
    x = "Log2FoldChange (HFD/Control)",
    y = "-Log10p-value(adjusted)",
    title = "",
    subtitle = ""
  ) +
  xlim(-5, 5) +
  ylim(0, 10) +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFDvCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD DA ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdDAVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD/Control)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-7.5,3.75) +
  ylim(0,27.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_DAvCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::CONTROL DA ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_conDAVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Control DA/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-8.75,3.75) +
  ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotCON_DAvCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


# VolcanoPlot::HFD Cpt1a ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdCpt1aVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD Cpt1aKO/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  ylim(0,42.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_Cpt1avCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD Cpt1a ISC vs CONTROL Cpt1a ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdCpt1aVconCpt1a0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD Cpt1aKO/Control Cpt1aKO)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  ylim(0,27.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_Cpt1avCON_Cpt1a_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


# VolcanoPlot::CONTROL Cpt1a ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_conCpt1aVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Control Cpt1aKO/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  ylim(0,45) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotCON_Cpt1avCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::CONTROL DA ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_conDAVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Control DA-KO/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-8.75,3.75) +
  ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotCON_DAvCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD WT ISC vs HFD Cpt1a ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdWTVhfdCpt1a0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD WT/HFD Cpt1a-KO)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_WTvHFD_Cpt1a_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD WT ISC vs HFD DA ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdWTVhfdDA0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD WT/HFD DA-KO)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_WTvHFD_DA_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD APC ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdAPCVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD APC-KO/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_APCvCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::CONTROL APC ISC vs CONTROL WT ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_conAPCVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Control APC-KO/Control WT)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotCON_APCvCON_WT_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


# VolcanoPlot::HFD APC ISC vs CONTROL APC ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdAPCVconAPC0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD APC-KO/Control APC-KO)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_APCvCON_APC_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::FASTED ISC vs CONTROL ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_fastedVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Fasted (24hr)/Control)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotFastedvCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::FASTED ISC vs HFD ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_fastedVhfdWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Fasted (24hr)/HFD)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-8.75,3.75) +
  # ylim(0,30) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotFastedvHFD_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD ISC vs HFD TAC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdISCVhfdTAC0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD ISC/HFD TAC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-10,7.5) +
  ylim(0,40) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_ISCvHFD_TAC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::CONTROL ISC vs CONTROL TAC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_conISCVconTAC0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD ISC/HFD TAC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-10,7.5) +
  ylim(0,27.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotControl_ISCvControl_TAC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD TAC vs CONTROL TAC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdTACVconTAC0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD TAC/Control TAC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-10,10) +
  ylim(0,12.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD_TACvControl_TAC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD TAC vs CONTROL TAC (Does not Plot) ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfdTACVconTAC_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD TAC/Control TAC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  # xlim(-10,10) +
  # ylim(0,12.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )

# VolcanoPlot::HFD-1wk vs CONTROL ISC ---------------------------------------
EnhancedVolcano(
  annotatedPeak_hfd1wkVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD-1wk/Control ISC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-6,6) +
  ylim(0,8) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD-1wkvControl_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# VolcanoPlot::HFD-4wk vs CONTROL ISC ---------------------------------------

EnhancedVolcano(
  annotatedPeak_hfd4wkVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (HFD-4wk/Control ISC)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-6,6) +
  ylim(0,8) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotHFD-4wkvControl_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


# VolcanoPlot::HFD-4wk vs CONTROL ISC ---------------------------------------

EnhancedVolcano(
  annotatedPeak_rereversalVconWT0_df,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  pointSize = 5,
  pCutoff = 0.1,
  FCcutoff = 1
) +
  labs(x = "Log2FoldChange (Rereversal/Control)", y = "-Log10p-value(adjusted)", title = "", subtitle = "") +
  xlim(-6,6) +
  ylim(0,12.5) +
  theme_cowplot() +
  background_grid(minor = "xy") +
  theme(
    aspect.ratio = 1 / 1,
    legend.position = "",
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename="VolcanoPlotRerev_vControl_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)




# KEGGPathwayEnrich:: HFD WT ISC vs CON WT ISC ----------------------------
annotatedHF_promoter <-
  annotatedPeak_ISC_hfdVcon00_df[which(
    abs(annotatedPeak_ISC_hfdVcon00_df$distanceToTSS) < 5000 &
      annotatedPeak_ISC_hfdVcon00_df$pvalue < 0.05 &
      annotatedPeak_ISC_hfdVcon00_df$log2FoldChange > 0
  ), ]
keggHF <-
  kegga(
    de = as.character(annotatedHF_promoter$geneId),
    species = "Mm",
    species.KEGG = "mmu",
    plot = T
  )
annotatedCON_promoter <-
  annotatedPeak_ISC_hfdVcon00_df[which(
    abs(annotatedPeak_ISC_hfdVcon00_df$distanceToTSS) < 5000 &
      annotatedPeak_ISC_hfdVcon00_df$pvalue < 0.05 &
      annotatedPeak_ISC_hfdVcon00_df$log2FoldChange < 0
  ), ]
keggCON <-
  kegga(
    de = as.character(annotatedCON_promoter$geneId),
    species = "Mm",
    species.KEGG = "mmu",
    plot = T
  )
keggHF_sig <- keggHF[which(keggHF$P.DE < 0.05), ]
keggCON_sig <- keggCON[which(keggCON$P.DE < 0.05), ]
keggHF_sig$ES <- -log10(keggHF_sig$P.DE)
keggCON_sig$ES <- log10(keggCON_sig$P.DE)

# Find pathways that are unique to df1
pathways_unique_to_HF <-
  setdiff(keggHF_sig$Pathway, keggCON_sig$Pathway)

# Filter df1 to exclude shared pathways
keggHF_sig_filtered <-
  keggHF_sig[keggHF_sig$Pathway %in% pathways_unique_to_HF,]

# Find pathways that are unique to df2
pathways_unique_to_CON <-
  setdiff(keggCON_sig$Pathway, keggHF_sig$Pathway)

# Filter df2 to exclude shared pathways
keggCON_sig_filtered <-
  keggCON_sig[keggCON_sig$Pathway %in% pathways_unique_to_CON,]

kegg_analysis <- rbind(keggHF_sig_filtered, keggCON_sig_filtered)
kegg_analysis_unfiltered <- rbind(keggHF_sig, keggCON_sig)
pathwaystoplot <- c(
  "PPAR signaling pathway",
  "cAMP signaling pathway",
  "cGMP-PKG signaling pathway",
  "Fatty acid degradation",
  "Fatty acid metabolism",
  "ABC transporters",
  "Fanconi anemia pathway",
  "Hippo signaling pathway - multiple species",
  "Maturity onset diabetes of the young",
  "Taurine and hypotaurine metabolism"
)

kegg_analysis$newPathway = str_wrap(kegg_analysis$Pathway, width = 25)

ggplot(data = kegg_analysis[kegg_analysis$Pathway %in% pathwaystoplot,],
       mapping = aes(
         x = ES,
         y = reorder(newPathway, ES),
         fill = ES
       )) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(
    low = "#E69F00",
    mid = "lightblue1",
    high = "#56B4E9",
    midpoint = median.default(kegg_analysis$ES)
  ) +
  labs(title = "Pathway Enrichment Analysis", x = "Enrichment Score", y = "Pathway") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(
    aspect.ratio = 2 / 1,
    legend.position = "right",
    axis.title = element_text(face = "bold", size = 25),
    axis.text.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename = "KEGGenrichmentHFDvCON_ISC.png" ,
        path = "~/Documents/ATAC/crunchtime/",
        dpi = 600)

# KEGGPathwayEnrich:: HFD WT ISC vs CON WT ISC ----------------------------
annotatedHF_promoter <-
  annotatedPeak_ISC_hfdVcon00_df[which(
    abs(annotatedPeak_ISC_hfdVcon00_df$distanceToTSS) < 5000 &
      annotatedPeak_ISC_hfdVcon00_df$pvalue < 0.05 &
      annotatedPeak_ISC_hfdVcon00_df$log2FoldChange > 0
  ), ]
keggHF <-
  kegga(
    de = as.character(annotatedHF_promoter$geneId),
    species = "Mm",
    species.KEGG = "mmu",
    plot = T
  )
annotatedCON_promoter <-
  annotatedPeak_ISC_hfdVcon00_df[which(
    abs(annotatedPeak_ISC_hfdVcon00_df$distanceToTSS) < 5000 &
      annotatedPeak_ISC_hfdVcon00_df$pvalue < 0.05 &
      annotatedPeak_ISC_hfdVcon00_df$log2FoldChange < 0
  ), ]
annotatedCON_enhancer <-
  annotatedPeak_ISC_hfdVcon00_df[which(
    annotatedPeak_ISC_hfdVcon00_df$distanceToTSS < -5000 &
      annotatedPeak_ISC_hfdVcon00_df$pvalue < 0.05 &
      annotatedPeak_ISC_hfdVcon00_df$log2FoldChange < 0
  ), ]
keggCON <-
  kegga(
    de = as.character(annotatedCON_promoter$geneId),
    species = "Mm",
    species.KEGG = "mmu",
    plot = T
  )
keggCONe <-
  kegga(
    de = as.character(annotatedCON_enhancer$geneId),
    species = "Mm",
    species.KEGG = "mmu",
    plot = T
  )
keggHF_sig <- keggHF[which(keggHF$P.DE < 0.05), ]
keggCON_sig <- keggCON[which(keggCON$P.DE < 0.05), ]
keggHF_sig$ES <- -log10(keggHF_sig$P.DE)
keggCON_sig$ES <- log10(keggCON_sig$P.DE)

# Find pathways that are unique to df1
pathways_unique_to_HF <-
  setdiff(keggHF_sig$Pathway, keggCON_sig$Pathway)

# Filter df1 to exclude shared pathways
keggHF_sig_filtered <-
  keggHF_sig[keggHF_sig$Pathway %in% pathways_unique_to_HF,]

# Find pathways that are unique to df2
pathways_unique_to_CON <-
  setdiff(keggCON_sig$Pathway, keggHF_sig$Pathway)

# Filter df2 to exclude shared pathways
keggCON_sig_filtered <-
  keggCON_sig[keggCON_sig$Pathway %in% pathways_unique_to_CON,]

kegg_analysis <- rbind(keggHF_sig_filtered, keggCON_sig_filtered)
kegg_analysis_unfiltered <- rbind(keggHF_sig, keggCON_sig)
pathwaystoplot <- c(
  "PPAR signaling pathway",
  "cAMP signaling pathway",
  "cGMP-PKG signaling pathway",
  "Fatty acid degradation",
  "Fatty acid metabolism",
  "ABC transporters",
  "Fanconi anemia pathway",
  "Hippo signaling pathway - multiple species",
  "Maturity onset diabetes of the young",
  "Taurine and hypotaurine metabolism"
)

kegg_analysis$newPathway = str_wrap(kegg_analysis$Pathway, width = 25)

ggplot(data = kegg_analysis[kegg_analysis$Pathway %in% pathwaystoplot,],
       mapping = aes(
         x = ES,
         y = reorder(newPathway, ES),
         fill = ES
       )) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(
    low = "#E69F00",
    mid = "lightblue1",
    high = "#56B4E9",
    midpoint = median.default(kegg_analysis$ES)
  ) +
  labs(title = "Pathway Enrichment Analysis", x = "Enrichment Score", y = "Pathway") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(
    aspect.ratio = 2 / 1,
    legend.position = "right",
    axis.title = element_text(face = "bold", size = 25),
    axis.text.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.line = element_line(linewidth = 2),
    plot.background = element_rect(fill = "transparent")
  )
ggsave2(filename = "KEGGenrichmentHFDvCON_ISC.png" ,
        path = "~/Documents/ATAC/crunchtime/",
        dpi = 600)



# Reversal analysis -------------------------------------------------------

# Generate regions that are retained after 1 week off of high fat diet
# CON_iscVtac -------------------------------------------------------------

# Subset based on the significance for control ISCs
annotatedPeak_conISCVconTAC_df_sigISCs <-
  annotatedPeak_conISCVconTAC_df[which(
    annotatedPeak_conISCVconTAC_df$padj < 0.2 &
      annotatedPeak_conISCVconTAC_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for control TACs
annotatedPeak_conISCVconTAC_df_sigTACs <-
  annotatedPeak_conISCVconTAC_df[which(
    annotatedPeak_conISCVconTAC_df$padj < 0.2 &
      annotatedPeak_conISCVconTAC_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
CON_iscVtac_regionsISC <-
  GRanges(annotatedPeak_conISCVconTAC_df_sigISCs)
CON_iscVtac_regionsTAC <-
  GRanges(annotatedPeak_conISCVconTAC_df_sigTACs)


# HFD_iscVtac -------------------------------------------------------------

# Subset based on the significance for HFD ISCs
annotatedPeak_hfdISCVhfdTAC_df_sigISCs <-
  annotatedPeak_hfdISCVhfdTAC_df[which(
    annotatedPeak_hfdISCVhfdTAC_df$padj < 0.2 &
      annotatedPeak_hfdISCVhfdTAC_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for HFD TACs
annotatedPeak_hfdISCVhfdTAC_df_sigTACs <-
  annotatedPeak_hfdISCVhfdTAC_df[which(
    annotatedPeak_hfdISCVhfdTAC_df$padj < 0.2 &
      annotatedPeak_hfdISCVhfdTAC_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
HFD_iscVtac_regionsISC <-
  GRanges(annotatedPeak_hfdISCVhfdTAC_df_sigISCs)
HFD_iscVtac_regionsTAC <-
  GRanges(annotatedPeak_hfdISCVhfdTAC_df_sigTACs)

# ISC_hfdVcon -------------------------------------------------------------

# Generate annotated results page
# Subset based on the significance for HFD ISCs
annotatedPeak_ISC_hfdVcon_df_sigHFDs <-
  annotatedPeak_ISC_hfdVcon_df[which(
    annotatedPeak_ISC_hfdVcon_df$padj < 0.2 &
      annotatedPeak_ISC_hfdVcon_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for Control ISCs
annotatedPeak_ISC_hfdVcon_df_sigCONs <-
  annotatedPeak_ISC_hfdVcon_df[which(
    annotatedPeak_ISC_hfdVcon_df$padj < 0.2 &
      annotatedPeak_ISC_hfdVcon_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
ISC_hfdVcon_regionsHFD <-
  GRanges(annotatedPeak_ISC_hfdVcon_df_sigHFDs)
ISC_hfdVcon_regionsCON <-
  GRanges(annotatedPeak_ISC_hfdVcon_df_sigCONs)

# TAC_hfdVcon -------------------------------------------------------------

# Subset based on the significance for HFD TACs
annotatedPeak_hfdTACVconTAC_df_sigHFDs <-
  annotatedPeak_hfdTACVconTAC_df[which(
    annotatedPeak_hfdTACVconTAC_df$padj < 0.2 &
      annotatedPeak_hfdTACVconTAC_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for HFD CONs
annotatedPeak_hfdTACVconTAC_df_sigCONs <-
  annotatedPeak_hfdTACVconTAC_df[which(
    annotatedPeak_hfdTACVconTAC_df$padj < 0.2 &
      annotatedPeak_hfdTACVconTAC_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
TAC_hfdVcon_regionsHFD <-
  GRanges(annotatedPeak_hfdTACVconTAC_df_sigHFDs)
TAC_hfdVcon_regionsCON <-
  GRanges(annotatedPeak_hfdTACVconTAC_df_sigCONs)


# hfd1wkVcon -------------------------------------------------------------

# Subset based on the significance for HFD TACs
annotatedPeak_hfd1wkVcon_df_sigHFD1wk <-
  annotatedPeak_hfd1wkVconWT_df[which(
    annotatedPeak_hfd1wkVconWT_df$padj < 0.2 &
      annotatedPeak_hfd1wkVconWT_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for HFD CONs
annotatedPeak_hfd1wkVcon_df_sigCON <-
  annotatedPeak_hfd1wkVconWT_df[which(
    annotatedPeak_hfd1wkVconWT_df$padj < 0.2 &
      annotatedPeak_hfd1wkVconWT_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
hfd1wkVcon_regionsHFD1wk <-
  GRanges(annotatedPeak_hfd1wkVcon_df_sigHFD1wk)
hfd1wkVcon_regionsCON <-
  GRanges(annotatedPeak_hfd1wkVcon_df_sigCON)


# hfd4wkVcon -------------------------------------------------------------

# Subset based on the significance for HFD TACs
annotatedPeak_hfd4wkVcon_df_sigHFD4wk <-
  annotatedPeak_hfd4wkVconWT_df[which(
    annotatedPeak_hfd4wkVconWT_df$padj < 0.2 &
      annotatedPeak_hfd4wkVconWT_df$log2FoldChange > 0
  ), ]

# Subset based on the significance for HFD CONs
annotatedPeak_hfd4wkVcon_df_sigCON <-
  annotatedPeak_hfd4wkVconWT_df[which(
    annotatedPeak_hfd4wkVconWT_df$padj < 0.2 &
      annotatedPeak_hfd4wkVconWT_df$log2FoldChange < 0
  ), ]

# Generate GRanges object to compare result sections
hfd4wkVcon_regionsHFD4wk <-
  GRanges(annotatedPeak_hfd4wkVcon_df_sigHFD4wk)
hfd4wkVcon_regionsCON <-
  GRanges(annotatedPeak_hfd4wkVcon_df_sigCON)

HFD4wkbed <- annotatedPeak_hfd4wkVcon_df_sigHFD4wk[,c(1:3,6)]
HFD4wkbed$blank <- ""
HFD4wkbed$strand <- "+"
HFD1wkbed <- annotatedPeak_hfd1wkVcon_df_sigHFD1wk[,c(1:3,6)]
HFD1wkbed$blank <- ""
HFD1wkbed$strand <- "+"

write.table(x = HFD4wkbed,
            file = "HFD4wkHOMER.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(x = HFD1wkbed,
            file = "HFD1wkHOMER.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

# Progenitor Analysis -----------------------------------------------------

# Generate any shared regions in ISCs and TACs that change as an effect of diet
common_change_CON <-
  GenomicRanges::intersect(x = ISC_hfdVcon_regionsCON, y = TAC_hfdVcon_regionsCON)
common_change_CON_df <- as.data.frame(common_change_CON)
common_change_HFD <-
  GenomicRanges::intersect(x = ISC_hfdVcon_regionsHFD, y = TAC_hfdVcon_regionsHFD)
common_change_HFD_df <- as.data.frame(common_change_HFD)

# Regions that pass the cutoff of padj < 0.1 are exclusive as an effect of a high-fat diet and all but 1 region in control

# Generate cell-type exclusive differentially accessible regions
ISC_regionsCON <-
  GenomicRanges::setdiff(x = ISC_hfdVcon_regionsCON, y = TAC_hfdVcon_regionsCON)
ISC_regionsCON_df <- as.data.frame(ISC_regionsCON)

TAC_regionsCON <-
  GenomicRanges::setdiff(x = TAC_hfdVcon_regionsCON, y = ISC_hfdVcon_regionsCON)
TAC_regionsCON_df <- as.data.frame(TAC_regionsCON)

# There are cell specific differences as an effect of diet
ISC_regionsHFD <-
  GenomicRanges::setdiff(x = ISC_hfdVcon_regionsHFD, y = TAC_hfdVcon_regionsHFD)
ISC_regionsHFD_df <- as.data.frame(ISC_regionsHFD)

TAC_regionsHFD <-
  GenomicRanges::setdiff(x = TAC_hfdVcon_regionsHFD, y = ISC_hfdVcon_regionsHFD)
TAC_regionsHFD_df <- as.data.frame(TAC_regionsHFD)

# There are cell specific differences as an effect of diet
common_regionsISC <-
  GenomicRanges::intersect(x = HFD_iscVtac_regionsISC, y = CON_iscVtac_regionsISC)
common_regionsISC_df <- as.data.frame(common_regionsISC)

common_regionsTAC <-
  GenomicRanges::intersect(x = HFD_iscVtac_regionsTAC, y = CON_iscVtac_regionsTAC)
common_regionsTAC_df <- as.data.frame(common_regionsTAC)

ISC_regionsCON_df$names <- paste0("peak", 1:nrow(ISC_regionsCON_df))
ISC_regionsCON_df$blank <- ""
ISC_regionsCON_df <- ISC_regionsCON_df[, c(1:3, 6, 7, 5)]
ISC_regionsCON_df$strand <- "+"
write.table(x = ISC_regionsCON_df,
            file = "ISC_regionsCON.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

ISC_regionsHFD_df$names <- paste0("peak", 1:nrow(ISC_regionsHFD_df))
ISC_regionsHFD_df$blank <- ""
ISC_regionsHFD_df <- ISC_regionsHFD_df[, c(1:3, 6, 7, 5)]
ISC_regionsHFD_df$strand <- "+"
write.table(x = ISC_regionsHFD_df,
            file = "ISC_regionsHFD.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

TAC_regionsCON_df$names <- paste0("peak", 1:nrow(TAC_regionsCON_df))
TAC_regionsCON_df$blank <- ""
TAC_regionsCON_df <- TAC_regionsCON_df[, c(1:3, 6, 7, 5)]
TAC_regionsCON_df$strand <- "+"
write.table(x = TAC_regionsCON_df,
            file = "TAC_regionsCON.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

TAC_regionsHFD_df$names <- paste0("peak", 1:nrow(TAC_regionsHFD_df))
TAC_regionsHFD_df$blank <- ""
TAC_regionsHFD_df <- TAC_regionsHFD_df[, c(1:3, 6, 7, 5)]
TAC_regionsHFD_df$strand <- "+"
write.table(x = TAC_regionsHFD_df,
            file = "TAC_regionsHFD.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

plotPCA(object = ntd, intgroup = "group")
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=F)
df <- as.data.frame(colData(dds)[,c("Diet","group")])
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=FALSE,kmeans_k = 10,
         cluster_cols=T, annotation_col=df)

de_genes <- subset(fastedVconWT, padj < 0.001 & abs(log2FoldChange) > 1)

# Extract expression matrix for differentially expressed genes
de_expr <- assay(dds)[rownames(de_genes), ]
de_expr <- as.data.frame(de_expr)
# Group samples based on conditions
conditions <- colData(dds)$group
de_expr_t <- t(de_expr)
# Split expression matrix by conditions
split_expr <- aggregate(de_expr, list(dds$group), FUN = rowMeans)

# Calculate average expression for each group
grouped_expr_means <- lapply(as.array(grouped_expr), conditions, FUN = colMeans, na.rm = TRUE)

# Combine the group means into a matrix
grouped_expr_matrix <- do.call(rbind, grouped_expr_means)

# Generate heatmap
pheatmap(grouped_expr_matrix, cluster_rows = TRUE, cluster_cols = TRUE)
# Determining Specific Regions --------------------------------------------

HFD_isc <-
  setdiff(x = HFD_iscVtac_regionsISC, y = CON_iscVtac_regionsISC)
HFD_isc_df <- as.data.frame(HFD_isc)
CON_isc <-
  setdiff(x = CON_iscVtac_regionsISC, y = HFD_iscVtac_regionsISC)
CON_isc_df <- as.data.frame(CON_isc)
NONspec_isc <-
  intersect(x = HFD_iscVtac_regionsISC, y = CON_iscVtac_regionsISC)
NONspec_isc_df <- as.data.frame(NONspec_isc)

HFD_tac <-
  setdiff(x = HFD_iscVtac_regionsTAC, y = CON_iscVtac_regionsTAC)
HFD_tac_df <- as.data.frame(HFD_tac)
CON_tac <-
  setdiff(x = CON_iscVtac_regionsTAC, y = HFD_iscVtac_regionsTAC)
CON_tac_df <- as.data.frame(CON_tac)
NONspec_tac <-
  intersect(x = HFD_iscVtac_regionsTAC, y = CON_iscVtac_regionsTAC)
NONspec_tac_df <- as.data.frame(NONspec_tac)

HFD_isc_spec <- setdiff(x = HFD_isc, y = CON_isc)
HFD_isc_spec_df <- as.data.frame(HFD_isc_spec)
HFD_isc_spec_df$names <- paste0("peak", 1:nrow(HFD_isc_spec_df))
HFD_isc_spec_df$blank <- ""
HFD_isc_spec_df <- HFD_isc_spec_df[, c(1:3, 6, 7, 5)]
HFD_isc_spec_df$strand <- "+"

# write.table(x = HFD_isc_spec_df,
#             file = "HFD_isc_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

CON_isc_spec <- setdiff(x = CON_isc, y = HFD_isc)
CON_isc_spec_df <- as.data.frame(CON_isc_spec)
CON_isc_spec_df$names <- paste0("peak", 1:nrow(CON_isc_spec_df))
CON_isc_spec_df$blank <- ""
CON_isc_spec_df <- CON_isc_spec_df[, c(1:3, 6, 7, 5)]
CON_isc_spec_df$strand <- "+"

# write.table(x = CON_isc_spec_df,
#             file = "CON_isc_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

HFD_tac_spec <- setdiff(x = HFD_tac, y = CON_tac)
HFD_tac_spec_df <- as.data.frame(HFD_tac_spec)
HFD_tac_spec_df$names <- paste0("peak", 1:nrow(HFD_tac_spec_df))
HFD_tac_spec_df$blank <- ""
HFD_tac_spec_df <- HFD_tac_spec_df[, c(1:3, 6, 7, 5)]
HFD_tac_spec_df$strand <- "+"

# write.table(x = HFD_tac_spec_df,
#             file = "HFD_tac_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

CON_tac_spec <- setdiff(x = CON_tac, y = HFD_tac)
CON_tac_spec_df <- as.data.frame(CON_tac_spec)
CON_tac_spec_df$names <- paste0("peak", 1:nrow(CON_tac_spec_df))
CON_tac_spec_df$blank <- ""
CON_tac_spec_df <- CON_tac_spec_df[, c(1:3, 6, 7, 5)]
CON_tac_spec_df$strand <- "+"

# write.table(x = CON_tac_spec_df,
#             file = "CON_tac_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

consensus_bed <- countMatrix[, c(1:4)]
consensus_bed$blank <- ""
consensus_bed$strand <- "+"

# write.table(x = consensus_bed,
#             file = "consensus.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

TAC_spec <-
  intersect(x = HFD_iscVtac_regionsTAC, y = CON_iscVtac_regionsTAC)
TAC_spec_df <- as.data.frame(TAC_spec)
TAC_spec_df$names <- paste0("peak", 1:nrow(TAC_spec_df))
TAC_spec_df$blank <- ""
TAC_spec_df <- TAC_spec_df[, c(1:3, 6, 7, 5)]
TAC_spec_df$strand <- "+"

# write.table(x = TAC_spec_df,
#             file = "TAC_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

ISC_spec <-
  intersect(x = HFD_iscVtac_regionsISC, y = CON_iscVtac_regionsISC)
ISC_spec_df <- as.data.frame(ISC_spec)
ISC_spec_df$names <- paste0("peak", 1:nrow(ISC_spec_df))
ISC_spec_df$blank <- ""
ISC_spec_df <- ISC_spec_df[, c(1:3, 6, 7, 5)]
ISC_spec_df$strand <- "+"

# write.table(x = ISC_spec_df,
#             file = "ISC_spec_df.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

temp <-
  annotatedPeak_hfdISCVhfdTAC_df[which(annotatedPeak_hfdISCVhfdTAC_df$padj >= 0.05), ]
temp2 <-
  annotatedPeak_conISCVconTAC_df[which(annotatedPeak_conISCVconTAC_df$padj >= 0.05), ]
temp_hfd <- setdiff(x = GRanges(temp), y = GRanges(temp2))
temp_hfd_df <- as.data.frame(temp_hfd)
temp_hfd_df$names <- paste0("peak", 1:nrow(temp_hfd_df))
temp_hfd_df$blank <- ""
temp_hfd_df <- temp_hfd_df[, c(1:3, 6, 7, 5)]
temp_hfd_df$strand <- "+"
temp_con <- setdiff(x = GRanges(temp2), y = GRanges(temp))
temp_con_df <- as.data.frame(temp_con)
temp_con_df$names <- paste0("peak", 1:nrow(temp_con_df))
temp_con_df$blank <- ""
temp_con_df <- temp_con_df[, c(1:3, 6, 7, 5)]
temp_con_df$strand <- "+"
write.table(
  x = temp_hfd_df,
  file = "common_hfd_df.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  x = temp_con_df,
  file = "common_con_df.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

library(msigdbr)
library(clusterProfiler)
library(ReactomePA)

annotatedPeak_ISC_hfdVcon_df_promo <-
  annotatedPeak_ISC_hfdVcon_df[which(abs(annotatedPeak_ISC_hfdVcon_df$distanceToTSS) < 5000),]

control_kegg <-
  enrichKEGG(gene = CON_iscVtac_regionsISC$geneId,
             organism = "mmu",
             keyType = "kegg")
dotplot(control_kegg)
controlTA_kegg <-
  enrichKEGG(gene = CON_iscVtac_regionsTAC$geneId,
             organism = "mmu",
             keyType = "kegg")
dotplot(controlTA_kegg)
HFD_kegg <-
  enrichKEGG(gene = HFD_iscVtac_regionsISC$geneId,
             organism = "mmu",
             keyType = "kegg")
dotplot(HFD_kegg)
HFDTA_kegg <-
  enrichKEGG(gene = HFD_iscVtac_regionsTAC$geneId,
             organism = "mmu",
             keyType = "kegg")
dotplot(HFDTA_kegg)

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# we want the log2 fold change
original_gene_list <- annotatedPeak_ISC_hfdVcon_df_promo$log2FoldChange

# name the vector
names(original_gene_list) <- annotatedPeak_ISC_hfdVcon_df_promo$ENSEMBL

# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = annotatedPeak_ISC_hfdVcon_df_promo[annotatedPeak_ISC_hfdVcon_df_promo$ENSEMBL %in% dedup_ids$ENSEMBL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- annotatedPeak_ISC_hfdVcon_df_promo$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- annotatedPeak_ISC_hfdVcon_df_promo$geneId

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
dedup <- kegg_gene_list[]

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kk2 <- clusterProfiler::gseKEGG(geneList     = kegg_gene_list,
               organism     = "mmu",
               minGSSize    = 1,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")
head(gse$Description)

# Filter values that start with "Gm"
filtered_values <-
  annotatedPeak_hfdDAVconWT_df[grepl("^Gm", annotatedPeak_hfdDAVconWT_df$SYMBOL),]
filtered_values <-
  filtered_values[which(filtered_values$padj < 0.1), ]
filtered_valuesUP <-
  filtered_values[which(filtered_values$log2FoldChange > 0 &
                          abs(filtered_values$distanceToTSS) < 5000), ]
filtered_valuesDOWN <-
  filtered_values[which(filtered_values$log2FoldChange < 0 &
                          abs(filtered_values$distanceToTSS) < 5000), ]

# Print the filtered values
print(filtered_valuesUP$SYMBOL)
# Making PCA --------------------------------------------------------------

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$padj < 0.25),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1,3,5,7,9,11,13,15,17:22)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$Diet[c(1,3,5,7,9,11,13,15,17:22)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = `meta_data$Diet[c(1, 3, 5, 7, 9, 11, 13, 15, 17:22)]`)) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = `meta_data$Diet[c(1, 3, 5, 7, 9, 11, 13, 15, 17:22)]`)) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD ISC vs Control ISC", x = "PC1 (45.69%)", y = "PC2 (12.88%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFDvCON_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

# Generate data frame with normalized counts
df <- as.data.frame(counts(dds0, normalized = TRUE))

# Generate a list of peak names that are in the top 25% of most comparable regions
x <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$padj < 0.25),]
df <- df[rownames(df) %in% x$name,]

# Subset to just ISC WT samples
df <- df[,c(1:22)]

# Perform log2 transformation
log.pl <- log2(df + 1)

# Transpose dataframe
t.log.pl <- t(log.pl)

# Perform principle component analysis
prin_comp <- prcomp(t.log.pl, rank. = 4)

# Assign components
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components$PC4 <- -components$PC4
components = cbind(components, meta_data$group[c(1:22)])

# 2D PCA
summary(prin_comp) # Run to adjust % variation per principle component

ggplot(data = components, mapping = aes(x = PC1, y = PC2, color = `meta_data$group[c(1:22)]`)) +
  #geom_text(label = rownames(components)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = `meta_data$group[c(1:22)]`)) +  # Add geometric ellipses with fill colors
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values = c("#E69F00", "#DB5A00", "#56B4E9", "#54D98C")) +
  scale_fill_manual(values = c("#E69F00", "#DB5A00", "#56B4E9", "#54D98C")) +  # Set fill colors manually
  # xlim(-20, 20) + ylim(-5, 10) +
  labs(title = "PCA", subtitle = "HFD & Control ISC & TAC", x = "PC1 (32.18%)", y = "PC2 (18.30%)") +
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold",size = 25), 
        axis.text = element_text(face = "bold",size = 20), axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent"))
ggsave2(filename="PCA_top25percent_HFDvCON_ISCandTAC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)


a <- annotatedPeak_ISC_hfdVcon0_df[which(annotatedPeak_ISC_hfdVcon0_df$padj < 0.1 & annotatedPeak_ISC_hfdVcon0_df$log2FoldChange < 0),]
b <- annotatedPeak_hfd1wkVconWT0_df[which(annotatedPeak_hfd1wkVconWT0_df$padj < 0.1 & annotatedPeak_hfd1wkVconWT0_df$log2FoldChange < 0),]
c <- intersect(x = GRanges(a), y = GRanges(b))
c_df <- as.data.frame(c)
d <- annotatedPeak_hfd4wkVconWT0_df[which(annotatedPeak_hfd4wkVconWT0_df$padj < 0.1 & annotatedPeak_hfd4wkVconWT0_df$log2FoldChange < 0),]
e <- intersect(x = GRanges(a), y = GRanges(d))
e_df <- as.data.frame(e)


hfdr <- annotatedPeak_ISC_hfdVcon00_df[,c(6,17,20,24)]
hfdr <- cbind(hfdr,annotatedPeak_hfd1wkVconWT00_df[,c(20,24)])
hfdr <- cbind(hfdr,annotatedPeak_hfd4wkVconWT00_df[,c(20,24)])
colnames(hfdr) <- c("name","SYMBOL","FC0","P0","FC1","P1","FC4","P4")
hfdpeaknames <- annotatedPeak_ISC_hfdVcon00_df[which(annotatedPeak_ISC_hfdVcon00_df$padj < 0.1 & annotatedPeak_ISC_hfdVcon00_df$log2FoldChange > 0),]
hfdr_filt <- hfdr %>% filter(name %in% hfdpeaknames$name)
hfdr_filt0 <- hfdr_filt
hfdr_filt0$group <- ifelse(hfdr_filt0$P4 < 0.1, "LT",
                           ifelse(hfdr_filt0$P1 < 0.1, "ST", "NM"))
# Reshape the data frame from wide to long format
df_long <- tidyr::pivot_longer(
  hfdr_filt0,
  cols = c(P0, P1, P4),
  names_to = "Timepoint",
  values_to = "P"
)
ggplot(df_long, aes(x = Timepoint, y = -log10(P), group = name)) +
  geom_smooth(method = "loess", se = F, aes(color = group)) +
  geom_point(aes(color = group), size = 3) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", linewidth = 2) +
  labs(x = "Timepoint", y = "-Log10p-value(adjusted)", title = "Significance at Different Timepoints") +
  scale_color_manual(values = c("ST" = "gold", "LT" = "green3", "NM" = "red3")) + 
  theme_cowplot() + background_grid(minor = "xy") +
  theme(aspect.ratio = 1/1, legend.position = "",
        axis.title = element_text(face = "bold", size = 25),
        axis.text = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 2),
        plot.background = element_rect(fill = "transparent")) +
  scale_x_discrete(expand = expansion(add = c(0.175, 0.175)),labels=c("P0" = "HFD", "P1" = "HFD-1wk", "P4" = "HFD-4wk"))
ggsave2(filename="Memory_at1wk&4wk_ISC.png" ,path = "~/Documents/ATAC/crunchtime/", dpi = 600)

LT <- hfdr_filt0 %>% filter(group == "LT")
ST <- hfdr_filt0 %>% filter(group == "ST")
NM <- hfdr_filt0 %>% filter(group == "NM")

hfd_lt <- annotatedPeak_ISC_hfdVcon00_df %>% filter(name %in% LT$name)
hfd_lt <- hfd_lt[,c(1:3,6)]
hfd_lt$blank <- ""
hfd_lt$strand <- "+"
write.table(x = hfd_lt,
            file = "hfd_lt.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

hfd_st <- annotatedPeak_ISC_hfdVcon00_df %>% filter(name %in% ST$name)
hfd_st <- hfd_st[,c(1:3,6)]
hfd_st$blank <- ""
hfd_st$strand <- "+"
write.table(x = hfd_st,
            file = "hfd_st.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
hfd_nm <- annotatedPeak_ISC_hfdVcon00_df %>% filter(name %in% NM$name)
hfd_nm <- hfd_nm[,c(1:3,6)]
hfd_nm$blank <- ""
hfd_nm$strand <- "+"
write.table(x = hfd_nm,
            file = "hfd_nm.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

x <- annotatedPeak_ISC_hfdVcon00_df %>% filter(padj < 0.1)
c <- x %>% filter(log2FoldChange < 0)
h <- x %>% filter(log2FoldChange > 0)
h0 <- h %>% filter(distanceToTSS < 1000 & distanceToTSS > -5000)

x1 <- annotatedPeak_hfd1wkVconWT00_df %>% filter(pvalue < 0.05)
x1 <- x1 %>% filter(name %in% x$name)
c1 <- x1 %>% filter(log2FoldChange < 0)
h1 <- x1 %>% filter(log2FoldChange > 0)

h <- h[,c(1:3,6)]
h$blank <- ""
h$strand <- "+"
h0 <- h0[,c(1:3,6)]
h0$blank <- ""
h0$strand <- "+"

write.table(x = h,
            file = "hfdpadj05.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(x = h0,
            file = "hfdpadj05_cisreg.bed",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

# annotatedPeak_ISC_hfdVcon00_df <- annotatedPeak_ISC_hfdVcon00_df[,c(1:3,6)]
# annotatedPeak_ISC_hfdVcon00_df$blank <- ""
# annotatedPeak_ISC_hfdVcon00_df$strand <- "+"
# 
# write.table(x = annotatedPeak_ISC_hfdVcon00_df,
#             file = "conc.bed",
#             quote = FALSE,
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)
ppartarget <- read.delim("~/Dropbox (ASU)/Mac (2)/Downloads/ppartarget.txt")
`pparDtarget.(1)` <- read.delim("~/Dropbox (ASU)/Mac (2)/Downloads/pparDtarget (1).txt")
target <- cbind(ppartarget,`pparDtarget.(1)`[,22])
x <- annotatedPeak_ISC_hfdVcon00_df %>% filter(name %in% target$PeakID)
ppartarget_sorted <- target[order(target$PeakID),]
x_sorted <- x[order(x$name),]
x_fin <- cbind(x_sorted, ppartarget_sorted[,c(20:25)])

write.table(x = x_fin,
            file = "pparTargets.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = T)
