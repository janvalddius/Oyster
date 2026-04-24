
### DATA CLEAN OF OUTLIERS  ----
countData            <- read_excel("D:/Decicomp/Decicomp_matrix.xlsx", sheet = "matrix_outliers")
countData            <- as.data.frame(countData)
row.names(countData) <- countData[, 1]
countData            <- countData[,-1]
countData            <- countData %>% #dplyr::select(matches(grep("T0|T0|T6|T12", colnames(countData), value = TRUE))) %>% 
  dplyr::select(-c("F14R_C3_T24_S5", 
                   "F14R_C3_T24_S6",
                   
                   "F14R_C7_T0_S3",  
                   "F14R_C7_T24_S4",
                   
                   "F14R_C8_T6_S3",  
                   "F14R_C8_T24_S1",
                   
                   "H2D_C3_T24_S4")) 

coldata <- colnames(countData) %>% 
  as.data.frame() %>%
  dplyr::rename("sample"=".") %>%
  dplyr::mutate(group = case_when(
    grepl("H2D_C3_T0",  sample) ~  "H2D_C3_T0",
    grepl("H2D_C3_T3",  sample) ~  "H2D_C3_T3",
    grepl("H2D_C3_T6",  sample) ~  "H2D_C3_T6",
    grepl("H2D_C3_T12", sample) ~  "H2D_C3_T12",
    grepl("H2D_C3_T24", sample) ~  "H2D_C3_T24",
    
    grepl("H2D_C7_T0",  sample) ~  "H2D_C7_T0",
    grepl("H2D_C7_T3",  sample) ~  "H2D_C7_T3",
    grepl("H2D_C7_T6",  sample) ~  "H2D_C7_T6",
    grepl("H2D_C7_T12", sample) ~  "H2D_C7_T12",
    grepl("H2D_C7_T24", sample) ~  "H2D_C7_T24",
    
    grepl("H2D_C8_T0",  sample) ~  "H2D_C8_T0", 
    grepl("H2D_C8_T3",  sample) ~  "H2D_C8_T3",
    grepl("H2D_C8_T6",  sample) ~  "H2D_C8_T6",
    grepl("H2D_C8_T12", sample) ~  "H2D_C8_T12",
    grepl("H2D_C8_T24", sample) ~  "H2D_C8_T24",
    
    grepl("F14R_C3_T0", sample) ~  "F14R_C3_T0",
    grepl("F14R_C3_T3", sample) ~  "F14R_C3_T3",
    grepl("F14R_C3_T6", sample) ~  "F14R_C3_T6",
    grepl("F14R_C3_T12",sample) ~  "F14R_C3_T12",
    grepl("F14R_C3_T24",sample) ~  "F14R_C3_T24",
    
    grepl("F14R_C7_T0", sample) ~  "F14R_C7_T0",
    grepl("F14R_C7_T3", sample) ~  "F14R_C7_T3",
    grepl("F14R_C7_T6", sample) ~  "F14R_C7_T6",
    grepl("F14R_C7_T12",sample) ~  "F14R_C7_T12",
    grepl("F14R_C7_T24",sample) ~  "F14R_C7_T24",
    
    grepl("F14R_C8_T0", sample) ~  "F14R_C8_T0", 
    grepl("F14R_C8_T3", sample) ~  "F14R_C8_T3",
    grepl("F14R_C8_T6", sample) ~  "F14R_C8_T6",
    grepl("F14R_C8_T12",sample) ~  "F14R_C8_T12", 
    grepl("F14R_C8_T24",sample) ~  "F14R_C8_T24", TRUE ~ NA_character_)) %>% 
  tibble::column_to_rownames(var = "sample")

coldata$group        <- as.factor(coldata$group)

dds                  <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ group)
keep                 <- rowSums(counts(dds)) >= 15 # dim
dds                  <- dds[keep,]
dds                  <- DESeq(dds)
sizeFactors(dds)
resultsNames(dds) # importante
normalizedCounts    <- counts(dds, normalized = TRUE) # dim(normalizedCounts)
sampleData          <- colData(dds)
dds_norm            <- vst(dds, blind = F) # Blind =F, to take into account the experimental design

norm_genes_T0 <- as.data.frame(assay(dds_norm)) %>% tibble::rownames_to_column("Gene")

#
F14R_C3_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S1)
write.table(F14R_C3_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C3_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S2)
write.table(F14R_C3_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C3_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S3)
write.table(F14R_C3_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C3_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S4)
write.table(F14R_C3_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C3_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S5)
write.table(F14R_C3_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C3_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C3_T0_S6)
write.table(F14R_C3_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C3_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

#
F14R_C7_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S1)
write.table(F14R_C7_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C7_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S2)
write.table(F14R_C7_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

#F14R_C7_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S3)
#write.table(F14R_C7_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C7_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S4)
write.table(F14R_C7_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C7_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S5)
write.table(F14R_C7_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C7_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C7_T0_S6)
write.table(F14R_C7_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C7_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

#
F14R_C8_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S1)
write.table(F14R_C8_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C8_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S2)
write.table(F14R_C8_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C8_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S3)
write.table(F14R_C8_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C8_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S4)
write.table(F14R_C8_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C8_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S5)
write.table(F14R_C8_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

F14R_C8_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, F14R_C8_T0_S6)
write.table(F14R_C8_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/F14R_C8_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")


#
H2D_C3_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S1)
write.table(H2D_C3_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C3_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S2)
write.table(H2D_C3_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C3_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S3)
write.table(H2D_C3_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C3_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S4)
write.table(H2D_C3_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C3_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S5)
write.table(H2D_C3_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C3_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C3_T0_S6)
write.table(H2D_C3_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C3_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

#
H2D_C7_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S1)
write.table(H2D_C7_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C7_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S2)
write.table(H2D_C7_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C7_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S3)
write.table(H2D_C7_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C7_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S4)
write.table(H2D_C7_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C7_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S5)
write.table(H2D_C7_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C7_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C7_T0_S6)
write.table(H2D_C7_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C7_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

#
H2D_C8_T0_S1 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S1)
write.table(H2D_C8_T0_S1, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S1.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C8_T0_S2 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S2)
write.table(H2D_C8_T0_S2, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S2.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C8_T0_S3 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S3)
write.table(H2D_C8_T0_S3, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S3.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C8_T0_S4 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S4)
write.table(H2D_C8_T0_S4, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S4.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C8_T0_S5 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S5)
write.table(H2D_C8_T0_S5, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S5.tsv", row.names = F, quote = F, col.names=F, sep = "\t")

H2D_C8_T0_S6 <- norm_genes_T0 %>% dplyr::select(Gene, H2D_C8_T0_S6)
write.table(H2D_C8_T0_S6, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/Celine/H2D_C8_T0_S6.tsv", row.names = F, quote = F, col.names=F, sep = "\t")




