# Load necessary libraries
library(minfi)  # For methylation data analysis
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # Annotation for EPIC array
library(IlluminaHumanMethylationEPICmanifest)  # Manifest for EPIC array
library(FlowSorted.Blood.EPIC)  # Reference for blood cell types
library(GeneNet)  # For network analysis

# Step 1: List IDAT files and create a sample sheet
# ------------------------------------------------
# List all IDAT files in the current directory
idat_files <- list.files(pattern = "idat")

# Extract basenames by removing "_Grn.idat" or "_Red.idat" suffixes
basenames <- unique(gsub("_(Grn|Red).idat", "", idat_files))

# Create a targets data frame with basenames and sample IDs
targets <- data.frame(Basename = basenames, ID = basenames)

# Save the targets data frame to a CSV file
write.csv(targets, "clean_targets.csv", row.names = FALSE)

# Step 2: Read the sample sheet and load raw data
# ------------------------------------------------
# Read the sample sheet from the current directory
targets <- read.csv("clean_targets.csv")

# Load raw methylation data from IDAT files into an RGChannelSet object
rgset <- read.metharray.exp(targets = targets, force = TRUE)

# Step 3: Quality control - Filter samples and probes based on detection p-values
# ------------------------------------------------------------------------------
# Calculate detection p-values for all probes and samples
detP <- detectionP(rgset)

# Keep samples with mean detection p-value < 0.05
keep_samples <- colMeans(detP) < 0.05

# Keep probes with mean detection p-value < 0.05
keep_probes <- rowMeans(detP) < 0.05

# Filter the RGChannelSet to retain only high-quality samples and probes
rgset_filtered <- rgset[keep_probes, keep_samples]

# Step 4: Normalize the data using functional normalization (Funnorm)
# -------------------------------------------------------------------
# Apply functional normalization to the filtered data
rgset_processed <- preprocessFunnorm(rgset_filtered)

# Step 5: Predict and plot sample sex information
# -----------------------------------------------
# Predict sample sex using methylation data
sex_info <- getSex(rgset_processed, cutoff = -2)

# Add sex information to the sample metadata
colData(rgset_processed)$predicted_sex <- sex_info$predictedSex
# Plot sex prediction results
plotSex(rgset_processed)

# Step 6: Map to genome and filter SNP-affected probes
# ---------------------------------------------------
# Map the normalized data to the genome to create a GenomicRatioSet
GRset <- mapToGenome(rgset_processed)

# Remove probes affected by SNPs (SBE and CpG sites with MAF > 0)
GRset <- dropLociWithSnps(GRset, snps = c("SBE", "CpG"), maf = 0)

# Step 7: Remove cross-reactive probes
# ------------------------------------
# Load a list of cross-reactive probes from a CSV file
cross_reactive <- read.csv("cross_reactive_probes.csv")$probeID

# Filter out cross-reactive probes from the GenomicRatioSet
GRset <- GRset[!rownames(GRset) %in% cross_reactive, ]

# Step 8: Extract SNP information, annotation, and beta values
# ------------------------------------------------------------
# Extract SNP information from the GenomicRatioSet
snpinfo <- getSnpInfo(GRset)

# Extract annotation information for the probes
Annotation <- getAnnotation(GRset)

# Extract beta values (methylation levels) from the GenomicRatioSet
Beta2 <- getBeta(GRset)

# Step 9: Save key objects for downstream analysis
# -----------------------------------------------
# Save SNP information to an RData file
save(snpinfo, file = "snpinfo_56789_full.rdata")

# Save annotation information to an RData file
save(Annotation, file = "Annotation_56789_full.rdata")

# Save beta values to an RData file
save(Beta2, file = "Beta_56789_full.rdata")

# Save the GenomicRatioSet object to an RData file
save(GRset, file = "GRset_56789_full.rdata")

# Step 10: Perform quality control (QC) on the data
# -------------------------------------------------
# Convert the RGChannelSet to a MethylSet using raw intensities
mset_filtered <- preprocessRaw(rgset_filtered)

# Map the MethylSet to the genome to create a GenomicMethylSet
gmset_filtered <- mapToGenome(mset_filtered)

# Perform QC on the GenomicMethylSet
qc_res <- getQC(gmset_filtered)

# Save QC results to a CSV file
write.csv(qc_res, "qc_res.csv")

# Plot QC results
plotQC(qc_res)

# Generate a QC report in PDF format
qcReport(rgset, pdf = "qcreport.pdf")

# Step 11: Network analysis using GeneNet
# ---------------------------------------
# Define a function to perform Gaussian Graphical Model (GGM) analysis
Get_GGM <- function(data, ggmfilename) {
  # Estimate partial correlations between probes
  pcorAll2 <- ggm.estimate.pcor(t(data))  # Transpose for features x samples
  
  # Test the significance of edges in the network
  pvalsAll2 <- network.test.edges(pcorAll2, direct = TRUE)
  
  # Apply Bonferroni correction for multiple testing
  benf <- 0.05 / (ncol(data) * (ncol(data) - 1) / 2)  # Bonferroni threshold
  
  # Keep only significant edges (p-value < Bonferroni threshold)
  AllPairsSig <- pvalsAll2[pvalsAll2[, 4] < benf, ]
  
  # Save the significant edges to an RData file
  save(AllPairsSig, file = ggmfilename)
  
  # Return the significant edges
  return(AllPairsSig)
}

# Step 12: Subset beta values for network analysis
# ------------------------------------------------
# Calculate the variance of beta values for each probe
probe_variance <- rowVars(Beta2)

# Select the top 10,000 most variable probes
top_probes <- order(probe_variance, decreasing = TRUE)[1:10000]

# Subset the beta values to include only the top probes
Beta2_subset <- Beta2[top_probes, ]
probe_variance_2 <- rowVars(Beta2_subset)

write.csv(probe_variance_2,"probe_variance_2.csv")  # to be entered as attribute in cytoscape (import file form table)
# Step 13: Perform GGM analysis on the subset of beta values
# ---------------------------------------------------------
# Run the GGM function on the subset of beta values
Pairs <- Get_GGM(Beta2_subset, "methylation_network_edges.rdata")

#OR (for faster work):
#load("methylation_network_edges.rdata")
#Pairs=AllPairsSig

# Step 14: Export network results for visualization in Cytoscape
# -------------------------------------------------------------
# Create a data frame of significant edges for Cytoscape
network_table <- data.frame(
  Source = rownames(Beta2_subset)[Pairs[, 2]],  # Source node (probe)
  Target = rownames(Beta2_subset)[Pairs[, 3]],  # Target node (probe)
  Weight = Pairs[, 1],  # Partial correlation value
  P.Value = Pairs[, 4]  # P-value for the edge
)
#dim(Pairs)
#max_pval=max(network_table$P.Value)  #0.00026315
network_table_more_strc=network_table[which(network_table$P.Value<0.000000001), ]
# Save the network table to a CSV file for Cytoscape
write.csv(network_table, "cytoscape_network.csv", row.names = FALSE)

write.csv(network_table_more_strc, "network_table_more_strc.csv", row.names = FALSE)  #import network form file in cytoscape
