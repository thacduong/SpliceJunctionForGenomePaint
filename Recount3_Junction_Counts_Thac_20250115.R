##### Install and Load Reqired Libraries #####
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("recount3", "GenomicRanges", "org.Hs.eg.db", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(recount3)
library(GenomicRanges)
library(org.Hs.eg.db)
library(dplyr)

##### Step 1: Set Working Directory #####
work_dir = "C:/Users/80026261/OneDrive - Moffitt Cancer Center/Desktop/Perna Leukemia Proj/Recount3_SpliceJxn_workflow_SLC12A7_Thac"
setwd(work_dir)

##### Step 2: Define Gene Region and Project IDs #####
# Region of Interest
splice_junction = "SLC12A7"
gene_chr = "5" # SLC12A7
gene_start = 1081895
gene_end = 1115816

project = "Leucegene"
project_ids = c(
  "SRP056295",
  "SRP048759",
  "SRP033266",
  "SRP028567",
  "SRP056197",
  "SRP032455",
  "SRP122507",
  "SRP028594",
  "SRP098651",
  "SRP028554"
) # Leucegene

# project = "TCGA_AML"
# project_ids = c(
#   "LAML" # TCGA_AML
# )

##### Step 3: Define Function to Analyze a Project #####
analyze_project <- function(project_id, gene_chr, gene_start, gene_end, output_file) {
  # Find available projects in recount3
  available_projects <- available_projects()
  
  # Filter for the desired project
  project_info <- subset(available_projects, project == project_id)
  
  if (nrow(project_info) == 0) {
    stop(paste("Project", project_id, "not found in recount3 database."))
  }
  
  # Create the RangedSummarizedExperiment (RSE) for junction-level data
  rse_jx <- create_rse(project_info, type = "jxn")
  
  # Extract junction ranges
  jx_gr <- rowRanges(rse_jx)
  
  # Extract motifs: Left and Right
  motif_left <- as.character(mcols(jx_gr)$left_motif)
  motif_right <- as.character(mcols(jx_gr)$right_motif)
  print(motif_right)
  print(motif_left)
  
  # Define genomic region for the gene
  gene_region <- GRanges(
    seqnames = paste0("chr", gene_chr),  # Ensure 'chr' prefix matches the dataset
    ranges = IRanges(
      start = gene_start,
      end = gene_end
    )
  )
  
  # Filter junctions overlapping the gene region
  gene_junctions <- jx_gr[overlapsAny(jx_gr, gene_region)]
  
  if (length(gene_junctions) == 0) {
    warning(paste("No junctions found for project", project_id, "in the specified region."))
    return(NULL)
  }
  
  # Subset the RSE object to only include filtered junctions
  rse_filtered <- rse_jx[overlapsAny(jx_gr, gene_region), ]
  
  # Calculate average junction counts across samples
  junction_counts <- assays(rse_filtered)$counts
  junction_sum <- as.integer(rowSums(junction_counts))
  # junction_averages <- rowMeans(junction_counts)
  
  # Combine junction coordinates with counts
  junction_table <- data.frame(
    Chr = seqnames(gene_junctions),
    Start = start(gene_junctions),
    End = end(gene_junctions),
    Strand = strand(gene_junctions),
    Motif_Left = motif_left[overlapsAny(jx_gr, gene_region)],
    Motif_Right = motif_right[overlapsAny(jx_gr, gene_region)],
    Sum_Count = junction_sum # sum of all counts
  )
  View(junction_table)
  
  # Save results to file
  write.csv(junction_table, output_file, row.names = FALSE)
  
  cat(paste("Results for project", project_id, "_",  splice_junction, "saved to", output_file, "\n"))
}

##### Step 4: Analyse Each Project #####
output_dir <- paste0(work_dir, "/", project, "/Junction_Counts_TD")
if (!dir.exists(output_dir)) dir.create(output_dir)

for (project_id in project_ids) {
  output_file <- file.path(output_dir, paste0(project_id, "_", splice_junction, "_junction_counts.csv"))
  analyze_project(project_id, gene_chr, gene_start, gene_end, output_file)
}

cat("Analysis completed for all projects.\n")

##### Step 5: Convert Junction Counts to BED Files #####
convert_to_bed <- function(jxn_count, bed_file, include_both_motifs = TRUE) {
  # Load the junction table
  junction_table <- read.csv(jxn_count)
  
  # Create the 'Name' field for BED
  if (include_both_motifs) {
    junction_table$Name <- paste(junction_table$Motif_Left, 
                                 junction_table$Motif_Right, 
                                 sep = "-")
  } else {
    # Use only the left motif
    junction_table$Name <- junction_table$Motif_Left
  }
  View(junction_table$Name)
  
  # Convert to BED format
  bed_data <- data.frame(
    Chr = junction_table$Chr,
    Start = junction_table$Start - 1,  # Convert to 0-based start for BED format
    End = junction_table$End,         # 1-based end remains as is
    Strand = junction_table$Strand,
    Name = junction_table$Name,
    Score = junction_table$Sum_Count  # Use average count as the score
  )
  View(bed_data)
  
  # Write BED file
  write.table(
    bed_data,
    bed_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  cat(paste("BED file created:", bed_file, "\n"))
}

##### Step 6: Summarize Junction Expression Across Samples #####
summarize_junction_expression <- function(jxn_count, cutoff = 3) {
  # Load the junction table
  junction_table <- read.csv(jxn_count)
  
  # Count how many samples express each junction above the cutoff
  junction_counts <- as.matrix(junction_table$Sum_Count)  # Exclude Chr, Start, End, Strand
  print(junction_counts)
  expressing_samples <- rowSums(junction_counts >= cutoff)
  
  # Add expression summary to the table
  junction_table$Expressing_Samples <- expressing_samples
  
  # Save the updated table
  summary_file <- sub(".csv", "_summary.csv", jxn_count)
  write.csv(junction_table, summary_file, row.names = FALSE)
  
  cat(paste("Summary saved:", summary_file, "\n"))
}

##### Step 6: Process All Projects #####
output_bed_dir <- file.path(output_dir, "BED_Files")
if (!dir.exists(output_bed_dir)) dir.create(output_bed_dir)

for (project_id in project_ids) {
  # Input CSV and output BED paths
  jxn_count <- file.path(output_dir, paste0(project_id, "_", splice_junction, "_junction_counts.csv"))
  bed_file <- file.path(output_bed_dir, paste0(project_id, "_", splice_junction, "_junctions.bed"))
  
  # Convert to BED format
  convert_to_bed(jxn_count, bed_file)
  
  # Summarize junction expression
  summarize_junction_expression(jxn_count, cutoff = 3)
}

cat("Conversion to BED and summarization completed for all projects.\n")

