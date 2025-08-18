library(fishpond)
library(readr)
library(tximport)
library(tximeta)
library(SummarizedExperiment)
library(tidyverse)

# read sample CSV
Taf1 <- read_csv("/Users/UCL-TAF1/code/coldata_tximeta_mouse1.csv",
                 col_types = cols())
# Add row name
Taf1 <- as.data.frame(Taf1)
rownames(Taf1) <- Taf1$sample
Taf1$names = Taf1$sample

# Use tximeta build SummarizedExperiment dataset
se_all <- tximeta(coldata = Taf1, type = 'oarfish')

# Scale inference replicates for normalization
norm_exp <- scaleInfReps(se_all,saveMeanScaled = TRUE)
# Remove transcripts with zero expression
norm_exp <- labelKeep(norm_exp)
norm_exp <- norm_exp[mcols(norm_exp)$keep,]

# Assign gene_id manually (all transcripts belong to Taf1)
mcols(norm_exp)$gene_id <- rep('Taf1', nrow(mcols(norm_exp)))

# Add data for Swish
# Convert factors and ensure proper data types
norm_exp$tissue <- as.factor(norm_exp$tissue)
norm_exp$sex <- as.factor(norm_exp$sex)

# Set seed for reproducibility
set.seed(120)
normalized_everything = assays(norm_exp)$meanScaled %>% as.data.frame()

# Generate DTE heatmap across all tissues(aggregated over sex)
normalized_everything %>%
  rownames_to_column('txp') %>% 
  reshape2::melt() %>% 
  left_join(Taf1, by = c("variable" = 'names')) %>% 
  mutate(txp = gsub("\\|.*","",txp)) %>% 
  group_by(tissue,txp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ggplot(aes(x = tissue, y = txp, fill = mean_exp)) + 
  geom_tile()

# Create transcript name mapping
transcript_mapping <- c(
  "ENSMUST00000249914.1" = "ENSMUST00000249914.1(Taf1-207)",
  "NCBI_NM_001290729.3_nTaf1_7168bp" = "NM_001290729.3(NCBI_nTaf1)",
  "ENSMUST00000143908.8" = "ENSMUST00000143908.8(Taf1-204)",
  "NM_001405959.2_Mus" = "NM_001405959.2(NCBI_cTaf1)"
)

# Create transcript ordering
# MSTRG transcripts in numerical order (MSTRG.1.1 to MSTRG.1.19)
mstrg_order <- paste0("MSTRG.1.", 1:19)
# Known transcripts in specified order
known_order <- c("NM_001405959.2(NCBI_cTaf1)", 
                 "NM_001290729.3(NCBI_nTaf1)", 
                 "ENSMUST00000249914.1(Taf1-207)", 
                 "ENSMUST00000143908.8(Taf1-204)")
# Combine orders
transcript_order <- c(mstrg_order, known_order)
desired_order <- c("striatum", "cerebellum", "hippocampus", "cortex", "heart", "muscle", "spleen")

# Prepare DTE data for heatmap construction
final_plot <- normalized_everything %>%
  rownames_to_column('txp') %>% 
  reshape2::melt() %>% 
  left_join(Taf1, by = c("variable" = 'names')) %>% 
  mutate(txp = gsub("\\|.*","",txp)) %>%
  
  # Filter out unwanted transcript
  filter(txp != "ENSMUST00000118878.9") %>%
  # Apply transcript name mapping
  mutate(txp = ifelse(txp %in% names(transcript_mapping), 
                      transcript_mapping[txp], 
                      txp)) %>%
  
  mutate(txp = factor(txp, levels = transcript_order)) %>%
  mutate(tissue = factor(tissue, levels = desired_order)) %>%
  group_by(tissue, txp, sex) %>% 
  summarise(mean_exp = mean(value, na.rm = TRUE)) %>%
  # Creat heatmap
  ggplot(aes(x = tissue, y = txp, fill = mean_exp)) + 
  geom_tile(color = NA) +
  facet_wrap(~sex) +
  scale_fill_gradientn(
    colors = c("#F5F5F5", "#27408B"),
    values = scales::rescale(c(0, 1)),
    na.value = "grey90",
    name = "Expression Level"
  ) +
  labs(x = "Tissue", y = "Transcript") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 15))

library(svglite)
svglite("/Users/lipeihang/Desktop/PaperFigure/heat_dte.svg", width = 8, height = 6, bg = "white")
print(final_plot)
dev.off()

