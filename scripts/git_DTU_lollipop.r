library(fishpond)
library(readr)
library(tximport)
library(tximeta)
library(SummarizedExperiment)
library(tidyverse)
library(Cairo)

## Read sample information and construct SE
Taf1 <- readr::read_csv(
  
  "/Users/UCL-TAF1/code/coldata_tximeta_mouse1.csv",
  col_types = cols()
)
rownames(Taf1) <- Taf1$sample  #set row name
Taf1$names       <- Taf1$sample      
se_all <- tximeta(Taf1, type = "oarfish")

## Normalization and filtering low expression
y0 <- scaleInfReps(se_all, saveMeanScaled = TRUE)
y0 <- labelKeep(y0)
y0 <- y0[mcols(y0)$keep, ]

## Analyze only Taf1 gene
mcols(y0)$gene_id <- rep("Taf1", nrow(y0))

## Set covariates as factors
colData(y0)$tissue <- factor(colData(y0)$tissue,
                             levels = c("striatum", "hippocampus", "cortex", "cerebellum", "heart"))
colData(y0)$sex    <- factor(colData(y0)$sex)

run_dtu_pair <- function(se_obj, t1, t2, cov = "sex",
                         q_cut = 0.10, seed = 120) {
## Take subset
  keep_cols <- colData(se_obj)$tissue %in% c(t1, t2)
  y <- se_obj[, keep_cols]
# Downgrade to prevent Swish from warning about factors containing unused levels
  colData(y)$tissue <- droplevels(colData(y)$tissue)
  colData(y)$sex    <- droplevels(colData(y)$sex)
  
## DTE Swish analysis
  set.seed(seed)
  y <- swish(y, x = "tissue", cov = cov)
  
## DTU Swish analysis
  iso <- isoformProportions(y)
  iso <- swish(iso, x = "tissue", cov = cov)
  
  
## Organize output;
# only DTU results are retured here；to obtain DTE results substitute mcols(iso) to mcols(y)
  dtu_res <- mcols(iso) |>
    tibble::as_tibble() |>
    mutate(
      transcript = rownames(assay(iso)),
      transcript = sub("\\|.*", "", transcript),
      comparison = paste(t1, "vs", t2, sep = "_")
    ) |>
    filter(transcript != "ENSMUST00000118878.9") %>%
    dplyr::filter(qvalue < q_cut)
  
  return(dtu_res)
}

## Define comparison
pairs <- list(
  c("heart", "striatum"),
  c("heart", "cortex"),
  c("heart", "hippocampus"),
  c("heart", "cerebellum")
  
)

## Run DTE Analysis Across All Pairs
library(purrr)
dtu_all <- map_dfr(pairs, ~ run_dtu_pair(y0, .x[1], .x[2]))

# Create transcript labels
create_transcript_labels <- function(transcript_ids) {
  # Define labels
  transcript_labels <- c(
    "NCBI_NM_001290729.3_nTaf1_7168bp" = "NM_001290729.3(NCBI_nTaf1)",
    "ENSMUST00000143908.8" = "ENSMUST00000143908.8(Taf1-204)",
    "NM_001405959.2_Mus" = "NM_001405959.2(NCBI_cTaf1)"
  )
  
  # Apply the mapping to each transcripts ID
  sapply(transcript_ids, function(id) {
    if (id %in% names(transcript_labels)) {
      return(transcript_labels[id])
    } else {
      return(id)  # if no mapping exists, keep original ID
    }
  }, USE.NAMES = FALSE)
}

## Run DTU Analysis Across All Pairs
dtu_all$transcript_label <- create_transcript_labels(dtu_all$transcript)

# Create comparison labels
dtu_all <- dtu_all %>%
  mutate(
    comparison_label = case_when(
      comparison == "heart_vs_striatum" ~ "heart vs striatum",
      comparison == "heart_vs_cortex" ~ "heart vs cortex", 
      comparison == "heart_vs_hippocampus" ~ "heart vs hippocampus",
      comparison == "heart_vs_cerebellum" ~ "heart vs cerebellum",
      TRUE ~ comparison
    )
  )


# Add significance level annotations
dtu_all <- dtu_all %>%
  mutate(
    significance = case_when(
      qvalue < 0.01 ~ "**",
      qvalue < 0.05 ~ "*",
      qvalue < 0.10 ~ "•",
      TRUE ~ ""
    ),
    direction = ifelse(log2FC > 0, "Up-regulated", "Down-regulated"),
    # Add absolute values for sorting
    abs_log2FC = abs(log2FC)
  )

## creat lollipop plot
library(scales)

p1 <- ggplot(dtu_all, aes(x = log2FC, y = reorder(transcript_label, log2FC))) +
  geom_segment(aes(x = 0, xend = log2FC, y = transcript_label, yend = transcript_label)) +
  #color = "gray60", size = 0.5) +
  geom_point(aes(fill = qvalue), size = 7,
             shape = 21, color = "white", stroke = 0.5) +
  geom_text(aes(label = significance), 
            hjust = 0.5, vjust = -0.05, size = 6.5, fontface = "bold") +
  #facet_wrap(~comparison_label, scales = "free_y", ncol = 1) +
  facet_wrap(~comparison_label, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "darkred", high = "lightcoral", 
                      name = "q-value", trans = "log10",
                      labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray95", size = 0.3),
    #panel.spacing.y = unit(1, "cm"),
    strip.text = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 16,face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    plot.title = element_text(size = 18, hjust = 0.5),
    #plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  labs(
    x = "log2 Fold Change",
    y = "Transcript",
    title = "Significant Differential Transcript Usage (DTU) Analysis",
    #subtitle = "Lollipop plot showing isoform proportion changes",
    caption = "*** p<0.001, ** p<0.01, * p<0.05, • p<0.10"
  )

print(p1)