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
rownames(Taf1) <- Taf1$sample #set row name
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
                             levels = c("striatum", "hippocampus", "cortex", "cerebellum", "muscle"))
colData(y0)$sex    <- factor(colData(y0)$sex)

run_dte_pair <- function(se_obj, t1, t2, cov = "sex",
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
  
  ## Organize output; only DTE results are retured here
  dte_res <- mcols(y) |>
    tibble::as_tibble() |>
    mutate(
      transcript = rownames(assay(y)),
      transcript = sub("\\|.*", "", transcript),
      comparison = paste(t1, "vs", t2, sep = "_")
    ) |>
    filter(transcript != "ENSMUST00000118878.9") %>%
    dplyr::filter(qvalue < q_cut)
  
  return(dte_res)
}

## Define comparison
pairs <- list(
  c("muscle", "striatum"),
  c("muscle", "cortex"),
  c("muscle", "hippocampus"),
  c("muscle", "cerebellum")
)

## Run DTE Analysis Across All Pairs
library(purrr)
dte_all <- map_dfr(pairs, ~ run_dte_pair(y0, .x[1], .x[2]))

# Create comparison labels
dte_all <- dte_all %>%
  mutate(
    comparison_label = case_when(
      comparison == "muscle_vs_striatum" ~ "muscle vs striatum",
      comparison == "muscle_vs_cortex" ~ "muscle vs cortex", 
      comparison == "muscle_vs_hippocampus" ~ "muscle vs hippocampus",
      comparison == "muscle_vs_cerebellum" ~ "muscle vs cerebellum",
      TRUE ~ comparison
    )
  )

# Add significance level annotations
dte_all <- dte_all %>%
  mutate(
    significance = case_when(
      qvalue < 0.01 ~ "**"
      qvalue < 0.05 ~ "*",
      qvalue < 0.10 ~ "•",
      TRUE ~ ""
    ),
    # Add categorical variables for color mapping
    direction = ifelse(log2FC > 0, "Up-regulated", "Down-regulated"),
    # Add absolute values for sorting
    abs_log2FC = abs(log2FC)
  )

## creat lollipop plot
library(scales)

p2 <- ggplot(dte_all, aes(x = log2FC, y = reorder(transcript_label, log2FC))) +
  geom_segment(aes(x = 0, xend = log2FC, y = transcript_label, yend = transcript_label)) +
  #color = "gray60", size = 0.5) +
  geom_point(aes(fill = qvalue), size = 7,
             shape = 21, color = "white", stroke = 0.5) +
  geom_text(aes(label = significance), 
            hjust = 0.5, vjust = -0.05, size = 6.5, fontface = "bold") +
  facet_wrap(~comparison_label, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "#27408B", high = "#87CEFF", 
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
    axis.text.y = element_text(size = 11, face = "bold"),
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
    title = "Significant Differential Transcript Expression (DTE) Analysis",
    #subtitle = "Lollipop plot showing isoform proportion changes",
    caption = "*** p<0.001, ** p<0.01, * p<0.05, • p<0.10"
  )

print(p2)
