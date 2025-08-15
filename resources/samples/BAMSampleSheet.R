args <- commandArgs(TRUE)

library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(jsonlite)
library(here)

# Read Excel file
file_path <- here("resources", "samples", args[1])
if (!file.exists(file_path)) stop("File does not exist: ", file_path)

metadata <- read_excel(file_path, .name_repair = "minimal")

cat("Initial metadata dimensions:", dim(metadata), "\n")
print(head(metadata))

# Assign tumor/normal using Tissue Type
metadata <- metadata %>%
  mutate(
    sample_name = str_extract(`Sample ID`, "^MMRF_\\d+"),
    tissue_group = case_when(
      `Tissue Type` == "Tumor" ~ "tumor",
      `Tissue Type` == "Normal" ~ "normal",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tissue_group))

cat("After tissue_group assignment:", dim(metadata), "\n")
table(metadata$tissue_group)

# Build sample pairs (row-wise pairing)
sample_pairs <- metadata %>%
  select(
    patient_id = `Case ID`,
    sample_name,
    tissue_group,
    bam = `File Name`
  ) %>%
  group_by(patient_id) %>%
  group_map(~ {
    tumors <- filter(.x, tissue_group == "tumor")
    normals <- filter(.x, tissue_group == "normal")
    
    if (nrow(normals) == 0 || nrow(tumors) == 0) return(NULL)
    
    n_pairs <- min(nrow(tumors), nrow(normals))
    
    data.frame(
      patient_id = rep(.y$patient_id, n_pairs),
      sample_name = tumors$sample_name[1:n_pairs],
      tumor_bam = tumors$bam[1:n_pairs],
      normal_bam = normals$bam[1:n_pairs],
      stringsAsFactors = FALSE
    )
  }) %>%
  compact() %>%
  bind_rows()


# Skip JSON writing if empty, otherwise use here()
if (nrow(sample_pairs) == 0) {
  message("No sample pairs generated; skipping JSON output")
} else {
  pyfile <- here("resources", "samples", "sample_pairs.py")
  cat("sample_pairs = [\n", file = pyfile)
  for (i in seq_len(nrow(sample_pairs))) {
    cat(
      "    {\n",
      sprintf('        "patient_id": "%s",\n', sample_pairs$patient_id[i]),
      sprintf('        "sample_name": "%s",\n', sample_pairs$sample_name[i]),
      sprintf('        "tumor_bam": "%s",\n', sample_pairs$tumor_bam[i]),
      sprintf('        "normal_bam": "%s"\n', sample_pairs$normal_bam[i]),
      if (i < nrow(sample_pairs)) "    },\n" else "    }\n",
      sep = "", file = pyfile, append = TRUE
    )
  }
  cat("]\n", file = pyfile, append = TRUE)
  message("sample_pairs.py written to: ", pyfile)
}

# Print preview
print(sample_pairs)