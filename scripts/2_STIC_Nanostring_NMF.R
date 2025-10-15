library(NMF)
root           <- "Location"
data_dir       <- normalizePath(file.path(root, "data", "raw"), mustWork = TRUE)

root           <- "/Users/sd_wo/Documents/projects/nanostring-paper-code"
data_dir       <- normalizePath(file.path(root, "data", "raw"), mustWork = TRUE)
data_processed <- normalizePath(file.path(root, "data", "processed"), mustWork = FALSE)
rdata_dir      <- file.path(data_processed, "Rdata")
fig_dir        <- file.path(root, "outputs", "figures", "02_analysis")
table_dir      <- file.path(root, "outputs", "tables")

dir.create(data_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,       recursive = TRUE, showWarnings = FALSE)

message("root:      ", root)
message("data_dir:  ", data_dir)
message("rdata_dir: ", rdata_dir)
message("fig_dir:   ", fig_dir)

geo_fn <- file.path(root, "script", "STIC_Nanostring_Custom_Functions.R")
if (file.exists(geo_fn)) source(geo_fn)
```

### Sample data adjustment
```{r}
# ---------------------------------------------------------
# Load data and annotation data
# ---------------------------------------------------------

# Load processed count data (vst_limma.xlsx) from processed/Rdata
count_data <- read_excel(file.path(data_dir, "count_data.xlsx")) %>%
  as.data.frame() %>%
  { rownames(.) <- .[[1]]; . } %>%
  dplyr::select(-1) %>%
  as.matrix()

# Load annotation metadata (Nanostring_metadata_submit.xlsx) from raw/Data
anno_data <- read_xlsx(file.path(data_dir, "Nanostring_metadata_submit.xlsx")) %>%
  as.data.frame()

colnames(anno_data) <- gsub(" ", "_", colnames(anno_data))

anno_data <- anno_data %>%
  mutate(Diagnosis = case_when(
    Diagnosis == "STIC" ~ "STIC",
    Diagnosis == "High grade serous carcinoma" ~ "HGSC",
    Diagnosis == "Normal fallopian tube epithelium" ~ "NFT",
    Diagnosis == "STIL" ~ "STIC",
    Diagnosis == "p53 signature" ~ "p53",
    Diagnosis == "Clear cell carcinoma" ~ "CCC",
    Diagnosis == "Endometrial endometrioid carcinoma" ~ "EC",
    Diagnosis == "Endometrioid carcinoma" ~ "EC",
    Diagnosis == "Normal mesothelium" ~ "Meso",
    TRUE ~ Diagnosis
  ),
  BRCA_category = recode(
    BRCA_category,
    "negative" = "Negative",
    .default   = BRCA_category
  ))

anno_data <- as.data.frame(anno_data)
rownames(anno_data) <- anno_data$SegmentDisplayName


# epithelial subset
anno_data.epi <- anno_data %>%
  dplyr::filter(
    type == "epithelium",
    !Diagnosis %in% c("CCC", "EC")  # exclude other carcinoma types
  )
count_data.epi <- count_data[, row.names(anno_data.epi)]

## Adjust data type
anno_data.epi <- anno_data.epi %>%
  mutate(
    type = as.factor(type),
    Diagnosis = as.factor(Diagnosis),
    Morphology_category = as.factor(Morphology_category),
    Ki67_percentage = as.integer(Ki67_percentage),
    p53_pattern = as.factor(p53_pattern),
    Stromal_lymphocyte = as.factor(Stromal_lymphocyte),
    BRCA_category = as.factor(BRCA_category),
    Molecular_Subtype = factor(
      Molecular_Subtype,
      levels = c("NFT", "Dormant", "Mixed", "Immunoreactive", "Proliferative", "HGSC")
    )
  )


anno_data.epi.subset<- anno_data.epi[!anno_data.epi$Diagnosis %in% c("CCC", "EC", "NFT", "HGSC"), ]
count_data.epi.subset <- count_data.epi[, row.names(anno_data.epi.subset)]

geoData <- count_data.epi.subset
geoData[geoData<0] <- 0

set.seed(12345)
res <- nmf(
  geoData,
  rank     = 2:6,
  method   = "brunet",
  nrun     = 100,
  seed     = 12345,
  .options = list(v = TRUE, P = TRUE, k = TRUE), 
)
