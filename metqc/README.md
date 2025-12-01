
# metqc

The goal of metqc is to provide Quartet-based QC tools for Metabolomics.

## Installation
You can install this package from github

```R
install.packages("devtools")
devtools::install_github("chinese-quartet/quartet-metabolism-qc", subdir="metqc")
```

## Example

```R
library(metqc)
### get example data and reprot template path
example_sample <- system.file("extdata", "sample_data.csv", package = "metqc")
example_metadata <- system.file("extdata", "sample_metadata.csv", package = "metqc")

## Get performance report for metabolomics data
get_performance(dt_file = sample_data, metadata_file = sample_metadata)

## Count SNR and plot related PCA plot
count_snr(dt_file = sample_data, metadata_file = sample_metadata)
#> [1] 4.511966

## Count RC and plot related scatter plot
count_rc(dt_file = sample_data, metadata_file = sample_metadata)
#> [1] 0.7637139

## Count Recall 
count_recall(dt_file = sample_data, metadata_file = sample_metadata)
#> [1] 0.1649485

#### generate report
### get reprot template path
report_template <- system.file("extdata", "quartet_template.docx", package = "metqc")

### get metabolomics metrics data and generate report
met_result = get_performance(dt_file = sample_data, metadata_file = sample_metadata)
generate_met_report(qc_result = met_result, report_template = report_template)
```
