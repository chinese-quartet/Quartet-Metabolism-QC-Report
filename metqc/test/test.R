devtools::load_all()

devtools::document()


library(metqc)

### get example data and reprot template path
example_sample <- system.file("extdata", "sample_data.csv", package = "metqc")
example_metadata <- system.file("extdata", "sample_metadata.csv", package = "metqc")
example_sample
example_metadata
## Get performance report for metabolomics data
get_performance(dt_file = example_sample, metadata_file = example_metadata)

get_performance()

report_template <- system.file("extdata", "quartet_template.docx", package = "metqc")

### get metabolomics metrics data and generate report
met_result = get_performance(dt_file = example_sample, metadata_file = example_metadata)
generate_met_report(qc_result = met_result, report_template = report_template,report_dir = './test/')
