library(readr)

turbo_assay_template_headerless.fn <-
  "build/turbo_assay_template_headerless.csv"

report_analyte_eligible.fn <- "build/report_analyte_eligible.csv"

turbo_assay_template_headerless <-
  read_csv(
    turbo_assay_template_headerless.fn,
    col_names = FALSE,
    col_types = cols(.default = "c")
  )

report_analyte_eligible <-
  read_csv(report_analyte_eligible.fn)

turbo_assay_template_headerless$X18[turbo_assay_template_headerless$X1 %in% report_analyte_eligible$assay] <-
  turbo_assay_template_headerless$X24[turbo_assay_template_headerless$X1 %in% report_analyte_eligible$assay]

turbo_assay_template_headerless$X24[turbo_assay_template_headerless$X1 %in% report_analyte_eligible$assay] <-
  ''

turbo_assay_template_headerless[is.na(turbo_assay_template_headerless)] <-
  ''

write.table(
  x = turbo_assay_template_headerless,
  file = turbo_assay_template_headerless.fn,
  sep = ",",
  col.names = FALSE,
  row.names = FALSE
)

write.table(
  x = turbo_assay_template_headerless,
  file = "build/turbo_assay_template_headerless.tsv",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# write.csv(turbo_assay_template_headerless,
#           file = turbo_assay_template_headerless.fn,
#           row.names = FALSE)