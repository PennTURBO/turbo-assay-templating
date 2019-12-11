names(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report)


# filter.form <-
#   formula(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$time == "24H")

partlist <- c("time", "scale", "prop", "meth")

part.tab <- function(partname) {
  # part <- "prop"
  part <- partname
  
  temp <-
    sort(
      table(
        PDS_labs_promising_LOINC_annotated_only_and_tabulation_report[, part],
        useNA = 'always'
      ),
      decreasing = TRUE
    )
  
  temp <- cbind.data.frame(names(temp), as.numeric(temp))
  
  names(temp) <- c(part, "count")
  
  print(temp)
}

placeholder <- lapply(partlist, part.tab)
