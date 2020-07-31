loinc_to_obo_mapping_reviewed <-
  read_csv("data/loinc_to_obo_mapping_reviewed.csv")

loinc_to_obo_mapping_reviewed <-
  left_join(x = loinc_to_obo_mapping_reviewed,
            y = Part,
            by  = c("PartNumber" = "PartNumber"))


loinc_to_obo_mapping_reviewed <-
  left_join(x = loinc_to_obo_mapping_reviewed,
            y = turbo.labels,
            by  = c("iris" = "s"))

write.csv(loinc_to_obo_mapping_reviewed, file = "data/loinc_to_obo_mapping_re_review.csv", row.names = FALSE)
