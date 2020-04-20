# dput(sort(c("COMPONENT", "analyte.core.vals", "LOINC", "SYSTEM", "LOINC_COUNT",
#   "METHOD", "PROPERTY", "SCALE", "TIME", "PartName.challenge",
#   "PartName.divisor", "PartName.suffix", "PartName.method", "method.count",
#   "system.iri", "time.iri", "search.counts", "analyte.vals", "COMPONENT.vals",
#   "search.vals", "target.term", "iri.from.search", "final.comp"
# )))


further.exploration <-
  c(
    "PartName.challenge",
    "PartName.divisor",
    "PartName.method",
    "PartName.suffix",
    "time.iri",
    "PROPERTY"
  )

further.exploration.frame <-
  pds.with.loinc.parts[, further.exploration]



# further.exploration <- c(further.exploration)

beaten.to.death <-
  setdiff(colnames(pds.with.loinc.parts), further.exploration)

beaten.to.death.frame <- pds.with.loinc.parts[, beaten.to.death]

lapply(further.exploration, function(current.further) {
  print(current.further)
  print(table(pds.with.loinc.parts[, current.further]))
  return()
})

prop.tab <- table(pds.with.loinc.parts$PROPERTY)
prop.tab <- cbind.data.frame(names(prop.tab), as.numeric(prop.tab))
colnames(prop.tab) <- c('PROPERTY', 'count')
prop.tab$PROPERTY <- as.character(prop.tab$PROPERTY)

# or could ahve gotten this from one of those pds reveiwed files
prop.dn <-
  Part[Part$PartTypeName == "PROPERTY", c("PartNumber", "PartDisplayName")]
# colnames(prop.text) <- c('PartNumber', 'PROPERTY.text')

prop.dn.tab <- base::merge(
  x = prop.tab,
  y = prop.dn,
  by.x = 'PROPERTY',
  by.y = 'PartNumber',
  all.x = TRUE
)
