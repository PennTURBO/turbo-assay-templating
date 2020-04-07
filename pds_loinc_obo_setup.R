library(config)
library(RJDBC)
library(readr)
library(reshape2)
# see also https://jangorecki.gitlab.io/data.cube/library/data.table/html/dcast.data.table.html
library(dplyr)
library(httr)
# different (bioconductor) install method!
# probably won't use anyway
#   library(rols)
library(tm)

config <- config::get(file = "loinc_in_obo.yaml")

### get mappings or search with bioportal
# start with public endpoint but eventually switch to appliance
api.base.uri <- "http://data.bioontology.org/ontologies"
api.ontology.name <- "LOINC"
term.ontology.name <- "LNC"
term.base.uri <-
  paste0("http://purl.bioontology.org/ontology",
         "/",
         term.ontology.name)
api.family <- "classes"
# source.term <- "http://purl.bioontology.org/ontology/LNC/LP17698-9"
api.method <- "mappings"
# what are the chances that a mapping query will return 0 ammpings, or that it will return multiple pages?

retreive.and.parse <- function(term.list) {
  outer <- lapply(term.list, function(current.term) {
    # current.term <- "LP102314-4"
    print(current.term)
    current.uri <- paste0(term.base.uri, "/", current.term)
    encoded.term <- URLencode(current.uri, reserved = TRUE)
    prepared.get <-
      paste(api.base.uri,
            api.ontology.name,
            api.family,
            encoded.term,
            api.method,
            sep = "/")
    mapping.res.list <-
      httr::GET(url = prepared.get,
                add_headers(
                  Authorization = paste0("apikey token=", config$bioportal.api.key)
                ))
    
    mapping.res.list <- rawToChar(mapping.res.list$content)
    
    mapping.res.list <- jsonlite::fromJSON(mapping.res.list)
    if (length(mapping.res.list) > 0) {
      # CUI, LOOM, "same URI", etc. Porbably only LOOM will be useful
      mapping.methods <- mapping.res.list$source
      
      source.target.details <-
        lapply(mapping.res.list$classes, function(current.mapping) {
          source.target.terms <- current.mapping$`@id`
          source.target.ontologies <- current.mapping$links$ontology
          return(c(
            rbind(source.target.terms, source.target.ontologies)
          ))
        })
      
      source.target.details <-
        do.call(rbind.data.frame, source.target.details)
      colnames(source.target.details) <-
        c("source.term",
          "source.ontology",
          "target.term",
          "target.ontology")
      
      source.target.details <-
        cbind.data.frame(source.target.details, mapping.methods)
      return(source.target.details)
    }
    
  })
}