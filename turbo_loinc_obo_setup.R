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

### get mappings with BioPortal
# or string-search somewhere?
# start with public endpoint but eventually switch to appliance

#### these are functioning like globals so they don't have to be passes to the function
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
# what are the chances that a mapping query will return 0 mappings, or that it will return multiple pages?

bp.map.retreive.and.parse <- function(term.list) {
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
      # CUI, LOOM, "same URI", etc. Probably only LOOM will be useful
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

#  http://data.bioontology.org/documentation#nav_search
#  http://data.bioontology.org/search?q=melanoma

# 9cf735c3-a44a-404f-8b2f-c49d48b2b8b2

# An API Key is required to access any API call. It can be provided in three ways:
#
#   Using the apikey query string parameter
# Providing an Authorization header: Authorization: apikey token=your_apikey (replace `your_apikey` with your actual key)

# Parameters

# ontologies={ontology_id1,ontology_id2,ontology_id3}

# require_exact_match={true|false} // default = false
# suggest={true|false} // default = false. Will perform a search specifically geared towards type-ahead suggestions.
# also_search_views={true|false} // Include ontology views in the search. default = false
# require_definitions={true|false} // default = false
# also_search_properties={true|false} // default = false
# also_search_obsolete={true|false} // default = false (exclude obsolete terms)
# cui={C0018787,C0225807,C0018787} // Allows restricting query by CUIs. default = null (no restriction)
# semantic_types={T023,T185,T061} // Allows restricting query by Semantic Types (TUI). default = null (no restriction)
# include={prefLabel, synonym, definition, notation, cui, semanticType} // default = (see Common Parameters section)
# page={integer representing the page number} // default = 1
# pagesize={integer representing the size of the returned page} // default = 50

bioportal.string.search <- function(current.string) {
  # current.string <- 'asthma'
  print(current.string)
  prepared.get <-
    paste0(
      'http://data.bioontology.org/search?q=',
      current.string  ,
      '&include=prefLabel,synonym',
      '&pagesize=999'
    )
  prepared.get <- URLencode(prepared.get, reserved = FALSE)
  search.res.list <-
    httr::GET(url = prepared.get,
              add_headers(Authorization = paste0(
                "apikey token=", config$bioportal.api.key
              )))
  
  search.res.list <- rawToChar(search.res.list$content)
  search.res.list <- jsonlite::fromJSON(search.res.list)
  search.res.list <- search.res.list$collection
  
  # print(search.res.list$links$ontology)
  
  if (is.data.frame(search.res.list)) {
    if (nrow(search.res.list) > 0) {
      ontology <- search.res.list$links$ontology
      #  , 'ontologyType'
      search.res.list <- search.res.list[, c('prefLabel', '@id')]
      colnames(search.res.list) <- c('prefLabel', 'iri')
      search.res.list <-
        cbind.data.frame(search.res.list, 'ontology' = ontology)
      search.res.list$rank <- 1:nrow(search.res.list)
      return(search.res.list)
    }
  }
}

# temp <- bioportal.string.search(current.string = 'asthma')
# temp <- bioportal.string.search(current.string = 'leukocyte')
# temp <- bioportal.string.search(current.string = 'leukocytes')

# > str(potential.cells)
# 'data.frame':	22 obs. of  2 variables:
#   $ PartTypeVal    : chr  "LP15429-1" "LP14328-6" "LP15100-8" "LP14539-8" ...
# $ PartDisplayName: chr  "Base excess" "Basophils" "Blasts" "Eosinophils" ...

# pendulum swings completely the other way here
# everything is an arguemt now

# ontology.filter = &ontology=uberon,cl,pr,chebi
# kept.row.count= 9
#  ir rank important?

# do singlular and plurals together in the same funtion?
# are alternative labels getting searched? yes, but ranked lower
#   acetaminophen paracetamol

# see https://www.ebi.ac.uk/ols/docs/api
ols.serch.term.labels.universal <-
  function(current.string,
           current.id,
           strip.final.s = FALSE,
           ontology.filter,
           kept.row.count = 9,
           req.exact = 'false') {
    if (strip.final.s) {
      current.string <-
        sub(pattern = "s$",
            replacement = "",
            x = current.string)
    }
    
    singular.lc <- current.string
    
    print(singular.lc)
    
    # or just try url encoding?
    
    current.string <-
      gsub(pattern = "[[:punct:] ]",
           replacement = ",",
           x = current.string)
    
    # current.string <-
    #   gsub(pattern = "[[:punct:]]",
    #        replacement = " ",
    #        x = current.string)
    
    #
    # current.string <-
    #   gsub(pattern = "^ +",
    #        replacement = "",
    #        x = current.string)
    #
    # current.string <-
    #   gsub(pattern = " +$",
    #        replacement = "",
    #        x = current.string)
    #
    # current.string <-
    #   gsub(pattern = " +",
    #        replacement = "+",
    #        x = current.string)
    
    # print(current.string)
    
    prepared.query <- paste0(
      "https://www.ebi.ac.uk/ols/api/search?q={",
      current.string,
      "}&type=class&local=true",
      ontology.filter ,
      "&rows=",
      kept.row.count,
      '&exact=',
      req.exact,
      "&fieldList=iri,label,synonym,short_form,obo_id,ontology_name,ontology_prefix"
      # " & ",
      # 'queryFields=label,synonym'
    )
    
    # print(prepared.query)
    
    
    # fieldList
    # Specifcy the fields to return, the defaults are {iri,label,short_form,obo_id,ontology_name,ontology_prefix,description,type}
    #
    # queryFields
    # Specifcy the fields to query, the defaults are {label, synonym, description, short_form, obo_id, annotations, logical_description, iri}
    
    # print(prepared.query)
    
    ols.attempt <-
      httr::GET(prepared.query)
    
    ols.attempt <- ols.attempt$content
    ols.attempt <- rawToChar(ols.attempt)
    ols.attempt <- jsonlite::fromJSON(ols.attempt)
    ols.attempt <- ols.attempt$response$docs
    if (is.data.frame(ols.attempt)) {
      if (nrow(ols.attempt) > 0) {
        ols.attempt$query <- singular.lc
        ols.attempt$loinc.part <- current.id
        
        ols.attempt$rank <- 1:nrow(ols.attempt)
        ols.attempt$label <- tolower(ols.attempt$label)
        ols.attempt$query <- tolower(ols.attempt$query)
        
        return(ols.attempt)
      }
    }
  }

update.accounting <- function(recent.successes,
                              source.id,
                              target.term,
                              authority) {
  print("before update")
  print(length(current.component.mapping.complete))
  print(length(current.needs.component.mapping))
  print(nrow(current.component.mapping.frame))
  
  print("update row count")
  print(nrow(recent.successes))
  
  temp <- unlist(recent.successes[, source.id])
  source.term <-
    paste0("http://purl.bioontology.org/ontology/LNC/", temp)
  
  authority <- rep(x = authority , nrow(recent.successes))
  
  current.component.mapping.complete <<-
    union(current.component.mapping.complete, temp)
  current.needs.component.mapping <<-
    setdiff(current.needs.component.mapping, temp)
  
  temp <-
    cbind(temp,
          source.term,
          recent.successes[, target.term],
          authority)
  
  colnames(temp) <-
    c('source.id', 'source.term', 'target.term', 'authority')
  
  current.component.mapping.frame <<-
    rbind(current.component.mapping.frame, temp)
  
  print("after update")
  print(length(current.component.mapping.complete))
  print(length(current.needs.component.mapping))
  print(nrow(current.component.mapping.frame))
  
  print(sort(table(
    current.component.mapping.frame$authority
  )))
  
}
