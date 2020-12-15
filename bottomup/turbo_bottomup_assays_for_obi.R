# todo
# could break the results out for analyte, cytometry and "other" assays

# include counts more pervasively in the result instances, result items, order and specimens dataframes ?

# requires that the robot jar and shell script are on the path
# robot morror commands are slow... use a makefile/check for new file strategy?
# could also make terms lists and do extracts for a smaller merged file and faster downstream processing?
# pr is by far the largest

# will likely require a VPN and ssh tunnel
# run from bottom folder
# setup script below requires a few "pre-commit" text files with version info, etc.

# requires that one or more graphdb servers are pre-populated with XXX

# requires that several OBO ontologies have been "robot mirrored" into a "purl.obolibrary.org/obo/" folder
# chebi, cl, go, obi, pato, pr
# and then "robot merged" into "purl.obolibrary.org/obo/obo_merged_for_assays.owl"

# this script writes tables to be annotated to the build folder
# and then assumes that they have been annoated
# so reads them from the build/annotated folder

# dbxrefs, like LOINC, are lost in uniqification at this point



source(
  "https://raw.githubusercontent.com/PennTURBO/turbo-globals/master/turbo_R_setup.R"
)

library(tidytext)
library(googlesheets4)

####

# where can we find OBO terms required for assay modeling?
primary.gdb.endpoint <- "http://localhost:7200"
primary.gdb.repo <- "bottom_up"

# what should we strip from the end of an ontology IRI or graph name
# as one step in finding the ontology abbreviation
owl.extension <- ".owl"
# what prefix should be removed?
obo.iri.prefix <- "http://purl.obolibrary.org/obo/"

# where can we find a TURBO graph with patients identified by EMPIS
# the assay modelling will be based on lab results avaialble for those patients
patient.scoping.repo <- "mam_drivetrain_prd"

# what's a small query we can run to check our conenction to the EHR?
db.con.check.q <- "select * from MDM.R_LOCATION"

# sometimes the ids.labels datfarame ends up empty
#  when this script is run beginning to end
# I added 10 second sleep breaks between the steps last time I ran it
#  and it turned out OK
# 5 seconds was ok at least once too
sleepsecs <- 2

# SQL setup for finding lab results and reference material
# there are methods for retrieving similar data from cached CSV fies below
pdsDriver <-
  JDBC(driverClass = "oracle.jdbc.OracleDriver",
       classPath = config$oracle.jdbc.path)

pds.con.string <- paste0(
  "jdbc:oracle:thin:@//",
  config$pds.host,
  ":",
  config$pds.port,
  "/",
  config$pds.database
)

# the EHR may have some housekeeping columns in all tables
# don't retrain these
unwanted.cols <-
  c("ACTIVE_YN",
    "MDM_LAST_UPDATE_DATE",
    "MDM_INSERT_UPDATE_FLAG")

lab.result.instances.res.fp <-
  "local/_SELECT_RLRI_PK_LAB_RESULT_ITEM_ID_ROI_PK_ORDER_ITEM_ID_LR_PK_LA_202012101858.csv"

roi.res.fp <- "local/R_ORDER_ITEM_202012101244.csv"

rlri.res.fp <- "local/R_LAB_RESULT_ITEM_202012111046.csv"

# two header rows that should be included for ROBOT templating
obi_assay_headers.fp <-  "data/obi_assay_headers.txt"

# how to create IRIs for the new assays?
id.prefix <- "http://purl.obolibrary.org/obo/OBI_"
id.start <- 3100001

# assuming there is a single OBI github issue that describes the content we'll be generating
term_tracker_id <- "https://github.com/obi-ontology/obi/issues/1253"

term_editor <-
  "Mark A. Miller, ORCID:0000-0001-9076-6066 & Chris Stoeckert, ORCID ORCID:0000-0002-5714-991X"
def_source <- "PennTURBO team"

# how should the results be written
write.robot.tsv <- TRUE

# with gsheets.auth.email,
# it's likely that the gsheet write will complete without user interaction
interactively.write.gsheet <- TRUE
gsheets.auth.email <- 'mamillerpa@gmail.com'

# these are probably unnecessary now
# instead of doing this with rules in the script,
# overrides can be asserted for all terms in the
# result instances, result items and orders lookup tables
analyte.normalization.rules <- c('', 'bicarbonate', 'oxygen', '')
names(analyte.normalization.rules) <-
  c(' molecular entity',
    'hydrogencarbonate',
    'dioxygen',
    ' specimen')

# dioxygen vs oxygen molent?
# move (human) suffix for proteins forward and drop parens ?

####

split.curie.col <-
  function(curie.col,
           prefix.name = "term.prefix",
           term.rhs.name = "term.rhs") {
    temp <- strsplit(curie.col, split = ":", fixed = TRUE)
    temp <- do.call(rbind.data.frame, temp)
    colnames(temp) <- c(prefix.name, term.rhs.name)
    return(temp)
  }

onto.iri.to.prefix <- function(onto.iri.col) {
  temp <- onto.iri.col
  temp <-
    sub(
      pattern = obo.iri.prefix,
      replacement = "",
      x = temp,
      fixed = TRUE
    )
  temp <-
    sub(
      pattern = owl.extension,
      replacement = "",
      x = temp,
      fixed = TRUE
    )
  return(temp)
}

# the orders may contain specimen and target entity knwoledge
# (the result items sometimes contain specimen and order knowledge, too,
# but that isn't being utlized yet)
compare.to.order <- function(data,
                             order.col.name,
                             comparitor.col.name) {
  temp  <- rep("OK",
               nrow(data))
  order.col.na.flag <- is.na(data[, order.col.name])
  comparitor.col.na.flag <- is.na(data[, comparitor.col.name])
  no.na <- order.col.na.flag + comparitor.col.na.flag == 0
  temp[no.na &
         data[, order.col.name] == data[, comparitor.col.name]] <-
    "duplication"
  temp[no.na &
         data[, order.col.name] != data[, comparitor.col.name]] <-
    "conflict"
  return(temp)
}

summarize.comparison <-
  function(data,
           comparison.col.name,
           reportable.levels = c("duplication", "conflict"),
           reported.columns) {
    comparison.temp <- pull(data, comparison.col.name)
    flag.it <- comparison.temp %in% reportable.levels
    data <- unique(data[flag.it , reported.columns])
  }

obo.curie.to.iri <- function(curies) {
  inner <- sub(pattern = ":",
               replacement = "_",
               x = curies)
  inner <- paste0(obo.iri.prefix, inner)
  return(inner)
}

obo.iri.to.curie <- function(iris) {
  inner <- sub(pattern = obo.iri.prefix,
               replacement = "",
               x = iris)
  inner <- sub(pattern = "_",
               replacement = ":",
               x = inner)
  return(inner)
}

# SPARQL for detemining whose lab resutls should be used for bootstrapping assays
get.empis.from.graphdb <-
  function(inner.endpoint = config$my.graphdb.base,
           inner.repo = config$my.selected.repo,
           inner.sa = saved.authentication) {
    empi.q <- "
PREFIX pmbb: <http://www.itmat.upenn.edu/biobank/>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX ontologies: <http://transformunify.org/ontologies/>
select
#distinct
#?d ?l
#*
?p ?r
where {
    graph pmbb:expanded {
        ?s a <http://purl.obolibrary.org/obo/IAO_0000028> ;
           obo:BFO_0000050 ?c ;
           <http://transformunify.org/ontologies/TURBO_0010094> ?r .
        ?c a <http://purl.obolibrary.org/obo/IAO_0000578> ;
           obo:BFO_0000051 ontologies:TURBO_0010295 ;
           obo:IAO_0000219 ?p .
        ?p a obo:NCBITaxon_9606 .
    }
}
"

# print(inner.endpoint)
# print(inner.repo)

empi.res <-
  q2j2df(
    query = empi.q,
    auth = inner.sa,
    endpoint = inner.endpoint,
    repo = inner.repo
  )

return(empi.res)

  }

# SQL *for the UPHS EHR* for getting lab results
# we're currently using a chced file instead
get.lab.res.for.empis.from.pds <- function(empis) {
  # this still doesn't seem to handle all of the possible cases
  # such as pdsConnection defined and still active
  # defined but expired
  # not defined yet
  if (exists("pdsConnection")) {
    print("exists")
    result <- tryCatch({
      temp.res <- dbGetQuery(pdsConnection, db.con.check.q)
      print(temp.res)
    }, error = function(e) {
      print("exists but inactive. Trying to reconnect.")
      pdsConnection <-
        dbConnect(pdsDriver,
                  pds.con.string,
                  config$pds.user,
                  config$pds.pw)
      temp.res <- dbGetQuery(pdsConnection, db.con.check.q)
      print(temp.res)
    })
  } else {
    print("connection does not exist yet.")
    pdsConnection <-
      dbConnect(pdsDriver,
                pds.con.string,
                config$pds.user,
                config$pds.pw)
    temp.res <- dbGetQuery(pdsConnection, db.con.check.q)
    print(temp.res)
  }
  
  lab.res.q <- paste0(
    "
                      SELECT
  LR.PK_LAB_RESULT_ID,
  ROI.PK_ORDER_ITEM_ID,
  RLRI.PK_LAB_RESULT_ITEM_ID,
  LR.SPECIMEN,
  LR.SPECIMEN_COLLECTION_METHOD
  FROM
  MDM.PATIENT_ENCOUNTER PE
  INNER JOIN MDM.LAB_RESULTS LR ON
  PE.PK_PATIENT_ENCOUNTER_ID = LR.FK_PATIENT_ENCOUNTER_ID
  INNER JOIN MDM.R_LAB_RESULT_ITEM RLRI ON
  LR.FK_LAB_RESULT_ITEM_ID = RLRI.PK_LAB_RESULT_ITEM_ID
  INNER JOIN MDM.R_ORDER_ITEM ROI ON
  LR.FK_ORDER_ITEM_ID = ROI.PK_ORDER_ITEM_ID
  WHERE
  PE.EMPI IN ( ",
    paste0("'", paste0(graph.empis$r, collapse = "', '"), "'") ,
    " )
  AND LR.RESULT_VALUE_NUM IS NOT NULL
  AND UNIT_OF_MEASURE IS NOT NULL
  AND RESULT_TYPE NOT IN ('ALPHA', 'FREETEXT')"
  )
  
  print(Sys.time())
  timed.system <- system.time(lab.result.instances.res <-
                                dbGetQuery(pdsConnection, lab.result.instances.q))
  print(Sys.time())
  print(timed.system)
  
  return(lab.result.instances.res)
  
  
}

get.lab.res.instrances.for.empis.from.file <-
  function(inner.filename) {
    inner.res <- read_csv(inner.filename)
    inner.res <-
      inner.res[, setdiff(colnames(inner.res), unwanted.cols)]
    return(inner.res)
  }

get.lab.orders.from.file <- function(inner.filename) {
  inner.res <- read_csv(inner.filename)
  inner.res <-
    inner.res[, setdiff(colnames(inner.res), unwanted.cols)]
  return(inner.res)
}

get.lab.res.items.for.empis.from.file <- function(inner.filename) {
  inner.res <- read_csv(inner.filename)
  inner.res <-
    inner.res[, setdiff(colnames(inner.res), unwanted.cols)]
  return(inner.res)
}

minimal.spaces <- function(input.vector) {
  input.vector <- gsub(pattern = "^ +",
                       replacement = "",
                       x = input.vector)
  
  input.vector <- gsub(pattern = " +",
                       replacement = " ",
                       x = input.vector)
  
  input.vector <- gsub(pattern = " +$",
                       replacement = "",
                       x = input.vector)
  
  return(input.vector)
}

####

graph.empis <-
  get.empis.from.graphdb(
    inner.endpoint = config$my.graphdb.base,
    inner.repo = patient.scoping.repo,
    inner.sa = saved.authentication
  )

lab.result.instances.res <-
  get.lab.res.instrances.for.empis.from.file(lab.result.instances.res.fp)

roi.res <- get.lab.orders.from.file(roi.res.fp)

rlri.res <-
  get.lab.res.items.for.empis.from.file(rlri.res.fp)

lab.result.instances.res.enriched <-
  left_join(x = lab.result.instances.res, y = roi.res)


# get targent from order if possible
lab.result.instances.res.enriched <-
  left_join(x = lab.result.instances.res.enriched, y = rlri.res)

####

# text mine orders
fortm <-
  lab.result.instances.res.enriched[, c("PK_ORDER_ITEM_ID",
                                        "ORDER_ITEM_CODE",
                                        "ORDER_ITEM_DESCRIPTION")]
fortm$text <-
  paste0(fortm$ORDER_ITEM_CODE, " ", fortm$ORDER_ITEM_DESCRIPTION)
fortm <- fortm[, c("PK_ORDER_ITEM_ID", "text")]

tidy_text <- fortm %>%
  unnest_tokens(word, text)

text_word_counts <- tidy_text %>%
  count(word, PK_ORDER_ITEM_ID, sort = TRUE)

dtm_from_tidy <-
  cast_dtm(
    data = text_word_counts,
    term = word,
    document = PK_ORDER_ITEM_ID,
    value = n
  )

dtm_from_tidy <- as.data.frame.matrix(dtm_from_tidy)
dtm_from_tidy <- dtm_from_tidy[colSums(dtm_from_tidy) > 1]
dtm_from_tidy <- as.matrix(dtm_from_tidy)
dtm_from_tidy[,] <-
  ifelse(dtm_from_tidy %in% c("0", "FALSE"), FALSE, TRUE)
dtm_from_tidy <-
  cbind.data.frame("PK_ORDER_ITEM_ID" = as.numeric(rownames(dtm_from_tidy)), dtm_from_tidy)
dtm_from_tidy <- dtm_from_tidy[,!is.na(colnames(dtm_from_tidy))]

dim(dtm_from_tidy)
# > dim(dtm_from_tidy)
# [1] 339 476

fortm <-
  unique(lab.result.instances.res.enriched[, c("PK_ORDER_ITEM_ID",
                                               "ORDER_ITEM_CODE",
                                               "ORDER_ITEM_DESCRIPTION")])

# additional uniqification unnecessary?
dtm_from_tidy <- unique(inner_join(x = fortm, y = dtm_from_tidy))

write.csv(x = dtm_from_tidy,
          file = "build/order_tm.csv",
          row.names = FALSE)

# text mine specimens
# REFACTOR

fortm <-
  unique(lab.result.instances.res.enriched[, c("SPECIMEN",
                                               "SPECIMEN_COLLECTION_METHOD")])

fortm$text <-
  paste0(fortm$SPECIMEN, " ", fortm$SPECIMEN_COLLECTION_METHOD)

fortm$space.sep.catted.index <-
  paste0(fortm$SPECIMEN, " ", fortm$SPECIMEN_COLLECTION_METHOD)

tidy_text <- fortm %>%
  unnest_tokens(word, text)

text_word_counts <- tidy_text %>%
  count(word, space.sep.catted.index, sort = TRUE)

dtm_from_tidy <-
  cast_dtm(
    data = text_word_counts,
    term = word,
    document = space.sep.catted.index,
    value = n
  )

dtm_from_tidy <- as.data.frame.matrix(dtm_from_tidy)
dtm_from_tidy <- dtm_from_tidy[colSums(dtm_from_tidy) > 1]
dtm_from_tidy <- as.matrix(dtm_from_tidy)
dtm_from_tidy[, ] <-
  ifelse(dtm_from_tidy %in% c("0", "FALSE"), FALSE, TRUE)

dtm_from_tidy <-
  cbind.data.frame("space.sep.catted.index" = rownames(dtm_from_tidy), dtm_from_tidy)

dtm_from_tidy <- dtm_from_tidy[, !is.na(colnames(dtm_from_tidy))]

write.csv(x = dtm_from_tidy,
          file = "build/specimen_tm.csv",
          row.names = FALSE)

# done, just requires recreating space-delimited concatenated
# for joining back with enhanced instances frame


# text mine result items/target entities/analytes
# REFACTOR

fortm <-
  lab.result.instances.res.enriched[, c("PK_LAB_RESULT_ITEM_ID",
                                        "LAB_RESULT_ITEM_CODE",
                                        "LAB_RESULT_ITEM_DESCRIPTION")]
fortm$text <-
  paste0(fortm$LAB_RESULT_ITEM_CODE,
         " ",
         fortm$LAB_RESULT_ITEM_DESCRIPTION)
fortm <- fortm[, c("PK_LAB_RESULT_ITEM_ID", "text")]

rlri.counts <-
  make.table.frame(lab.result.instances.res.enriched$PK_LAB_RESULT_ITEM_ID)

colnames(rlri.counts) <- c("PK_LAB_RESULT_ITEM_ID", "rlri.count")

rlri.counts$PK_LAB_RESULT_ITEM_ID <-
  as.numeric(rlri.counts$PK_LAB_RESULT_ITEM_ID)

fortm <- inner_join(x = fortm, y = rlri.counts)
fortm <- fortm[fortm$rlri.count > 1 ,]


tidy_text <- fortm %>%
  unnest_tokens(word, text)

text_word_counts <- tidy_text %>%
  count(word, PK_LAB_RESULT_ITEM_ID, sort = TRUE)

dtm_from_tidy <-
  cast_dtm(
    data = text_word_counts,
    term = word,
    document = PK_LAB_RESULT_ITEM_ID,
    value = n
  )

dtm_from_tidy <- as.data.frame.matrix(dtm_from_tidy)
dtm_from_tidy <- dtm_from_tidy[colSums(dtm_from_tidy) > 1]
dtm_from_tidy <- as.matrix(dtm_from_tidy)
dtm_from_tidy[, ] <-
  ifelse(dtm_from_tidy %in% c("0", "FALSE"), FALSE, TRUE)
dtm_from_tidy <-
  cbind.data.frame("PK_LAB_RESULT_ITEM_ID" = as.numeric(rownames(dtm_from_tidy)), dtm_from_tidy)
dtm_from_tidy <- dtm_from_tidy[, !is.na(colnames(dtm_from_tidy))]

# > dim(dtm_from_tidy)
# [1] 314 420
# with min 1: 253 385

rlri.counts <-
  make.table.frame(lab.result.instances.res.enriched$PK_LAB_RESULT_ITEM_ID)
colnames(rlri.counts) <- c("PK_LAB_RESULT_ITEM_ID", "rlri.count")
rlri.counts$PK_LAB_RESULT_ITEM_ID <-
  as.numeric(rlri.counts$PK_LAB_RESULT_ITEM_ID)

fortm <-
  lab.result.instances.res.enriched[, c("PK_LAB_RESULT_ITEM_ID",
                                        "LAB_RESULT_ITEM_CODE",
                                        "LAB_RESULT_ITEM_DESCRIPTION")]

fortm <- inner_join(x = fortm, y = rlri.counts)

# additional uniqification unnecessary?
dtm_from_tidy <- unique(inner_join(x = fortm, y = dtm_from_tidy))

write.csv(x = dtm_from_tidy,
          file = "build/result_item_tm.csv",
          row.names = FALSE)
####

obo.ids.labels.q <- "
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
select
distinct
?g ?s ?i ?lab ?lang
where {
    graph ?g {
        ?s a owl:Class .
        # optional {
            ?s rdfs:label ?lab .
            bind(lang(?lab) as ?lang)
        # }
        # optional {
            ?s oboInOwl:id ?i
        # }
    }
}"

obo.ids.labels.res <-
  q2j2df(query = obo.ids.labels.q,
         endpoint = primary.gdb.endpoint,
         repo = primary.gdb.repo)

# OBI and some other ontologies don't use CURIE style IDs
# REFACTOR/generalize
obi.labels.q <- '
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
select
distinct
("http://purl.obolibrary.org/obo/obi.owl" as ?g) ?s ?i ?lab ?lang
where {
    graph <http://purl.obolibrary.org/obo/obi.owl> {
        ?s a owl:Class .
        # optional {
            ?s rdfs:label ?lab .
            bind(lang(?lab) as ?lang)
        # }
    }
}'

obi.labels.res <-
  q2j2df(query = obi.labels.q,
         endpoint = primary.gdb.endpoint,
         repo = primary.gdb.repo)

# create curies
obi.labels.res$i <- obo.iri.to.curie(obi.labels.res$s)

# bind with flexibility regarding presence of columns
obo.labels.res <- bind_rows(obo.ids.labels.res, obi.labels.res)

obo.labels.res.cuire.cols <- split.curie.col(obo.labels.res$i)

obo.labels.res <-
  cbind.data.frame(obo.labels.res, obo.labels.res.cuire.cols)

obo.labels.res$obo.ont <-
  grepl(pattern = obo.iri.prefix,
        x = obo.labels.res$g,
        fixed = TRUE)

print(table(obo.labels.res$obo.ont, useNA = 'a'))

obo.labels.res <- obo.labels.res[obo.labels.res$obo.ont,]

obo.labels.res$onto.prefix <- onto.iri.to.prefix(obo.labels.res$g)

obo.labels.res$owner <-
  tolower(obo.labels.res$term.prefix) == obo.labels.res$onto.prefix

print(table(obo.labels.res$owner, useNA = 'a'))

obo.labels.res <- obo.labels.res[obo.labels.res$owner,]

#### 2020-12-15 11:15

# lab.result.instances <- read_csv(lab.result.instances.fp)
#
# rlri_all_pds_all_obo <- read_csv(rlri_all_pds_all_obo.fp)

rlri_all_pds_all_obo <-
  read_csv("build/annotated/result_item_tm_annotated.csv")

external.shuttle <-
  unique(lab.result.instances.res.enriched[, c("PK_LAB_RESULT_ITEM_ID", "LOINC")])
rlri_all_pds_all_obo <-
  left_join(x = rlri_all_pds_all_obo, y = external.shuttle)

####

obo.to.external <-
  unique(rlri_all_pds_all_obo[, c("LOINC", "result.item.targent.obo")])

obo.to.external <-
  obo.to.external[complete.cases(obo.to.external), ]

promiscuous.externals <- make.table.frame(obo.to.external$LOINC)
promiscuous.externals <-
  promiscuous.externals[promiscuous.externals$count > 1 ,]

obo.to.external <-
  obo.to.external[!(obo.to.external$LOINC %in% promiscuous.externals$value) ,]
#
# rlri.complete <- read_csv(rlri.complete.fp)
# rlri.complete <-
#   rlri.complete[, setdiff(colnames(rlri.complete), unwanted.cols)]
#
# length(unique(lab.result.instances$PK_LAB_RESULT_ITEM_ID))
# length(unique(rlri_all_pds_all_obo$PK_LAB_RESULT_ITEM_ID))

instances.no.analyte.mapping.pks <-
  setdiff(
    lab.result.instances.res$PK_LAB_RESULT_ITEM_ID,
    rlri_all_pds_all_obo$PK_LAB_RESULT_ITEM_ID
  )

discover.analyte.from.external <-
  rlri_all_pds_all_obo[rlri_all_pds_all_obo$PK_LAB_RESULT_ITEM_ID %in% instances.no.analyte.mapping.pks , ]

discover.analyte.from.external <-
  discover.analyte.from.external[!(is.na(discover.analyte.from.external$LOINC)),]

discover.analyte.from.external <-
  inner_join(x = discover.analyte.from.external, y = obo.to.external)

# IF the nrow is > 0, then there are lab result instances for which the obo analyte can be discovered
# transitively trough some external code with previously curated, apparently 1:1 mappings with OBO

####

# roi_all_pds_all_obo <- read_csv(roi_all_pds_all_obo.fp)
roi_all_pds_all_obo <-
  read_csv("build/annotated/order_tm_annotated.csv")

# specimen_with_method_table_annotated_all <-
#   read_csv(specimen_with_method_table_annotated_all.fp)
#
# specimen_with_method_table_annotated_all <-
#   specimen_with_method_table_annotated_all[!is.na(specimen_with_method_table_annotated_all$obo) ,]

specimen_with_method_table_annotated_all <-
  read_csv("build/annotated/specimen_tm_annotated.csv")

####

lab.result.instances.res.enriched <-
  inner_join(x = lab.result.instances.res.enriched, y = rlri_all_pds_all_obo[c(
    "PK_LAB_RESULT_ITEM_ID",
    "result.item.order",
    "result.item.specimen.obo",
    "result.item.specimen.label.override",
    "result.item.targent.obo",
    "result.item.targent.label.override"
  )])

lab.result.instances.res.enriched$space.sep.catted.index <-
  paste0(
    lab.result.instances.res.enriched$SPECIMEN,
    " ",
    lab.result.instances.res.enriched$SPECIMEN_COLLECTION_METHOD
  )

lab.result.instances.res.enriched <-
  left_join(x = lab.result.instances.res.enriched, y = specimen_with_method_table_annotated_all[, c(
    "space.sep.catted.index",
    "specimen.specimen.obo",
    "specimen.specimen.label.override"
  )])

lab.result.instances.res.enriched <-
  left_join(x = lab.result.instances.res.enriched, y = roi_all_pds_all_obo[, c(
    "PK_ORDER_ITEM_ID",
    "order.order",
    "order.specimen.obo",
    "order.specimen.label.override",
    "order.targent.obo",
    "order.targent.label.override"
  )])

####

no.result.item.targ.ent <-
  is.na(lab.result.instances.res.enriched$result.item.targent.obo)

has.order.tag.ent <-
  !is.na(lab.result.instances.res.enriched$result.item.targent.obo)

table(no.result.item.targ.ent, has.order.tag.ent, useNA = 'a')

lab.result.instances.res.enriched$use.order.targent <-
  has.order.tag.ent

lab.result.instances.res.enriched$result.item.targent.obo[has.order.tag.ent] <-
  lab.result.instances.res.enriched$result.item.targent.obo[has.order.tag.ent]

lab.result.instances.res.enriched$result.item.targent.label.override[has.order.tag.ent] <-
  lab.result.instances.res.enriched$result.item.targent.label.override[has.order.tag.ent]

lab.result.instances.res.enriched <-
  lab.result.instances.res.enriched[!(is.na(lab.result.instances.res.enriched$result.item.targent.obo)) , ]


# # should we check if there are already some analyte only assays?
# # how would we know whether the order portion is important?
# already.targent.only <-
#   lab.result.instances.res.enriched[is.na(lab.result.instances.res.enriched$obo.specimen.explicit) &
#                                   is.na(lab.result.instances.res.enriched$obo) ,]
lab.result.instances.res.enriched$targent.only <- FALSE

# do we need to elect the analytes and specimens (order vs. X) first? just for conflicts?
# may already have something like this
analyte.list <-
  unique(lab.result.instances.res.enriched$result.item.targent.obo)

analyte.list <-
  cbind.data.frame(result.item.targent.obo = analyte.list, targent.only = TRUE)


# check this with a left join sometimes.
# not getting a label for chebi globulin
analyte.list <-
  inner_join(
    x = analyte.list,
    y = obo.labels.res[, c("i", "lab")],
    by = c("result.item.targent.obo" = "i")
  )

colnames(analyte.list) <-
  c("result.item.targent.obo",
    "targent.only",
    "result.item.targent.label.override")

lab.result.instances.res.enriched <-
  bind_rows(lab.result.instances.res.enriched, analyte.list)

####

analyte.candidates <-
  sort(unique(lab.result.instances.res.enriched$result.item.targent.obo))
ac.concat <- paste(analyte.candidates, collapse = '" "')
ac.concat <- paste0('"', ac.concat, '"')

mol.ent.q <- paste0(
  "select * where {
    values ?s { " ,
  ac.concat ,
  " }
    ?c <http://www.geneontology.org/formats/oboInOwl#id> ?s ;
       <http://www.w3.org/2000/01/rdf-schema#subClassOf>+ <http://purl.obolibrary.org/obo/CHEBI_23367>
}"
)

mol.ent.res <-
  q2j2df(query = mol.ent.q,
         endpoint = "http://localhost:7200",
         repo = "bottom_up")

lab.result.instances.res.enriched$mol.ent.targent <- FALSE
lab.result.instances.res.enriched$mol.ent.targent[lab.result.instances.res.enriched$result.item.targent.obo %in% mol.ent.res$s] <-
  TRUE

print(table(lab.result.instances.res.enriched$mol.ent.targent, useNA = 'a'))

####

cell.q <- paste0(
  "select * where {
    values ?s { " ,
  ac.concat ,
  " }
    ?c <http://www.geneontology.org/formats/oboInOwl#id> ?s ;
       <http://www.w3.org/2000/01/rdf-schema#subClassOf>* <http://purl.obolibrary.org/obo/CL_0000000>
}"
)

cell.res <-
  q2j2df(query = cell.q,
         endpoint = "http://localhost:7200",
         repo = "bottom_up")

lab.result.instances.res.enriched$cell.targent <- FALSE
lab.result.instances.res.enriched$cell.targent[lab.result.instances.res.enriched$result.item.targent.obo %in% cell.res$s] <-
  TRUE

print(table(lab.result.instances.res.enriched$cell.targent, useNA = 'a'))

####

# http://purl.obolibrary.org/obo/GO_0032991

complex.q <- paste0(
  "select * where {
    values ?s { " ,
  ac.concat ,
  " }
    ?c <http://www.geneontology.org/formats/oboInOwl#id> ?s ;
       <http://www.w3.org/2000/01/rdf-schema#subClassOf>* <http://purl.obolibrary.org/obo/GO_0032991>
}"
)

complex.res <-
  q2j2df(query = complex.q,
         endpoint = "http://localhost:7200",
         repo = "bottom_up")

lab.result.instances.res.enriched$complex.targent <- FALSE
lab.result.instances.res.enriched$complex.targent[lab.result.instances.res.enriched$result.item.targent.obo %in% complex.res$s] <-
  TRUE

print(table(lab.result.instances.res.enriched$complex.targent, useNA = 'a'))

####

# http://purl.obolibrary.org/obo/GO_0032991

property.q <- paste0(
  "select * where {
    values ?s { " ,
  ac.concat ,
  " }
    ?c <http://www.geneontology.org/formats/oboInOwl#id> ?s ;
       <http://www.w3.org/2000/01/rdf-schema#subClassOf>+ <http://purl.obolibrary.org/obo/PATO_0000001>
}"
)

property.res <-
  q2j2df(query = property.q,
         endpoint = "http://localhost:7200",
         repo = "bottom_up")

lab.result.instances.res.enriched$property.targent <- FALSE
lab.result.instances.res.enriched$property.targent[lab.result.instances.res.enriched$result.item.targent.obo %in% property.res$s] <-
  TRUE

print(table(lab.result.instances.res.enriched$property.targent, useNA = 'a'))

####

process.q <- paste0(
  "select * where {
    values ?s { " ,
  ac.concat ,
  " }
    ?c <http://www.geneontology.org/formats/oboInOwl#id> ?s ;
       <http://www.w3.org/2000/01/rdf-schema#subClassOf>+ <http://purl.obolibrary.org/obo/GO_0008150>
}"
)

process.res <-
  q2j2df(query = process.q,
         endpoint = "http://localhost:7200",
         repo = "bottom_up")

lab.result.instances.res.enriched$process.targent <- FALSE
lab.result.instances.res.enriched$process.targent[lab.result.instances.res.enriched$result.item.targent.obo %in% process.res$s] <-
  TRUE

print(table(lab.result.instances.res.enriched$process.targent, useNA = 'a'))

# print(lab.result.instances.res.enriched$result.item.targent.label.override[lab.result.instances.res.enriched$process.targent])

####

temp <-
  lab.result.instances.res.enriched[!(
    lab.result.instances.res.enriched$mol.ent.targent |
      lab.result.instances.res.enriched$cell.targent |
      lab.result.instances.res.enriched$complex.targent |
      lab.result.instances.res.enriched$property.targent |
      lab.result.instances.res.enriched$process.targent
  ) &
    (!(
      is.na(lab.result.instances.res.enriched$result.item.targent.obo)
    )),]

####

lab.result.instances.res.enriched$analyte.compared <-
  compare.to.order(
    lab.result.instances.res.enriched,
    "order.targent.obo",
    "result.item.targent.obo"
  )

# print(table(lab.result.instances.res.enriched$analyte.compared))

temp <- summarize.comparison(
  data = lab.result.instances.res.enriched,
  comparison.col.name = "analyte.compared",
  reportable.levels = c("conflict"),
  reported.columns = c(
    "ORDER_ITEM_CODE",
    "ORDER_ITEM_DESCRIPTION",
    "order.targent.obo",
    "order.targent.label.override",
    "LAB_RESULT_ITEM_CODE",
    "LAB_RESULT_ITEM_DESCRIPTION",
    "result.item.targent.obo",
    "result.item.targent.label.override",
    "analyte.compared"
  )
)

# conflicts between orders and result items
print(temp)


# temp <- summarize.comparison(
#   data = lab.result.instances.res.enriched,
#   comparison.col.name = "analyte.compared",
#   reportable.levels = c("duplication"),
#   reported.columns = c(
#     "ORDER_ITEM_CODE",
#     "ORDER_ITEM_DESCRIPTION",
#     "order.targent.obo",
#     "order.targent.label.override",
#     "LAB_RESULT_ITEM_CODE",
#     "LAB_RESULT_ITEM_DESCRIPTION",
#     "result.item.targent.obo",
#     "result.item.targent.label.override",
#     "analyte.compared"
#   )
# )
#
# # redundancies between orders and result items
# # note the presence of NAs
# print(temp)

lab.result.instances.res.enriched$specimen.compared <-
  compare.to.order(lab.result.instances.res.enriched,
                   "order.specimen.obo",
                   "specimen.specimen.obo")

temp <- summarize.comparison(
  data = lab.result.instances.res.enriched,
  comparison.col.name = "specimen.compared",
  reportable.levels = c("conflict"),
  reported.columns = c(
    "ORDER_ITEM_CODE",
    "ORDER_ITEM_DESCRIPTION",
    "order.specimen.obo",
    "order.specimen.label.override",
    "SPECIMEN",
    "SPECIMEN_COLLECTION_METHOD",
    "specimen.specimen.obo",
    "specimen.specimen.label.override",
    "specimen.compared"
  )
)


print(temp)


###

# result PK doesn't go into model, but keep track of the PK/IRI associations

robot.col.names <- c(
  "ID",
  "label",
  "curation_status",
  "alt_term",
  "def",
  "def_source",
  "usage_example",
  "editor_note",
  "term_editor",
  "term_requestor",
  "term_tracker_id",
  "logical_type",
  "logical_note",
  "parent_class",
  "mat_proc_tech",
  "det_tech",
  "evaluant",
  "analyte",
  "device",
  "reagent",
  "mol_lab",
  "input",
  "output",
  "targ_ent",
  "objective",
  "assoc_axes",
  "dbxref"
)

# any shared column names?
print(intersect(
  colnames(lab.result.instances.res.enriched),
  robot.col.names
))

robot.portion <-
  matrix(
    data = "",
    nrow = nrow(lab.result.instances.res.enriched),
    ncol = length(robot.col.names)
  )
robot.portion <- as.data.frame(robot.portion)
colnames(robot.portion) <- robot.col.names

lab.result.instances.res.enriched <-
  cbind.data.frame(lab.result.instances.res.enriched, robot.portion)

lab.result.instances.res.enriched$term_editor <- term_editor
lab.result.instances.res.enriched$def_source <- def_source
lab.result.instances.res.enriched$term_tracker_id <- term_tracker_id

lab.result.instances.res.enriched$logical_type <- "subclass"
lab.result.instances.res.enriched$logical_type[lab.result.instances.res.enriched$targent.only] <-
  "equivalent"

print(table(lab.result.instances.res.enriched$logical_type, useNA = 'a'))


####

# [29] "mol.ent.targent"                     "cell.targent"                        "complex.targent"                     "property.targent"
# [33] "process.targent"

lab.result.instances.res.enriched$parent_class <- "assay"

lab.result.instances.res.enriched$parent_class[lab.result.instances.res.enriched$mol.ent.targent] <-
  "'analyte assay'"

lab.result.instances.res.enriched$parent_class[lab.result.instances.res.enriched$cell.targent] <-
  "'cytometry assay'"

print(table(lab.result.instances.res.enriched$parent_class, useNA = 'a'))

####

lab.result.instances.res.enriched$dbxref[!(is.na(lab.result.instances.res.enriched$LOINC))] <-
  lab.result.instances.res.enriched$LOINC[!(is.na(lab.result.instances.res.enriched$LOINC))]

# print(table(lab.result.instances.res.enriched$dbxref, useNA = 'a'))

####

# # LESS RELEVANT NOW THAT WE HAVE LABEL OVERRIDE COLUMNS?
# # no, still have to put something in the anayte column, label or ID
#
# lab.result.instances.res.enriched$any.quoted.quote <- FALSE

temp <-
  match(lab.result.instances.res.enriched$result.item.targent.obo,
        obo.labels.res$i)
temp <- obo.labels.res$lab[temp]
temp[is.na(temp)] <- ""

lab.result.instances.res.enriched$analyte[lab.result.instances.res.enriched$mol.ent.targent] <-
  temp[lab.result.instances.res.enriched$mol.ent.targent]
lab.result.instances.res.enriched$analyte[lab.result.instances.res.enriched$analyte == "''"] <-
  ""

sq.flag <-
  grepl(pattern = "'", x = lab.result.instances.res.enriched$analyte)

print(table(sq.flag, useNA = 'a'))
print(lab.result.instances.res.enriched$analyte[sq.flag])

# lab.result.instances.res.enriched$any.quoted.quote <-
#   lab.result.instances.res.enriched$any.quoted.quote + qq.flag
# print(lab.result.instances.res.enriched$analyte[lab.result.instances.res.enriched$any.quoted.quote > 0])
#
# qq.flag <-
#   grepl(pattern = '"', x = lab.result.instances.res.enriched$analyte)
# print(table(qq.flag, useNA = 'a'))
# print(lab.result.instances.res.enriched$analyte[qq.flag])
# lab.result.instances.res.enriched$any.quoted.quote <-
#   lab.result.instances.res.enriched$any.quoted.quote + qq.flag
# print(lab.result.instances.res.enriched$analyte[lab.result.instances.res.enriched$any.quoted.quote > 0])

lab.result.instances.res.enriched$analyte[nchar(lab.result.instances.res.enriched$analyte) > 0 &
                                            !sq.flag] <-
  paste0("'", lab.result.instances.res.enriched$analyte[nchar(lab.result.instances.res.enriched$analyte) > 0 &
                                                          !sq.flag], "'")

lab.result.instances.res.enriched$analyte[nchar(lab.result.instances.res.enriched$analyte) > 0 &
                                            sq.flag] <-
  paste0('"', lab.result.instances.res.enriched$analyte[nchar(lab.result.instances.res.enriched$analyte) > 0 &
                                                          sq.flag], '"')


# would robot tolerate the "''" values?

lab.result.instances.res.enriched$targ_ent[!lab.result.instances.res.enriched$mol.ent.targent] <-
  temp[!lab.result.instances.res.enriched$mol.ent.targent]
lab.result.instances.res.enriched$targ_ent[lab.result.instances.res.enriched$targ_ent == "''"] <-
  ""

# qq.flag <-
#   grepl(pattern = "'", x = lab.result.instances.res.enriched$targ_ent)
# print(table(qq.flag, useNA = 'a'))
# print(lab.result.instances.res.enriched$targ_ent[qq.flag])
# lab.result.instances.res.enriched$any.quoted.quote <-
#   lab.result.instances.res.enriched$any.quoted.quote + qq.flag
# print(lab.result.instances.res.enriched$targ_ent[lab.result.instances.res.enriched$any.quoted.quote > 0])
#
# qq.flag <-
#   grepl(pattern = '"', x = lab.result.instances.res.enriched$targ_ent)
# print(table(qq.flag, useNA = 'a'))
# print(lab.result.instances.res.enriched$targ_ent[qq.flag])
# lab.result.instances.res.enriched$any.quoted.quote <-
#   lab.result.instances.res.enriched$any.quoted.quote + qq.flag
# print(lab.result.instances.res.enriched$targ_ent[lab.result.instances.res.enriched$any.quoted.quote > 0])

lab.result.instances.res.enriched$targ_ent[nchar(lab.result.instances.res.enriched$targ_ent) > 0] <-
  paste0("'", lab.result.instances.res.enriched$targ_ent[nchar(lab.result.instances.res.enriched$targ_ent) > 0], "'")

temp <-
  match(lab.result.instances.res.enriched$specimen.specimen.obo,
        obo.labels.res$i)
temp <- obo.labels.res$lab[temp]
temp[is.na(temp)] <- ""

lab.result.instances.res.enriched$evaluant <- paste0("'", temp, "'")
lab.result.instances.res.enriched$evaluant[lab.result.instances.res.enriched$evaluant == "''"] <-
  ""


# analyte... check "''"s
# targ_ent... check "''"s
# evaluant... all well characterized specimens now that alpha and freetext results are excluded

# todo... get components term labels! may override for assay labels
# label
# def

# # lab.result.instances.res.enriched$authoritative.analyte
#
# # order first
# # could just get rid of tidy orders that only contain an analyte?
# label.working <- ""
# label.working[!is.na(lab.result.instances.res.enriched$tidied)] <-
#   lab.result.instances.res.enriched$tidied[!is.na(lab.result.instances.res.enriched$tidied)]
#
# # tempX <- match(lab.result.instances.res.enriched$)
#
# # then specimen, preferring ordered specimen
# # if the order mentions a specimen,
# # and it is the same as the specimen form the result table
# # don't duplicate
# lab.result.instances.res.enriched$authoritative.labeling.specimen <- ""
# lab.result.instances.res.enriched$authoritative.labeling.specimen[lab.result.instances.res.enriched$specimen.compared == "OK" |
#                                                                 lab.result.instances.res.enriched$specimen.compared == "conflict"] <-
#   lab.result.instances.res.enriched$tidied[!is.na(lab.result.instances.res.enriched$tidied)]
#
#
# temp2 <- ""
# temp2[!is.na(lab.result.instances.res.enriched$tidied)] <-
#   lab.result.instances.res.enriched$tidied[!is.na(lab.result.instances.res.enriched$tidied)]

####


# template.for.robot <-
#   unique(lab.result.instances.res.enriched[lab.result.instances.res.enriched$any.quoted.quote < 1, c(robot.col.names, "specimen.compared", "analyte.compared")])

# try to find rows where the order has a specimen or target entity even though the lab result doesn't have a specimen
#   or the lab result item doesn't map to anything in OBO?

# lab.result.instances.res.enriched$specimen.specimen.obo[lab.result.instances.res.enriched$specimen.compared == "OK" &
#                                                           !is.na(lab.result.instances.res.enriched$order.specimen.obo)]

template.for.robot <-
  as.matrix(unique(lab.result.instances.res.enriched[, c(
    robot.col.names,
    "specimen.compared",
    "analyte.compared",
    "order.order",
    "order.specimen.label.override",
    "order.targent.label.override",
    "result.item.specimen.label.override",
    "result.item.targent.label.override",
    "specimen.specimen.label.override",
    "mol.ent.targent",
    "cell.targent",
    "complex.targent",
    "property.targent",
    "process.targent"
  )]))

template.for.robot[is.na(template.for.robot)] <- ""
template.for.robot <- as.data.frame(template.for.robot)

# # NAs never allowed
# template.for.robot <-
#   template.for.robot[complete.cases(template.for.robot), ]

# template.for.robot <-
#   template.for.robot[template.for.robot$analyte != "''" |
#                        template.for.robot$targ_ent != "''" , ]


template.for.robot <-
  template.for.robot[(template.for.robot$analyte != "''" &
                        template.for.robot$analyte != "") |
                       (template.for.robot$targ_ent != "''" &
                          template.for.robot$targ_ent != ""),]

template.for.robot$def_source <-
  rep(unique(def_source), nrow(template.for.robot))

####

# template.for.robot$label <- template.for.robot$ID

template.for.robot$label <-
  paste(
    template.for.robot$order.order,
    template.for.robot$specimen.specimen.label.override,
    template.for.robot$result.item.targent.label.override,
    "assay"
  )

# fix targent only assays
# I thought I had this under control when the targent only assays were first created

template.for.robot$label <- minimal.spaces(template.for.robot$label)

####

template.for.robot$def <- paste("An assay that measures ",
                                template.for.robot$result.item.targent.label.override)

# logicals got matrixed to strings
template.for.robot$def[template.for.robot$cell.targent == "TRUE"] <-
  paste0(
    "A cytometry assay that measures the count of ",
    template.for.robot$result.item.targent.label.override[template.for.robot$cell.targent == "TRUE"],
    's '
  )

template.for.robot$def[template.for.robot$mol.ent.targent == "TRUE"] <-
  paste0(
    "An analyte  assay that measures the abundance of ",
    template.for.robot$result.item.targent.label.override[template.for.robot$mol.ent.targent == "TRUE"]
  )

# specimen
current.flag <-
  (template.for.robot$specimen.specimen.label.override != "")

template.for.robot$def[current.flag] <-
  paste(
    template.for.robot$def[current.flag],
    "in",
    template.for.robot$specimen.specimen.label.override[current.flag]
  )

# order
current.flag <- (template.for.robot$order.order != "")
template.for.robot$def[current.flag] <-
  paste(template.for.robot$def[current.flag],
        "as part of a",
        template.for.robot$order.order[current.flag],
        "order")

template.for.robot$def <- minimal.spaces(template.for.robot$def)

# template.for.robot$def <- minimal.spaces(template.for.robot$def)

####

# causing duplicate labels?

# template.for.robot <- unique(template.for.robot[, c(robot.col.names,
#                                                     "specimen.compared",
#                                                     "analyte.compared")])

template.for.robot <-
  unique(template.for.robot[, setdiff(robot.col.names, "dbxref")])

label.check <- make.table.frame(template.for.robot$label)

template.for.robot$ID <-
  paste0(id.prefix, id.start + seq(0:(nrow(template.for.robot) - 1)))

dim(template.for.robot)
colnames(template.for.robot)

template.for.robot <- as.matrix(template.for.robot)
colnames(template.for.robot) <- NULL
colnames(template.for.robot)

####

obi_assay_headers <-
  read_delim(
    obi_assay_headers.fp,
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE
  )

obi_assay_headers <- as.matrix(obi_assay_headers)
obi_assay_headers[is.na(obi_assay_headers)] <- ""
colnames(obi_assay_headers) <- NULL
dim(obi_assay_headers)
colnames(obi_assay_headers)

template.for.robot <-
  as.data.frame(rbind(obi_assay_headers, template.for.robot))

if (write.robot.tsv) {
  write_delim(
    x = template.for.robot,
    file = "build/turbo_bottomup_assay_template.tsv",
    delim = "\t",
    col_names = FALSE
  )
}

if (interactively.write.gsheet) {
  gs4_auth(
    email = gsheets.auth.email,
    path = NULL,
    scopes = "https://www.googleapis.com/auth/spreadsheets",
    cache = gargle::gargle_oauth_cache(),
    use_oob = gargle::gargle_oob_default(),
    token = NULL
  )
  
  Sys.sleep(sleepsecs)
  
  # interactive? how to script authentication? OOB?
  # I haven't had much luck naming the workbook in advance... API calls time out
  # fix it manually afterwards?
  googlewb <-
    gs4_create()
  
  Sys.sleep(sleepsecs)
  
  #
  clinical.assays.from.TURBO <- template.for.robot
  colnames(clinical.assays.from.TURBO) <-
    clinical.assays.from.TURBO[1,]
  clinical.assays.from.TURBO <- clinical.assays.from.TURBO[-1,]
  
  # write
  sheet_write(data = clinical.assays.from.TURBO, ss = googlewb)
  
  # try to read the google sheet back into a dataframe
  clinical.assays.from.TURBO <-
    read_sheet(ss = googlewb, sheet = "clinical.assays.from.TURBO")
}

# REFACTOR
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/chebi.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/cl.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/foodon.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/go.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/obi.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/pato.owl", intern = TRUE)
print(result)
result <- system(command = "robot mirror --input-iri http://purl.obolibrary.org/obo/pr.owl", intern = TRUE)
print(result)

result <-
  system(command = 'cd purl.obolibrary.org/obo && export ROBOT_JAVA_ARGS=-8mx4G && robot merge --inputs "*.owl" --output obo_merged_for_turbo_assays.owl', intern = TRUE)
print(result)

result <-
  system(command = 'ROBOT_JAVA_ARGS=-Xmx8G && robot template --template build/turbo_bottomup_assay_template.tsv --output build/turbo_bottomup_assays.ttl --input purl.obolibrary.org/obo/obo_merged_for_turbo_assays.owl --ancestors', intern = TRUE)
print(result)

result <-
  system(command = 'export ROBOT_JAVA_ARGS=-Xmx8G && robot merge --input purl.obolibrary.org/obo/obi.owl --input build/turbo_bottomup_assays.ttl --output build/turbo_bottomup_assays_in_obi.ttl')
print(result)

# ~ 15 minutes
result <-
  system(command = 'export ROBOT_JAVA_ARGS=-Xmx8G && robot reason --reasoner hermit  --input build/turbo_bottomup_assays_in_obi.ttl --output build/turbo_bottomup_assays_in_obi_reasoned.ttl')
print(result)

# ROBOT_JAVA_ARGS=-8mx4G &&
# may not be necessary... did I hard code 24G into the shell script?
# uberon not necessary now that we're using precomposed OBI specimens

