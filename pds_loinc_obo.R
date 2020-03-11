# will probably require an existing VPN and tunnel

# eliminate/postpone loinc challenge parts? look for "^".

# is glomerular filtration rate still in there?

library(config)
library(RJDBC)
library(readr)
library(reshape2)
# see also https://jangorecki.gitlab.io/data.cube/library/data.table/html/dcast.data.table.html
library(dplyr)
library(httr)

# find driver from config file
config <- config::get(file = "loinc_in_obo.yaml")

LoincPartLink <-
  read_csv(config$LoincPartLink.file)

Part <-
  read_csv(config$Part.file)

# this is for finding common lab results,
# which have already been annotated with a  LOINC code by PDS

# decomposing unannotated PDS lab results is a separate task

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

pdsConnection <-
  dbConnect(pdsDriver,
            pds.con.string,
            config$pds.user,
            config$pds.pw)

pds.loinc.frequency.query <- "
SELECT
RLRI.LOINC,
COUNT(1) as loinc_count
FROM
mdm.LAB_RESULTS lr
JOIN mdm.R_LAB_RESULT_ITEM rlri ON
lr.FK_LAB_RESULT_ITEM_ID = RLRI.PK_LAB_RESULT_ITEM_ID
WHERE
RLRI.LOINC IS NOT NULL
GROUP BY
RLRI.LOINC"

# ~ 300 seconds for ~ 4000 rows

print(Sys.time())
timed.system <- system.time(
  pds.loinc.frequency.resuts <-
    dbGetQuery(pdsConnection, pds.loinc.frequency.query)
)
print(Sys.time())
print(timed.system)

# Close connection
dbDisconnect(pdsConnection)

common.loincs <-
  pds.loinc.frequency.resuts$LOINC[pds.loinc.frequency.resuts$LOINC_COUNT >= 2]

relevant.primary.part.cast <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs &
                  LoincPartLink$LinkTypeName == "Primary", ]

relevant.primary.part.cast <-
  relevant.primary.part.cast[, c("LoincNumber", "PartNumber", "PartTypeName")]

relevant.primary.part.cast <-
  dcast(data = relevant.primary.part.cast,
        formula = LoincNumber ~ PartTypeName,
        value.var = "PartNumber")

###

pds.with.loinc.parts <-
  left_join(x = pds.loinc.frequency.resuts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

###

pds_loinc_properties_reviewed <-
  read_csv(config$reviewed.properties.file)
pds_loinc_properties_reviewed$use[is.na(pds_loinc_properties_reviewed$use)] <-
  'n'
pds_loinc_properties_reviewed <-
  pds_loinc_properties_reviewed$PartTypeVal[pds_loinc_properties_reviewed$use == 'y']

# manually curate LOINC system-OBO term mappings
# other options: bioportal mappings or term lookups,
# OLS,
# local Solr
# ontobee API?
# loinc provided ChEBI mappings

pds_loinc_systems_reviewed <- read_csv(config$reviewed.systems.file)
pds_loinc_systems_reviewed$use[is.na(pds_loinc_systems_reviewed$use)] <-
  'n'
pds_loinc_systems_reviewed <-
  pds_loinc_systems_reviewed$PartTypeVal[pds_loinc_systems_reviewed$use == 'y']


PartRelatedCodeMapping <-
  read_csv(config$PartRelatedCodeMapping.file)


# see below
pds.with.loinc.parts <-
  pds.with.loinc.parts[pds.with.loinc.parts$SCALE == "LP7753-9" &
                         (pds.with.loinc.parts$TIME == "LP6960-1" |
                            pds.with.loinc.parts$TIME == "LP6924-7") &
                         pds.with.loinc.parts$PROPERTY %in% pds_loinc_properties_reviewed &
                         pds.with.loinc.parts$SYSTEM %in% pds_loinc_systems_reviewed ,]

pds.with.loinc.parts.melt <-
  melt(
    data = pds.with.loinc.parts,
    measure.vars = c("COMPONENT", "METHOD", "PROPERTY",
                     "SCALE", "SYSTEM", "TIME"),
    id.vars = c("LOINC", "LOINC_COUNT")
  )

pds.prominent.loinc.parts <-
  aggregate(
    pds.with.loinc.parts.melt$LOINC_COUNT,
    by = list(
      pds.with.loinc.parts.melt$variable,
      pds.with.loinc.parts.melt$value
    ),
    FUN = mean,
    na.rm = TRUE
  )

colnames(pds.prominent.loinc.parts) <-
  c("PartTypeName", "PartTypeVal", "pds.count")

pds.prominent.loinc.parts <-
  left_join(x = pds.prominent.loinc.parts,
            y = Part,
            by = c("PartTypeVal" = "PartNumber"))

# temp <- pds.with.loinc.parts

# use Qn SCALE only (5x more entires than next most commons, Nom and Ord)
# Pt TIME MUCH more common than next time (24H). XXX and * more common then 24H.
# ignore METHOD for now

# LP7753-9
# SCALE
# Qn
# Quantitative
#
# LP6960-1
# TIME
# Pt
# Point in time (spot)
#
# LP6924-7
# TIME
# 24H
# 24 hours

# component

# write.csv(pds.prominent.loinc.parts[pds.priminent.loinc.parts$PartTypeName == "PROPERTY",], "pds_loinc_properties.csv")

needs.component.mappings <- unique(pds.with.loinc.parts$COMPONENT)
loinc.provided.component.mappings <-
  PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% needs.component.mappings , ]

table(
  loinc.provided.component.mappings$ExtCodeSystem,
  loinc.provided.component.mappings$Equivalence
)

loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem %in%
                                      c("https://www.ebi.ac.uk/chebi",
                                        "http://www.nlm.nih.gov/research/umls/rxnorm") &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" , ]

loinc.provided.soemthing <-
  unique(loinc.provided.component.mappings$PartNumber)
loinc.unmapped.components <-
  setdiff(needs.component.mappings, loinc.provided.soemthing)
loinc.unmapped.components <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% loinc.unmapped.components , ]

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

retreived.and.parsed <-
  retreive.and.parse(sort(unique(loinc.unmapped.components$PartTypeVal)))

retreived.and.parsed <-
  do.call(rbind.data.frame, retreived.and.parsed)

rap.tab <- table(retreived.and.parsed$target.ontology)
rap.tab <- cbind.data.frame(names(rap.tab), as.numeric(rap.tab))
names(rap.tab) <- c("target.ontology", "map.count")

acceptable.ontologies <-
  c(
    'http://data.bioontology.org/ontologies/CHEBI',
    'http://data.bioontology.org/ontologies/CL',
    'http://data.bioontology.org/ontologies/CLO',
    'http://data.bioontology.org/ontologies/DRON',
    'http://data.bioontology.org/ontologies/FMA',
    'http://data.bioontology.org/ontologies/GO-PLUS',
    'http://data.bioontology.org/ontologies/PR',
    'http://data.bioontology.org/ontologies/UBERON',
    'http://data.bioontology.org/ontologies/UPHENO',
    'http://data.bioontology.org/ontologies/VO'
  )

acceptable.retreived.and.parsed <-
  unique(retreived.and.parsed[retreived.and.parsed$target.ontology %in% acceptable.ontologies , c("source.term", "target.term")])
acceptable.retreived.and.parsed$target.term <-
  as.character(acceptable.retreived.and.parsed$target.term)
acceptable.retreived.and.parsed$source.id <-
  gsub(pattern = "http://purl.bioontology.org/ontology/LNC/",
       replacement = "",
       x = acceptable.retreived.and.parsed$source.term)


loinc.still.unmapped.components <-
  setdiff(
    loinc.unmapped.components$PartTypeVal,
    acceptable.retreived.and.parsed$source.id
  )
loinc.still.unmapped.components <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% loinc.still.unmapped.components , ]
