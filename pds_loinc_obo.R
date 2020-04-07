# will probably require an existing VPN and tunnel

# eliminate/postpone loinc challenge parts? look for "^".

# is glomerular filtration rate still in there?

source('pds_loinc_obo_setup.R')

LoincPartLink <-
  read_csv(config$LoincPartLink.file)

Part <-
  read_csv(config$Part.file)

# this is for finding common lab results,
# which have already been annotated with a  LOINC code in PDS

# lexically decomposing and mapping unannotated PDS lab results is a separate task

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

# ~ 300 seconds for ~ 3800 rows, with as few as one unique EMPI

print(Sys.time())
timed.system <- system.time(
  pds.loinc.frequency.resuts <-
    dbGetQuery(pdsConnection, config$pds.loinc.frequency.query)
)
print(Sys.time())
print(timed.system)

# Close connection
dbDisconnect(pdsConnection)

# ~ 3450 r_lab_orders ordered for >= (2) unique EMPIs
common.loincs <-
  pds.loinc.frequency.resuts$LOINC[pds.loinc.frequency.resuts$LOINC_COUNT >= config$common.threshold]

relevant.primary.part.cast <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs &
                  LoincPartLink$LinkTypeName == "Primary", ]

relevant.primary.part.cast <-
  relevant.primary.part.cast[, c("LoincNumber", "PartNumber", "PartTypeName")]

relevant.primary.part.cast <-
  dcast(data = relevant.primary.part.cast,
        formula = LoincNumber ~ PartTypeName,
        value.var = "PartNumber")

# down to ~ 3400

# bad or old LOINCs in PDS?
# setdiff(common.loincs, relevant.primary.part.cast$LoincNumber)
# [1] "11575-5"  "41613"    "5876-7"   "103331-7" "5778-7"   "280602-2" "74217-8"  "39478-6"  "SOLOINC"  "18185-8"  "34695-0"
# [12] "6690-3"   "786-5"    "777-4"    "X29142-7" "UPE18"    "3148-9"   "4544-4"   "789-9"    "0051531"  "788-1"    "787-3"
# [23] "2890-3"   "L4985"    "718-8"    "785-7"    "L4613"

####

pds.with.loinc.parts <-
  left_join(x = pds.loinc.frequency.resuts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

# 3836
# method can legitimately be NA, but when any other column is NA,
# either the PDS count should be below the threshold
# or the LOINC code should eb in the bad/old list above

###

pds_loinc_properties_reviewed <-
  read_csv(config$reviewed.properties.file)
pds_loinc_properties_reviewed$use[is.na(pds_loinc_properties_reviewed$use)] <-
  'n'

# add TURBO datum terms to reviewed spreadsheet

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

# all single UBERON terms except
# LP7576-4 Ser/Plas http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969
# LP7579-8 Ser/Plas/Bld http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969|http://purl.obolibrary.org/obo/UBERON_0000178
# LP7536-8 RBC http://purl.obolibrary.org/obo/CL_0000232
# LP7720-8 WBC http://purl.obolibrary.org/obo/CL_0000738

pds_loinc_systems_reviewed <-
  pds_loinc_systems_reviewed$PartTypeVal[pds_loinc_systems_reviewed$use == 'y']


PartRelatedCodeMapping <-
  read_csv(config$PartRelatedCodeMapping.file)


# see scale and time filtering decisions below
pds.with.loinc.parts <-
  pds.with.loinc.parts[pds.with.loinc.parts$PROPERTY %in% pds_loinc_properties_reviewed &
                         pds.with.loinc.parts$SCALE == "LP7753-9" &
                         pds.with.loinc.parts$SYSTEM %in% pds_loinc_systems_reviewed &
                         (pds.with.loinc.parts$TIME == "LP6960-1" |
                            pds.with.loinc.parts$TIME == "LP6924-7") ,]

# ~ 1800

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

# ~ 1350
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

table(loinc.provided.component.mappings$ExtCodeSystem)

loinc.provided.rxnorm <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem ==
                                      "http://www.nlm.nih.gov/research/umls/rxnorm" , ]

loinc.provided.chebi <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem == "https://www.ebi.ac.uk/chebi" , ]

loinc.provided.rxnorm.only <-
  setdiff(loinc.provided.rxnorm$PartNumber,
          loinc.provided.chebi$PartNumber)
loinc.provided.rxnorm.only <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$PartNumber %in% loinc.provided.rxnorm.only,]

# adalimumab (HUMIRA) http://purl.obolibrary.org/obo/DRON_00018971
# opiates [confirmed not present verbatim in dron or chebi]


# SO JUST USE THESE
loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem  == "https://www.ebi.ac.uk/chebi" &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" , ]


####
# recap: we have made decisive choices about
#   ignoring method for now
#   28 handpicked concenrtation/count etc. properites, which will have to be mapped to TURBO or OBO terms (possibly new)
#   quantitiave scale (OBO/TURBO term  ?)
#   20 handpicked bodily fluid systems... shouldn't be too hard to map to Uberon or Cell Ontology (BioPortal?)
#   PT or 24H time  (OBO/TURBO directive terms =  ?)

# Still need to map components
# keep new successful mappings and current failures seperate

loinc.provided.soemthing <-
  unique(loinc.provided.component.mappings$PartNumber)

# this should get smaller over the course of executing this script
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
    "http://data.bioontology.org/ontologies/RXNORM",
    'http://data.bioontology.org/ontologies/UBERON',
    'http://data.bioontology.org/ontologies/UPHENO',
    'http://data.bioontology.org/ontologies/VO'
  )


acceptable.ontologies <-
  c(
    'http://data.bioontology.org/ontologies/CHEBI',
    'http://data.bioontology.org/ontologies/DRON',
    "http://data.bioontology.org/ontologies/RXNORM"
  )

acceptable.retreived.and.parsed <-
  unique(retreived.and.parsed[retreived.and.parsed$target.ontology %in% acceptable.ontologies , c("source.term", "target.term")])
acceptable.retreived.and.parsed$target.term <-
  as.character(acceptable.retreived.and.parsed$target.term)
acceptable.retreived.and.parsed$source.id <-
  gsub(pattern = "http://purl.bioontology.org/ontology/LNC/",
       replacement = "",
       x = acceptable.retreived.and.parsed$source.term)
acceptable.retreived.and.parsed$source.term <-
  as.character(acceptable.retreived.and.parsed$source.term)

# USE THESE
# retreive labels from loinc and external sides

loinc.still.unmapped.components <-
  setdiff(
    loinc.unmapped.components$PartTypeVal,
    acceptable.retreived.and.parsed$source.id
  )
loinc.still.unmapped.components <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% loinc.still.unmapped.components , ]

loinc.still.unmapped.components$s.end <-
  grepl(pattern = "s$",
        x = loinc.still.unmapped.components$PartName,
        ignore.case = TRUE)


loinc.still.unmapped.components$has.non.alpha <-
  grepl(pattern = "[^a-zA-Z ]",
        x = loinc.still.unmapped.components$PartName,
        ignore.case = TRUE)

potential.cells <-
  loinc.still.unmapped.components[loinc.still.unmapped.components$s.end &
                                    !(loinc.still.unmapped.components$has.non.alpha),
                                  c("PartTypeVal", "PartDisplayName")]
potential.cells <-
  potential.cells[order(potential.cells$PartDisplayName), ]

# modify the input so that it's easoer to merge abck in with LOINC codes

ols.attempts <-
  apply(
    X = potential.cells,
    MARGIN = 1,
    FUN = function(current.row) {
      #ols.attempts <- lapply(sort(potential.cells), function(current.potential) {
      current.potential <- tolower(current.row[["PartDisplayName"]])
      current.potential <-
        sub(pattern = "s$",
            replacement = "",
            x = current.potential)
      singular.lc <- current.potential
      current.potential <-
        gsub(pattern = " ",
             replacement = ",",
             x = current.potential)
      print(current.potential)
      
      ols.attempt <-
        httr::GET(
          paste0(
            "https://www.ebi.ac.uk/ols/api/search?q={",
            current.potential,
            "}&type=class&local=true&ontology=uberon,cl,pr,chebi&rows=9"
          )
        )
      ols.attempt <- ols.attempt$content
      ols.attempt <- rawToChar(ols.attempt)
      ols.attempt <- jsonlite::fromJSON(ols.attempt)
      ols.attempt <- ols.attempt$response$docs
      if (is.data.frame(ols.attempt)) {
        ols.attempt$query <- singular.lc
        ols.attempt$loinc.part <- current.row[["PartTypeVal"]]
        return(ols.attempt)
      }
    }
  )

ols.attempts <- do.call(rbind.data.frame, ols.attempts)
ols.successes <-
  ols.attempts[ols.attempts$label == ols.attempts$query ,]

# how to pick which hits are acceptable?
# lowercased query with trialing s removed should be an exact match for lowercased label?
# &queryFields={label,synonym}
# &ontology=

loinc.still.unmapped.components <-
  loinc.still.unmapped.components[(!(
    loinc.still.unmapped.components$PartTypeVal %in% ols.successes$loinc.part
  )), ]

# PR considers LP15333-5 "Alanine aminotransferase" underspecified...
# should contain a species and numerical subtype like "human Alanine aminotransferase 1"
# try to find other LOINC protein components
# subclassof "Enzymes" http://purl.bioontology.org/ontology/LNC/LP31392-1 <- CHAL
#  <- Chamistry and Chemistry challenge <- Lab <- LOINCPARTS

# LP6118-6 Albumin sco Protein http://purl.bioontology.org/ontology/LNC/LP15838-3
#  <- UA <- Lab ...

# LP17698-9 Erythrocyte distribution width
#   variability of etrythorcyte sizes... a quality of the distribution/population

# LP30809-5 Anion gap... difference between concentration of common cations and common anions

# LP14492-0 Urea nitrogen: overspecified as far as chebi is concerned
# LP15441-6 bicarbonate: underspecified/synonymy (hydrogencarbonate anion ChEBI 17544)

# LP286653-3 Neutrophils/100 leukocytes

# write.csv(loinc.still.unmapped.components, "loinc_components_obo-unmappable_by_loinc_bioportal_ols.csv")


# find common ancestors of unmapped LOINC components
# rdflib, rrdf, (sparql... against what endpoint), igraph (from loinc tabular files)

MultiAxialHierarchy <-
  read_csv("~/Loinc_2.67/AccessoryFiles/MultiAxialHierarchy/MultiAxialHierarchy.csv")

MultiAxialHierarchy.unmapped <-
  MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% loinc.still.unmapped.components$PartTypeVal , ]

split.paths <-
  strsplit(MultiAxialHierarchy.unmapped$PATH_TO_ROOT,
           split = ".",
           fixed = TRUE)

unlisted.paths <- unlist(split.paths)

node.appearances <- table(unlisted.paths)

node.appearances <-
  cbind.data.frame(names(node.appearances), as.numeric(node.appearances))

colnames(node.appearances) <- c("part", "appearances")

node.appearances$part <- as.character(node.appearances$part)

common.nodes <-
  node.appearances[node.appearances$appearances >= 5 , ]

common.nodes <-
  left_join(common.nodes, Part, by = c("part" = "PartNumber"))

common.nodes <- common.nodes[complete.cases(common.nodes),]

####

# lots of creatinine normalized chemicals, aas

####

# refactor

# chem.unmapped <-
#   grepl(pattern = "LP7786-9", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# chem.unmapped <- MultiAxialHierarchy.unmapped[chem.unmapped,]
# chem.unmapped$non.parent <- "LP7786-9"
#
# aa.unmapped <-
#   grepl(pattern = "LP18033-8", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# aa.unmapped <- MultiAxialHierarchy.unmapped[aa.unmapped,]
# aa.unmapped$non.parent <- "LP18033-8"
#
# drug.tox.unmapped <-
#   grepl(pattern = "LP7790-1", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# drug.tox.unmapped <-
#   MultiAxialHierarchy.unmapped[drug.tox.unmapped,]
# drug.tox.unmapped$non.parent <- "LP18033-8"
#
# drug.unmapped <-
#   grepl(pattern = "LP18046-0", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# drug.unmapped <- MultiAxialHierarchy.unmapped[drug.unmapped,]
# drug.unmapped$non.parent <- "LP18046-0"
#
# fa.unmapped <-
#   grepl(pattern = "LP14488-8", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# fa.unmapped <- MultiAxialHierarchy.unmapped[fa.unmapped,]
# fa.unmapped$non.parent <- "LP14488-8"
#
# electro.sing.val.unmapped <-
#   grepl(pattern = "LP19403-2", x = MultiAxialHierarchy.unmapped$PATH_TO_ROOT)
# electro.sing.val.unmapped <-
#   MultiAxialHierarchy.unmapped[electro.sing.val.unmapped,]
# electro.sing.val.unmapped$non.parent <- "LP19403-2"
#
# creatinine.opportunities <-
#   rbind.data.frame(
#     electro.sing.val.unmapped,
#     fa.unmapped,
#     drug.unmapped,
#     drug.tox.unmapped,
#     aa.unmapped,
#     chem.unmapped
#   )

creatinine.opportunities <-
  loinc.still.unmapped.components[grepl(pattern = "/Creatinine$", loinc.still.unmapped.components$PartName),]

print(length(unique(creatinine.opportunities$PartTypeVal)))

write.csv(creatinine.opportunities,
          "loinc_creatinine_normalized_components.csv")

creatinine.opportunities$analyte <-
  sub(pattern = "/Creatinine$",
      replacement = "",
      x = creatinine.opportunities$PartName)

# refactor, see almost identical code above for cell types
ols.attempts <-
  apply(
    X = creatinine.opportunities,
    MARGIN = 1,
    FUN = function(current.row) {
      #ols.attempts <- lapply(sort(potential.cells), function(current.potential) {
      current.potential <- tolower(current.row[["analyte"]])
      current.potential <-
        sub(pattern = "s$",
            replacement = "",
            x = current.potential)
      singular.lc <- current.potential
      current.potential <-
        gsub(pattern = " ",
             replacement = ",",
             x = current.potential)
      print(current.potential)
      
      ols.attempt <-
        httr::GET(
          paste0(
            "https://www.ebi.ac.uk/ols/api/search?q={",
            current.potential,
            "}&type=class&local=true&ontology=uberon,cl,pr,chebi&rows=9"
          )
        )
      ols.attempt <- ols.attempt$content
      ols.attempt <- rawToChar(ols.attempt)
      ols.attempt <- jsonlite::fromJSON(ols.attempt)
      ols.attempt <- ols.attempt$response$docs
      if (is.data.frame(ols.attempt)) {
        if (ncol(ols.attempt) == 10) {
          ols.attempt$query <- singular.lc
          ols.attempt$loinc.part <- current.row[["PartTypeVal"]]
          return(ols.attempt)
        }
        
      }
    }
  )

ols.no.creatinine <- do.call(rbind.data.frame, ols.attempts)
no.creatinine.successes <-
  ols.no.creatinine[ols.no.creatinine$label == ols.no.creatinine$query ,]

# USE THESE (from ChEBI)

###

tokens.from.unmapped <-
  gsub("[[:punct:]]",
       " ",
       tolower(loinc.still.unmapped.components$PartName))

tokens.from.unmapped <-
  termFreq(tokens.from.unmapped)

tokens.from.unmapped <-
  cbind.data.frame(names(tokens.from.unmapped),
                   as.numeric(tokens.from.unmapped))

names(tokens.from.unmapped) <- c("token", "count")
tokens.from.unmapped <-
  tokens.from.unmapped[tokens.from.unmapped$count >= 5 ,]

# ~ 300 rows in loinc.still.unmapped.components matching 'Ab.Ig'
# ~ 50 "Ab$"
# ~ 20 'Ag'
# ~ 20 "/100 leukocytes"
# ~ 30 "Cells."
# "^^" adjustment
# look all of these patterns up in LoincPartLink for the different semantic contexts (LinkPartName text and "Property" URI)

# ~ 100 '/Creatinine' see above
# ~ 10 "+" cases using drug or chemical names, not cell types
# "^" could indicate a challenge, but see also 4x "^peak" and 4x "^trough"
# use chebi appings except for

uses.plus.sign <-
  loinc.still.unmapped.components[grepl(pattern = "\\+", x = loinc.still.unmapped.components$PartName) &
                                    (!(
                                      grepl(pattern = "Cells", x = loinc.still.unmapped.components$PartName)
                                    )), ]

uses.plus.sign$left <-
  tolower(sub(pattern = "\\+.*$", "", uses.plus.sign$PartName))
uses.plus.sign$right <-
  tolower(sub(pattern = "^.*\\+", "", uses.plus.sign$PartName))

uses.plus.sign.queries <-
  unique(sort(c(
    uses.plus.sign$left, uses.plus.sign$right
  )))

# refactor, see almost identical code above for cell types
ols.attempts <-
  lapply(uses.plus.sign.queries, function(current.row) {
    #ols.attempts <- lapply(sort(potential.cells), function(current.potential) {
    current.potential <- current.row
    current.potential <-
      gsub(pattern = " ",
           replacement = ",",
           x = current.potential)
    print(current.potential)
    
    ols.attempt <-
      httr::GET(
        paste0(
          "https://www.ebi.ac.uk/ols/api/search?q={",
          current.potential,
          "}&type=class&local=true&rows=9&ontology=chebi,dron,rxnorm&rows=9"
        )
      )
    ols.attempt <- ols.attempt$content
    ols.attempt <- rawToChar(ols.attempt)
    ols.attempt <- jsonlite::fromJSON(ols.attempt)
    ols.attempt <- ols.attempt$response$docs
    if (is.data.frame(ols.attempt)) {
      if (ncol(ols.attempt) == 10) {
        ols.attempt$query <- current.row
        # ols.attempt$loinc.part <- current.row[["PartTypeVal"]]
        return(ols.attempt)
      }
      
    }
  })

ols.uses.plus <- do.call(rbind.data.frame, ols.attempts)
ols.uses.plus$label <- tolower(ols.uses.plus$label)
ols.uses.plus$query <- tolower(ols.uses.plus$query)
uses.plus.successes <-
  ols.uses.plus[ols.uses.plus$label == ols.uses.plus$query ,]

uses.plus.sign$left <- tolower(uses.plus.sign$left)
uses.plus.sign$right <- tolower(uses.plus.sign$right)

uses.plus.sign.remerge <-
  left_join(
    uses.plus.sign,
    uses.plus.successes,
    by = c("left" = "label"),
    suffix = c("", ".left")
  )
uses.plus.sign.remerge <-
  uses.plus.sign.remerge[, setdiff(colnames(uses.plus.sign.remerge), "description")]

uses.plus.sign.remerge <-
  left_join(
    uses.plus.sign.remerge,
    uses.plus.successes,
    by = c("right" = "label"),
    suffix = c("", ".right")
  )

uses.plus.sign.remerge <-
  uses.plus.sign.remerge[, setdiff(colnames(uses.plus.sign.remerge), "description")]

uses.plus.sign.remerge <-
  uses.plus.sign.remerge[complete.cases(uses.plus.sign.remerge), ]

setdiff(
  uses.plus.successes$label,
  union(uses.plus.sign.remerge$left, uses.plus.sign.remerge$right)
)


# vitamin d2
# chebi
# CHEBI
# class
# TRUE
# A vitamin D supplement and has been isolated from alfalfa.
# calciferol
#
# isovaleryl-l-carnitine
# chebi
# CHEBI
# class
# TRUE
# An O-isovalerylcarnitine that is the 3-methylbutanoyl (isovaleryl) derivative of L-carnitine.
# isovalerylcarnitine
#
# desmethyldoxepin
# chebi
# CHEBI
# class
# TRUE
# A dibenzooxepine resulting from the demethylation of the antidepressant doxepin. It is the active metabolite of doxepin.
# nordoxepin
#
# 2-methylbutyrylcarnitine
# chebi
# CHEBI
# class
# TRUE
# A C5-acylcarnitine having 2-methylbutyryl as the acyl substituent.
# methylbutyrylcarnitine (c5)
#
# --- doesn't look promsing below ---
#
# norclomipramine
#
# methsuximide
# normethsuximide

# USE THESE
# as long as left and right match

write.csv(uses.plus.sign, "loinc_plus_sign_components.csv")

###   ###   ###

# ^ pk components

###   ###   ###

any.antibody <-
  loinc.still.unmapped.components[grepl(pattern = "Ab$", x = loinc.still.unmapped.components$PartName),]
