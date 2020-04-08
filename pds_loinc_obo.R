# currently eliminating hard-coded expectation of querying PDS for common LOINC code usage

# idea:  eliminate/postpone LOINC challenge parts? look for "^".

# check: is glomerular filtration rate still in there?

source('pds_loinc_obo_setup.R')

#### CSV reads
### some from LOINC, some manually curated

# long-formatted relations between LOINC assay/order/result numbers and parts
#   "PartNumber"     "PartName" "PartTypeName"
# with parts mentioned multiple times with different contexts
#   "LinkTypeName"   "Property"
LoincPartLink <-
  read_csv(config$LoincPartLink.file)

# additional part details
Part <-
  read_csv(config$Part.file)

# get LOINC provided component mappings
PartRelatedCodeMapping <-
  read_csv(config$PartRelatedCodeMapping.file)

MultiAxialHierarchy <-
  read_csv(config$MultiAxialHierarchy.file)

####

# retrieve common lab results,
# which have already been annotated with a LOINC codes
# lexically decomposing and mapping unannotated PDS lab results is a separate task

# see count_loinc_annotated_pds_results.R for an example

ehr_loinc_counts <-
  as.data.frame(read_csv(config$loinc.count.readpath))

###
# add TURBO datum terms to reviewed spreadsheet
pds_loinc_properties_reviewed <-
  read_csv(config$reviewed.properties.file)
pds_loinc_properties_reviewed$use[is.na(pds_loinc_properties_reviewed$use)] <-
  'n'

pds_loinc_properties_reviewed <-
  pds_loinc_properties_reviewed$PartTypeVal[pds_loinc_properties_reviewed$use == 'y']

##

# manually curate LOINC system-OBO term mappings
# other options: BioPortal mappings or term lookups,
# OLS,
# local Solr
# OntoBee API?
# LOINC provided ChEBI mappings


# all single UBERON terms except
# LP7576-4 Ser/Plas http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969
# LP7579-8 Ser/Plas/Bld http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969|http://purl.obolibrary.org/obo/UBERON_0000178
# LP7536-8 RBC http://purl.obolibrary.org/obo/CL_0000232
# LP7720-8 WBC http://purl.obolibrary.org/obo/CL_0000738


pds_loinc_systems_reviewed <- read_csv(config$reviewed.systems.file)
pds_loinc_systems_reviewed$use[is.na(pds_loinc_systems_reviewed$use)] <-
  'n'
pds_loinc_systems_reviewed <-
  pds_loinc_systems_reviewed$PartTypeVal[pds_loinc_systems_reviewed$use == 'y']

####

# ~ 3450 r_lab_orders ordered for >= (2) unique EMPIs
common.loincs <-
  ehr_loinc_counts$LOINC[ehr_loinc_counts$LOINC_COUNT >= config$common.threshold]

#### put these in config?
# get part values for the "common", LOINC-annotated PDS lab results
relevant.primary.part.cast <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs &
                  LoincPartLink$LinkTypeName == "Primary",]

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
  left_join(x = ehr_loinc_counts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

# 3836
# method can legitimately be NA, but when any other column is NA,
# either the PDS count should be below the threshold
# or the LOINC code should eb in the bad/old list above

###

# make some manual decisions

# PT or 24H time

# scale choice:
#   use Qn SCALE only (5x more entries than next most commons, Nom and Ord)

# LP7753-9
# SCALE
# Qn
# Quantitative

# time choice:
#   Pt TIME MUCH more common than next time (24H). XXX and * more common than 24H.

# LP6960-1
# TIME
# Pt
# Point in time (spot)
#
# LP6924-7
# TIME
# 24H
# 24 hours

# ignore METHOD for now

####


# start merging in reviewed time, system, scale and property choices
# see csv reads at top of script
# put these in config?
pds.with.loinc.parts <-
  pds.with.loinc.parts[pds.with.loinc.parts$PROPERTY %in% pds_loinc_properties_reviewed &
                         pds.with.loinc.parts$SCALE == "LP7753-9" &
                         pds.with.loinc.parts$SYSTEM %in% pds_loinc_systems_reviewed &
                         (pds.with.loinc.parts$TIME == "LP6960-1" |
                            pds.with.loinc.parts$TIME == "LP6924-7") , ]

# ~ 1800

pds.with.loinc.parts.melt <-
  melt(
    data = pds.with.loinc.parts,
    measure.vars = c("COMPONENT", "METHOD", "PROPERTY",
                     "SCALE", "SYSTEM", "TIME"),
    id.vars = c("LOINC", "LOINC_COUNT")
  )

# why does mean aggregation get some fractional count results...
#   must be some rows that are duplicated
pds.prominent.loinc.parts <-
  aggregate(
    pds.with.loinc.parts.melt$LOINC_COUNT,
    by = list(
      pds.with.loinc.parts.melt$variable,
      pds.with.loinc.parts.melt$value
    ),
    FUN = sum,
    na.rm = TRUE
  )

colnames(pds.prominent.loinc.parts) <-
  c("PartTypeName", "PartTypeVal", "pds.count")

pds.prominent.loinc.parts$PartTypeName <-
  as.character(pds.prominent.loinc.parts$PartTypeName)

pds.prominent.loinc.parts <-
  left_join(
    x = pds.prominent.loinc.parts,
    y = Part,
    by = c("PartTypeVal" = "PartNumber", "PartTypeName" = "PartTypeName")
  )

####
# recap: we have made decisive choices about
#   ignoring method for now
#   28 handpicked concentration/count etc. properties, which will have to be mapped to TURBO or OBO terms (possibly new)
#   quantitative  scale (OBO/TURBO term  ?)
#   20 handpicked bodily fluid systems... shouldn't be too hard to map to Uberon or Cell Ontology (BioPortal?)
#   PT or 24H time  (OBO/TURBO directive terms =  ?)

#### components, esp. mapping

# PS how did I make those pds_loinc_*_reviewed files above?
# write.csv(pds.prominent.loinc.parts[pds.priminent.loinc.parts$PartTypeName == "PROPERTY",], "pds_loinc_properties.csv")

# keep new successful mappings and current failures separate
# TO-DO: use consistent names

# ~ 1350
original.needs.component.mapping <-
  unique(pds.with.loinc.parts$COMPONENT)
# # moot
# original.component.mapping.complete <- c()

current.needs.component.mapping <- original.needs.component.mapping
current.component.mapping.complete <- c()
current.component.mapping.frame <- matrix(ncol = 4)
colnames(current.component.mapping.frame) <-
  c('source.id', 'source.term', 'target.term', 'authority')

####

# LOINC asserts some mappings of their own

loinc.provided.component.mappings <-
  PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% current.needs.component.mapping ,]

table(
  loinc.provided.component.mappings$ExtCodeSystem,
  loinc.provided.component.mappings$Equivalence
)

loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem %in%
                                      c("https://www.ebi.ac.uk/chebi",
                                        "http://www.nlm.nih.gov/research/umls/rxnorm") &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" ,]

table(loinc.provided.component.mappings$ExtCodeSystem)

loinc.provided.rxnorm <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem ==
                                      "http://www.nlm.nih.gov/research/umls/rxnorm" ,]

loinc.provided.chebi <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem == "https://www.ebi.ac.uk/chebi" ,]

loinc.provided.rxnorm.only <-
  setdiff(loinc.provided.rxnorm$PartNumber,
          loinc.provided.chebi$PartNumber)
loinc.provided.rxnorm.only <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$PartNumber %in% loinc.provided.rxnorm.only, ]

print(loinc.provided.rxnorm.only)

# adalimumab (HUMIRA) http://purl.obolibrary.org/obo/DRON_00018971
# opiates [confirmed not present verbatim in DrOn or ChEBI]

# SO JUST USE THESE
loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem  == "https://www.ebi.ac.uk/chebi" &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" ,]

# accounting
current.component.mapping.complete <-
  union(
    current.component.mapping.complete,
    loinc.provided.component.mappings$PartNumber
  )
current.needs.component.mapping <-
  setdiff(current.needs.component.mapping,
          loinc.provided.component.mappings$PartNumber)

dput(colnames(loinc.provided.component.mappings))
# c("PartNumber", "PartName", "PartTypeName", "ExtCodeId", "ExtCodeDisplayName",
#   "ExtCodeSystem", "Equivalence", "ContentOrigin", "ExtCodeSystemVersion",
#   "ExtCodeSystemCopyrightNotice")

source.id <- loinc.provided.component.mappings$PartNumber
source.term <-
  paste0("http://purl.bioontology.org/ontology/LNC/", source.id)
# assume we're just using ChEBI for now
target.term <- loinc.provided.component.mappings$ExtCodeId
target.term <-
  sub(pattern = '^CHEBI:',
      replacement = 'http://purl.obolibrary.org/obo/CHEBI_',
      x = target.term)
authority <- 'LOINC provided ChEBI mappings'

temp <-
  cbind(source.id, source.term, target.term, authority)

current.component.mapping.frame <-
  rbind(current.component.mapping.frame, temp)

#### begin BioPortal mappings retrieval
# ~ 5 minutes vs remote BioPortal

# loinc.unmapped.components <-
#   pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping , ]

# is retreive.and.parse making more assumptions about globals than I realized/remembered?
retreived.and.parsed <-
  retreive.and.parse(sort(unique(current.needs.component.mapping)))

retreived.and.parsed <-
  do.call(rbind.data.frame, retreived.and.parsed)

# who did we get mapped to?
rap.tab <- table(retreived.and.parsed$target.ontology)
rap.tab <- cbind.data.frame(names(rap.tab), as.numeric(rap.tab))
names(rap.tab) <- c("target.ontology", "map.count")

acceptable.ontologies <- config$acceptable.map.onts.from.loinc

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
acceptable.retreived.and.parsed$authority <- 'BioPortal mappings'
acceptable.retreived.and.parsed <-
  acceptable.retreived.and.parsed[, c("source.id", "source.term", "target.term", "authority")]

current.component.mapping.complete <-
  union(current.component.mapping.complete,
        acceptable.retreived.and.parsed$source.id)
current.needs.component.mapping <-
  setdiff(current.needs.component.mapping,
          acceptable.retreived.and.parsed$source.id)

# SAF FALSE NOT HELPING... swtiched to matrix
current.component.mapping.frame <-
  rbind(current.component.mapping.frame ,
        as.matrix(acceptable.retreived.and.parsed))


####

# look for cell type mappings that BioPortal misses due to pluralization mismatches

# temp <-
#   setdiff(
#     loinc.unmapped.components$PartTypeVal,
#     acceptable.retreived.and.parsed$source.id
#   )

temp <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping , ]

temp$s.end <-
  grepl(pattern = "s$",
        x = temp$PartName,
        ignore.case = TRUE)


temp$has.non.alpha <-
  grepl(pattern = "[^a-zA-Z ]",
        x = temp$PartName,
        ignore.case = TRUE)

potential.cells <-
  temp[temp$s.end &
         !(temp$has.non.alpha),
       c("PartTypeVal", "PartDisplayName")]
potential.cells <-
  potential.cells[order(potential.cells$PartDisplayName), ]

# modify the input so that it's easier to merge back in with LOINC codes

ols.attempts <-
  apply(
    X = potential.cells,
    MARGIN = 1,
    FUN = function(current.row) {
      temp <-
        ols.serch.term.labels.universal(
          current.string = current.row[['PartDisplayName']],
          current.id = current.row[['PartTypeVal']],
          strip.final.s = TRUE,
          ontology.filter = "&ontology=cl",
          kept.row.count = 9
        )
      # print(is.data.frame(temp))
      # print(nrow(temp))
      if (is.data.frame(temp)) {
        if (nrow(temp) > 0) {
          # temp <- temp[, setdiff(colnames(temp), "description")]
          return(temp)
        }
      }
    }
  )

ols.attempts <- do.call(rbind.data.frame, ols.attempts)

# how to pick which hits are acceptable?
# lowercased query with trialing s removed should be an exact match for lowercased label?
# &queryFields={label,synonym}
# &ontology=

ols.successes <-
  ols.attempts[ols.attempts$label == ols.attempts$query , ]

ols.successes <-
  ols.successes[, setdiff(colnames(ols.successes), "description")]

# c("source.id", "source.term", "target.term", "authority")
source.id <- ols.successes$loinc.part
source.term <-
  paste0("http://purl.bioontology.org/ontology/LNC/", source.id)
# assume we're just using ChEBI for now
target.term <- ols.successes$iri
authority <- 'OLS cells'

temp <-
  cbind(source.id, source.term, target.term, authority)

current.component.mapping.frame <-
  rbind(current.component.mapping.frame, temp)

current.component.mapping.complete <-
  union(current.component.mapping.complete,
        ols.successes$loinc.part)
current.needs.component.mapping <-
  setdiff(current.needs.component.mapping,
          ols.successes$loinc.part)

# keep track of mapping sources
print(sort(table(current.component.mapping.frame[, 'authority'])))
# look for double-ammped LOINC parts
print(sort(table(current.component.mapping.frame[, 'source.id'])))

#### 2020 04 08 11 15

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

# find common ancestors of unmapped LOINC components
# rdflib, rrdf, (sparql... against what endpoint), igraph (from loinc tabular files)

MultiAxialHierarchy.unmapped <-
  MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% current.needs.component.mapping ,]

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

# add to config file
common.nodes <-
  node.appearances[node.appearances$appearances >= 5 ,]

common.nodes <-
  left_join(common.nodes, Part, by = c("part" = "PartNumber"))

common.nodes <- common.nodes[complete.cases(common.nodes), ]

####

# use common class nodes above to extract unmapped terms and look for patterns by eye
# => lots of creatinine normalized chemicals, amino acids, etc.

####

# # exploratory. refactor if needed to run again.
#
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
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping , ]

creatinine.opportunities <-
  creatinine.opportunities[grepl(pattern = '/Creatinine$', x = creatinine.opportunities$PartName) , ]

print(length(unique(creatinine.opportunities$PartTypeVal)))

# write.csv(creatinine.opportunities,
#           "loinc_creatinine_normalized_components.csv")

creatinine.opportunities$analyte <-
  sub(pattern = "/Creatinine$",
      replacement = "",
      x = creatinine.opportunities$PartName)

ols.attempts <-
  apply(
    X = creatinine.opportunities,
    MARGIN = 1,
    FUN = function(current.row) {
      temp <-
        ols.serch.term.labels.universal(
          current.string = current.row[['analyte']],
          current.id = current.row[['PartTypeVal']],
          strip.final.s = TRUE,
          ontology.filter = "&ontology=chebi",
          kept.row.count = 9
        )
      # print(is.data.frame(temp))
      # print(nrow(temp))
      if (is.data.frame(temp)) {
        if (nrow(temp) > 0) {
          temp <- temp[, setdiff(colnames(temp), "description")]
          return(temp)
        }
      }
    }
  )

# should probably state which columns are expected
#   and or do a column-mismatch toleratng bind
# temp <- sapply(ols.attempts, dim)
# temp <- sapply(ols.attempts, is.null)

ols.attempts <- do.call(rbind.data.frame, ols.attempts)

ols.attempts$atom.ok <-
  sub(pattern = ' atom$',
      replacement = '',
      x = ols.attempts$label)

no.creatinine.successes <-
  ols.attempts[ols.attempts$atom.ok == ols.attempts$query , ]

# c("source.id", "source.term", "target.term", "authority")
source.id <- no.creatinine.successes$loinc.part
source.term <-
  paste0("http://purl.bioontology.org/ontology/LNC/", source.id)
# assume we're just using ChEBI for now
target.term <- no.creatinine.successes$iri
authority <- 'OLS creatinine normalized chemicals'

temp <-
  cbind(source.id, source.term, target.term, authority)

current.component.mapping.frame <-
  rbind(current.component.mapping.frame, temp)

current.component.mapping.complete <-
  union(current.component.mapping.complete,
        no.creatinine.successes$loinc.part)
current.needs.component.mapping <-
  setdiff(current.needs.component.mapping,
          no.creatinine.successes$loinc.part)

# keep track of mapping sources
print(sort(table(current.component.mapping.frame[, 'authority'])))
# look for double-ammped LOINC parts
print(sort(table(current.component.mapping.frame[, 'source.id'])))

# diagnostics
temp <-
  creatinine.opportunities[!(creatinine.opportunities$PartTypeVal %in% no.creatinine.successes$loinc.part) ,]

###

# previously lookoed at common superclasses
# now look at common tokens

temp <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping , ]

tokens.from.unmapped <-
  gsub("[[:punct:]]",
       " ",
       tolower(temp$PartName))

tokens.from.unmapped <-
  termFreq(tokens.from.unmapped)

tokens.from.unmapped <-
  cbind.data.frame(names(tokens.from.unmapped),
                   as.numeric(tokens.from.unmapped))

names(tokens.from.unmapped) <- c("token", "count")
tokens.from.unmapped <-
  tokens.from.unmapped[tokens.from.unmapped$count >= 5 , ]

# ~ 300 rows in loinc.still.unmapped.components matching 'Ab.Ig'
# ~ 50 "Ab$"
# ~ 20 'Ag'
# ~ 20 "/100 leukocytes"
# ~ 30 "Cells."
# "^^" adjustment

# look all of these patterns up in LoincPartLink for the different semantic contexts (LinkPartName text and "Property" URI)
# 100 leukocytes
# http://loinc.org
# DIVISOR
# SyntaxEnhancement
# http://loinc.org/property/analyte-divisor



# > dput(colnames(LoincPartLink))
# c("LoincNumber", "LongCommonName", "PartNumber", "PartName",
#   "PartCodeSystem", "PartTypeName", "LinkTypeName", "Property")
# 10327-5	Eosinophils/100 leukocytes in Sputum by Manual count	LP14539-8	Eosinophils	http://loinc.org	COMPONENT	SyntaxEnhancement	http://loinc.org/property/analyte-core
# 10327-5	Eosinophils/100 leukocytes in Sputum by Manual count	LP31960-5	100 leukocytes	http://loinc.org	DIVISOR	SyntaxEnhancement	http://loinc.org/property/analyte-divisor

relevant.links <-
  LoincPartLink[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC ,]
common.divisors <-
  table(relevant.links$PartName[relevant.links$Property == 'http://loinc.org/property/analyte-divisor'])
common.divisors <-
  cbind.data.frame(names(common.divisors), as.numeric(common.divisors))
colnames(common.divisors) <- c('PartName', 'count')
common.divisors <- common.divisors[common.divisors$count > 2 ,]

common.core.analyte.names <-
  table(relevant.links$PartName[relevant.links$Property == 'http://loinc.org/property/analyte-core'])
common.core.analyte.names <-
  cbind.data.frame(names(common.core.analyte.names),
                   as.numeric(common.core.analyte.names))
colnames(common.core.analyte.names) <- c('PartName', 'count')
# common.core.analyte.names <- common.core.analyte.names[common.core.analyte.names$count > 2 , ]
# START by serching these
# searchable components below, too?
# bioportal, ols, or ... ?


# these should have already been looked up via bioportal (right?)
common.core.analyte.codes <-
  table(relevant.links$PartNumber[relevant.links$Property == 'http://loinc.org/property/analyte-core'])
common.core.analyte.codes <-
  cbind.data.frame(names(common.core.analyte.codes),
                   as.numeric(common.core.analyte.codes))
colnames(common.core.analyte.codes) <- c('PartName', 'count')
# common.core.analyte.codes <- common.core.analyte.codes[common.core.analyte.codes$count > 2 , ]

# Creatinine can also be a divisor... what does that leave?

# ~ 100 '/Creatinine' see above
# ~ 10 "+" cases using drug or chemical names, not cell types
# "^" could indicate a challenge, but see also 4x "^peak" and 4x "^trough"
# use chebi mappings except for... ?


####

temp <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping , ]

uses.plus.sign <-
  temp[grepl(pattern = "\\+", x = temp$PartName) &
         (!(grepl(
           pattern = "Cells", x = temp$PartName
         ))), ]


## see also

one.word.searchable.components <-
  LoincPartLink[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC &
                  LoincPartLink$Property == 'http://loinc.org/property/search' &
                  LoincPartLink$LinkTypeName == 'Search' &
                  LoincPartLink$PartTypeName == 'COMPONENT'  ,]
one.word.searchable.components$nospaces <-
  gsub(pattern =  ' +',
       replacement = '',
       x = one.word.searchable.components$PartName)
one.word.searchable.components$spacediff <-
  nchar(one.word.searchable.components$PartName) - nchar(one.word.searchable.components$nospaces)
one.word.searchable.components <-
  one.word.searchable.components[one.word.searchable.components$spacediff == 0 , ]
one.word.searchable.components <-
  table(one.word.searchable.components$LoincNumber)
one.word.searchable.components <-
  cbind.data.frame(
    names(one.word.searchable.components),
    as.numeric(one.word.searchable.components)
  )
colnames(one.word.searchable.components) <-
  c('LOINC.term', 'component.count')
hist(one.word.searchable.components$component.count)

#### 2020 04 08 5 27

# plus.but.not.cells

# COMPONENT
# Search
# http://loinc.org/property/search

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
  ols.uses.plus[ols.uses.plus$label == ols.uses.plus$query , ]

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
  uses.plus.sign.remerge[complete.cases(uses.plus.sign.remerge),]

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
  loinc.still.unmapped.components[grepl(pattern = "Ab$", x = loinc.still.unmapped.components$PartName), ]
