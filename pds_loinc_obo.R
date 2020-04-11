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
                            pds.with.loinc.parts$TIME == "LP6924-7") ,]

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

# ehr_loinc_counts: LOINC assay codes with counts of mentions in EHR

# LoincPartLink: details about each part. possibly many rows per LOINC assay/part combination
# head(LoincPartLink)
# LoincNumber LongCommonName              PartNumber PartName                 PartCodeSystem   PartTypeName LinkTypeName Property
# <chr>       <chr>                       <chr>      <chr>                    <chr>            <chr>        <chr>        <chr>
#   1 10000-8     R wave duration in lead AVR LP31088-5  R wave duration.lead AVR http://loinc.org COMPONENT    Primary      http://loinc.org/property/COMPONENT
# 2 10000-8     R wave duration in lead AVR LP6879-3   Time                     http://loinc.org PROPERTY     Primary      http://loinc.org/property/PROPERTY
# 3 10000-8     R wave duration in lead AVR LP6960-1   Pt                       http://loinc.org TIME         Primary      http://loinc.org/property/TIME_ASPCT
# 4 10000-8     R wave duration in lead AVR LP7289-4   Heart                    http://loinc.org SYSTEM       Primary      http://loinc.org/property/SYSTEM
# 5 10000-8     R wave duration in lead AVR LP7753-9   Qn                       http://loinc.org SCALE        Primary      http://loinc.org/property/SCALE_TYP
# 6 10000-8     R wave duration in lead AVR LP6244-0   EKG                      http://loinc.org METHOD       Primary      http://loinc.org/property/METHOD_TYP

# head(MultiAxialHierarchy)
# PATH_TO_ROOT                                      SEQUENCE IMMEDIATE_PARENT CODE       CODE_TEXT
# <chr>                                                <dbl> <chr>            <chr>      <chr>
# 1 NA                                                       1 NA               LP29693-6  Laboratory
# 2 LP29693-6                                                1 LP29693-6        LP343406-7 Microbiology and Antimicrobial susceptibility
# 3 LP29693-6.LP343406-7                                     1 LP343406-7       LP7819-8   Microbiology
# 4 LP29693-6.LP343406-7.LP7819-8                            1 LP7819-8         LP14559-6  Microorganism
# 5 LP29693-6.LP343406-7.LP7819-8.LP14559-6                  1 LP14559-6        LP98185-9  Bacteria
# 6 LP29693-6.LP343406-7.LP7819-8.LP14559-6.LP98185-9        1 LP98185-9        LP14082-9  Bacteria

# > head(Part)
# PartNumber PartTypeName PartName                     PartDisplayName              Status
# <chr>      <chr>        <chr>                        <chr>                        <chr>
#   1 LP101394-7 ADJUSTMENT   adjusted for maternal weight adjusted for maternal weight ACTIVE
# 2 LP101907-6 ADJUSTMENT   corrected for age            corrected for age            ACTIVE
# 3 LP115711-6 ADJUSTMENT   corrected for background     corrected for background     ACTIVE
# 4 LP147359-6 ADJUSTMENT   adjusted for body weight     adjusted for body weight     ACTIVE
# 5 LP173482-3 ADJUSTMENT   1st specimen                 1st specimen                 DEPRECATED
# 6 LP173483-1 ADJUSTMENT   post cyanocobalamin          post cyanocobalamin          ACTIVE

# PartRelatedCodeMapping
# external code mappings like ChEBI, RxNorm, taxonomy
# largely from SNOMED: systems, method, property, adjustment, time, scale
# genes (Entrez), multiple radlex fields, challenge (SNOMED or RxNorm)

# these are the counts for all of LOINC, not the assays mentioned within any EHR/CDW like PDS
# sort(table(PartRelatedCodeMapping$ExtCodeSystem[PartRelatedCodeMapping$PartTypeName == 'COMPONENT']))
#
# https://www.ncbi.nlm.nih.gov/clinvar             http://pubchem.ncbi.nlm.nih.gov           https://www.ncbi.nlm.nih.gov/gene
# 85                                         345                                         391
# https://www.ncbi.nlm.nih.gov/taxonomy http://www.nlm.nih.gov/research/umls/rxnorm                 https://www.ebi.ac.uk/chebi
# 831                                         970                                        2298
# http://snomed.info/sct
# 6027

# THIS is a tabulation of parts within the current EHR/CDW
# all part types are included
# > head(pds.prominent.loinc.parts)
# PartTypeName PartTypeVal pds.count                           PartName                            PartDisplayName Status
# 1    COMPONENT  LP100041-5      2344                         Lacosamide                                 Lacosamide ACTIVE
# 2    COMPONENT  LP100616-4         3                 1-Hydroxymidazolam                         1-Hydroxymidazolam ACTIVE
# 3    COMPONENT  LP102314-4       459                  PM-SCL-100 Ab.IgG                             PM-SCL-100 IgG ACTIVE
# 4       METHOD  LP102323-5    138620                                BCG Bromocresol green (BCG) dye binding method ACTIVE
# 5    COMPONENT  LP102618-8       415       Aquaporin 4 water channel Ab               Aquaporin 4 water channel Ab ACTIVE
# 6       METHOD  LP105134-3       114 Creatinine-based formula (CKD-EPI)         Creatinine-based formula (CKD-EPI) ACTIVE

# pds.with.loinc.parts: one row per assay/result type
# columns for the assay/result LOINC code and each LOINC parts code
# head(pds.with.loinc.parts)
# LOINC LOINC_COUNT COMPONENT   METHOD PROPERTY    SCALE   SYSTEM     TIME
# 2  6768-6     6615139 LP15346-7     <NA> LP6789-4 LP7753-9 LP7576-4 LP6960-1
# 4 10834-0      375012 LP14885-5 LP6165-7 LP6827-2 LP7753-9 LP7567-3 LP6960-1
# 5  2345-7    11723946 LP14635-4     <NA> LP6827-2 LP7753-9 LP7576-4 LP6960-1
# 6  3948-7       10807 LP14729-5     <NA> LP6827-2 LP7753-9 LP7576-4 LP6960-1
# 7  5693-7         764 LP14349-2     <NA> LP6827-2 LP7753-9 LP7576-4 LP6960-1
# 9 20661-5         295 LP15946-4     <NA> LP6860-3 LP7753-9 LP7576-4 LP6960-1

#### components, esp. mapping

# PS how did I make those pds_loinc_*_reviewed files above?
# write.csv(pds.prominent.loinc.parts[pds.priminent.loinc.parts$PartTypeName == "PROPERTY",], "pds_loinc_properties.csv")

# keep new successful mappings and current failures separate
# TO-DO: use consistent names

# ~ 1350

original.needs.component.mapping <-
  unique(pds.with.loinc.parts$COMPONENT)

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

# SO JUST USE ChEBI
loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem  == "https://www.ebi.ac.uk/chebi" &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" ,]

loinc.provided.component.mappings$ext.iri <-
  loinc.provided.component.mappings$ExtCodeId
# assume we're just using ChEBI for now
loinc.provided.component.mappings$ext.iri <-
  sub(pattern = '^CHEBI:',
      replacement = 'http://purl.obolibrary.org/obo/CHEBI_',
      x = loinc.provided.component.mappings$ext.iri)

update.accounting(
  loinc.provided.component.mappings,
  "PartNumber",
  "ext.iri",
  "LOINC provided ChEBI mappings"
)

#### begin BioPortal mappings retrieval
# ~ 5 minutes vs remote BioPortal

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

acceptable.retreived.and.parsed <-
  acceptable.retreived.and.parsed[, c("source.id", "source.term", "target.term")]

update.accounting(
  acceptable.retreived.and.parsed,
  "source.id",
  "target.term",
  "BioPortal IRI mappings"
)


####

# starting string searches now
# use core analytes and searchable components (separately)
# with and without trailing s
# strip match labels of 'atom$'
# or do cosine similarity
# ols?
# BioPortal?

# ~ 300 rows in loinc.still.unmapped.components matching 'Ab.Ig'
# ~ 50 "Ab$" antibodies
# ~ 20 'Ag'
# ~ 20 "/100 leukocytes"
# ~ 30 "Cells."
# "^^" adjustment

# see also IgX, Ag

# # ~ 100 '/Creatinine' see above
# # ~ 10 "+" cases using drug or chemical names, not cell types
# # "^" could indicate a challenge, but see also 4x "^peak" and 4x "^trough"
# # use ChEBI mappings except for... ?

# V frame of "CORE" unmapped component terms
# this should be added to the components that don't have CORE analytes
# and should be considered in the context of DIVISOR terms
temp <-
  pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping ,]
temp <- sort(unique(temp$PartTypeVal))
temp <- LoincPartLink[LoincPartLink$PartNumber %in% temp , ]
temp <-
  temp[temp$Property == 'http://loinc.org/property/analyte-core' , ]

names.and.numbers <- unique(temp[, c("PartNumber", "PartName")])

common.core.analyte.numbers <- table(temp$PartNumber)
common.core.analyte.numbers <-
  cbind.data.frame(names(common.core.analyte.numbers),
                   as.numeric(common.core.analyte.numbers))
colnames(common.core.analyte.numbers) <- c('PartNumber', 'count')
common.core.analyte.numbers$PartNumber <-
  as.character(common.core.analyte.numbers$PartNumber)

names.and.numbers <-
  base::merge(names.and.numbers, common.core.analyte.numbers)

# sums of free/bound forms, prodrug and metabolite, etc.
#   Cells.markers
# what were my "plus sign remerges" below for?
# "+" 30/170 core PDS core analytes

names.and.numbers$plus <-
  grepl(pattern = '+',
        x = names.and.numbers$PartName,
        fixed = TRUE)


# challenges
# pharmacokinetics (peak, trough) are considered challenges!
# NONE for these core analytes
# see also ^^
names.and.numbers$caret <-
  grepl(pattern = '^',
        x = names.and.numbers$PartName,
        fixed = TRUE)

# "." 17/170 core PDS core analytes
#   Anatomical/tissue location
#   Form/carrier/state
#   Cells.markers
names.and.numbers$period <-
  grepl(pattern = '.',
        x = names.and.numbers$PartName,
        fixed = TRUE)

# / very rare in core analytes.... that's the whole point
# but look for it in the "general" component terms
# Creatinine
# 100 leukocytes
# Hemoglobin.total
# 100 cells
# A few more
names.and.numbers$ratio <-
  grepl(pattern = '/',
        x = names.and.numbers$PartName,
        fixed = TRUE)

# "cells" 21/170 PDS core analytes
# Mostly Cells.markers
names.and.numbers$cells <-
  grepl(pattern = 'cell',
        x = names.and.numbers$PartName,
        ignore.case = TRUE)

# s$ + " ": almost all cells
# remember, CL does not use a final 's'
# along similar lines, ChEBI INCLUDES a final 'atom' for calcium, sodium, etc.
names.and.numbers$ends.s <-
  grepl(pattern = 's$',
        x = names.and.numbers$PartName,
        ignore.case = TRUE)

# ? mostly plus or period HERE
names.and.numbers$non.alpha <-
  grepl(pattern = '[^a-zA-Z ]',
        x = names.and.numbers$PartName)

names.and.numbers$whitespace <-
  grepl(pattern = ' ',
        x = names.and.numbers$PartName)

#### don't start searching without including
# system names and display names
# non-core names and display names
# divisor names

common.divisors <-
  LoincPartLink[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC &
                  LoincPartLink$Property == 'http://loinc.org/property/analyte-divisor' , ]

cd.names.numbers <-
  unique(common.divisors[, c("PartNumber", "PartName")])

common.divisors <-
  table(common.divisors$PartNumber)

common.divisors <-
  cbind.data.frame(names(common.divisors), as.numeric(common.divisors))

colnames(common.divisors) <- c('PartNumber', 'count')

common.divisors <- base::merge(cd.names.numbers, common.divisors)

# why bother
# common.divisors <- common.divisors[common.divisors$count > 2 , ]

####

# also loook for searchable components

# # COMPONENT
# # Search
# # http://loinc.org/property/search

####

# some more perspective/prioritization regarding
# what kind of patterns to look for

# why not start with all component parts or even all part terms
# as opposed to "salvage mode" below

# MultiAxialHierarchy.unmapped <-
#   MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% current.needs.component.mapping , ]

# MultiAxialHierarchy.unmapped <-
#   MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% original.needs.component.mapping , ]

MultiAxialHierarchy.relevant <-
  MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% pds.prominent.loinc.parts$PartTypeVal ,]

split.paths <-
  strsplit(MultiAxialHierarchy.relevant$PATH_TO_ROOT,
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
  node.appearances[node.appearances$appearances >= 2 ,]

common.nodes <-
  left_join(common.nodes, Part, by = c("part" = "PartNumber"))

common.nodes <- common.nodes[complete.cases(common.nodes), ]

####

# previously looked at common super classes
# now look at common tokens

# temp <-
#   pds.prominent.loinc.parts[pds.prominent.loinc.parts$PartTypeVal %in% current.needs.component.mapping ,]

# # no repeats
# temp <- table(pds.prominent.loinc.parts$PartName)
# temp <- cbind.data.frame(names(temp), as.numeric(temp))

# tabulated.tokens <-
#   gsub("[[:punct:]]",
#        " ",
#        tolower(temp$PartName))

tabulated.tokens <-
  gsub("[[:punct:]]",
       " ",
       tolower(pds.prominent.loinc.parts$PartName))

tabulated.tokens <-
  termFreq(tabulated.tokens)

tabulated.tokens <-
  cbind.data.frame(names(tabulated.tokens),
                   as.numeric(tabulated.tokens))

names(tabulated.tokens) <- c("token", "count")

# add to config
tabulated.tokens <-
  tabulated.tokens[tabulated.tokens$count > 3 ,]


####

all.part.types.in.source <- pds.with.loinc.parts$LOINC
all.part.types.in.source <-
  unique(LoincPartLink$PartNumber[LoincPartLink$LoincNumber %in% all.part.types.in.source])
all.part.types.in.source <-
  Part$PartTypeName[Part$PartNumber %in% all.part.types.in.source]

# # not necessary in our current task
# all.part.types.in.source <- all.part.types.in.source[!grepl(pattern = '^Rad\\.', x = all.part.types.in.source)]
# all.part.types.in.source <- all.part.types.in.source[!grepl(pattern = '^Document\\.', x = all.part.types.in.source)]

all.part.types.in.source <- table(all.part.types.in.source)
all.part.types.in.source <-
  cbind.data.frame(names(all.part.types.in.source),
                   as.numeric(all.part.types.in.source))
names(all.part.types.in.source) <- c('part.type', 'source.count')

# look for low0hanging fruit patterns
followup <-
  LoincPartLink[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC ,]

# # "LoincNumber", "LongCommonName",
# desperate <- unique(followup[, c(
#   "PartNumber",
#   "PartName",
#   "PartCodeSystem",
#   "PartTypeName",
#   "LinkTypeName",
#   "Property"
# )])

text.search.input <- unique(followup[, c("PartNumber",
                                         "PartName")])

colnames(text.search.input) <- c('PartNumber', 'original.name')

get.display.name <-
  Part[Part$PartNumber %in% text.search.input$PartNumber, c('PartNumber', 'PartDisplayName')]

colnames(get.display.name) <- c('PartNumber', 'original.name')

text.search.input <-
  unique(rbind.data.frame(text.search.input, get.display.name))

text.search.input$singularized <- sub(pattern = 's$',
                                      replacement = '',
                                      x = text.search.input$original.name)

original <- text.search.input[, c('PartNumber', 'original.name')]
singularized <- text.search.input[, c('PartNumber', 'singularized')]
colnames(singularized) <- colnames(original)
text.search.input <-
  unique(rbind.data.frame(original, singularized))
colnames(text.search.input) <- c('PartNumber', 'any.name')
print(nrow(text.search.input))

indices <- sort(unique(text.search.input$any.name))
print(length(indices))

# not enough of a difference to get fancy

text.search.input <-
  text.search.input[order(text.search.input$any.name),]

ols.attempts <-
  apply(
    X = text.search.input[1:99, ],
    MARGIN = 1,
    FUN = function(current.row) {
      temp <-
        ols.serch.term.labels.universal(
          current.string = current.row[['any.name']],
          current.id = current.row[['PartNumber']],
          strip.final.s = FALSE,
          # ontology.filter = "&ontology=cl",
          ontology.filter = "",
          kept.row.count = 5
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
# lowercased query with trialing s removed should be an exact match for lowercased label? (In cell ontology CL)
# or try both ways
# or check cosine similarity
# see also trailing 'atom' in ChEBI
# &queryFields={label,synonym}
# &ontology=

ols.successes <-
  ols.attempts[ols.attempts$label == ols.attempts$query ,]

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
# look for double-mapped LOINC parts
print(sort(table(current.component.mapping.frame[, 'source.id'])))

#### 2020 04 08 11 15

# PR considers LP15333-5 "Alanine aminotransferase" underspecified...
# should contain a species and numerical subtype like "human Alanine aminotransferase 1"
# try to find other LOINC protein components
# subclass of "Enzymes" http://purl.bioontology.org/ontology/LNC/LP31392-1 <- CHAL
#  <- Chemistry and Chemistry challenge <- Lab <- LOINCPARTS

# LP6118-6 Albumin subclass of Protein http://purl.bioontology.org/ontology/LNC/LP15838-3
#  <- UA <- Lab ...

# LP17698-9 Erythrocyte distribution width
#   variability of erythrocyte sizes... a quality of the distribution/population

# LP30809-5 Anion gap... difference between concentration of common cations and common anions

# LP14492-0 Urea nitrogen: over specified as far as ChEBI is concerned
# LP15441-6 bicarbonate: underspecified/synonymy (hydrogencarbonate anion ChEBI 17544)

# LP286653-3 Neutrophils/100 leukocytes

# find common ancestors of unmapped LOINC components
# rdflib, rrdf, (SPARQL... against what endpoint), igraph (from LOINC tabular files)

####


# log of terms split of of X+Y

# # vitamin d2
# # ChEBI
# # CHEBI
# # class
# # TRUE
# # A vitamin D supplement and has been isolated from alfalfa.
# # calciferol
# #
# # isovaleryl-l-carnitine
# # ChEBI
# # CHEBI
# # class
# # TRUE
# # An O-isovalerylcarnitine that is the 3-methylbutanoyl (isovaleryl) derivative of L-carnitine.
# # isovalerylcarnitine
# #
# # desmethyldoxepin
# # ChEBI
# # CHEBI
# # class
# # TRUE
# # A dibenzooxepine resulting from the demethylation of the antidepressant doxepin. It is the active metabolite of doxepin.
# # nordoxepin
# #
# # 2-methylbutyrylcarnitine
# # ChEBI
# # CHEBI
# # class
# # TRUE
# # A C5-acylcarnitine having 2-methylbutyryl as the acyl substituent.
# # methylbutyrylcarnitine (c5)
# #
# # --- doesn't look promising below ---
# #
# # norclomipramine
# #
# # methsuximide
# # normethsuximide
