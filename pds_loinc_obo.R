source('turbo_loinc_obo_setup.R')

### not using component display names any more, just core analyte?
### add systmes including display names (blood) back in

#### PUT FILENAME IN CONFIG
harmonized_ontology_rankings <-
  read.csv("~/loinc_in_obo/harmonized_ontology_rankings.csv",
           stringsAsFactors = FALSE)

### CSV reads
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

###   ###   ###

# retrieve common lab results,
# which have already been annotated with a LOINC codes
# lexically decomposing and mapping unannotated PDS lab results is a separate task

# see count_loinc_annotated_pds_results.R for an example

ehr_loinc_counts <-
  as.data.frame(read_csv(config$loinc.count.readpath))

###   ###   ###
# add TURBO datum terms to reviewed spreadsheet

turbo_loinc_properties_reviewed <-
  read_csv(config$reviewed.properties.file)

turbo_loinc_properties_reviewed$use[is.na(turbo_loinc_properties_reviewed$use)] <-
  'n'

turbo_loinc_properties_reviewed <-
  turbo_loinc_properties_reviewed[turbo_loinc_properties_reviewed$use == 'y' , ]

###   ###   ###

# manually curate LOINC system-OBO term mappings
# other options search: BioPortal mappings or term lookups,
# OLS,
# local Solr
# OntoBee API?
# LOINC provided ChEBI mappings

# all single UBERON terms except
# LP7576-4 Ser/Plas http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969
# LP7579-8 Ser/Plas/Bld http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969|http://purl.obolibrary.org/obo/UBERON_0000178
# LP7536-8 RBC http://purl.obolibrary.org/obo/CL_0000232
# LP7720-8 WBC http://purl.obolibrary.org/obo/CL_0000738

# bootstrapped from bulk.systems below
# but that's not working anymore
# could easily be reactivated
turbo_loinc_systems_reviewed <-
  read_csv(config$reviewed.systems.file)
turbo_loinc_systems_reviewed$use[is.na(turbo_loinc_systems_reviewed$use)] <-
  'n'

turbo_loinc_systems_reviewed <-
  turbo_loinc_systems_reviewed[turbo_loinc_systems_reviewed$use == 'y' ,]

###   ###   ###

print(nrow(ehr_loinc_counts))
common.loincs <-
  ehr_loinc_counts$LOINC[ehr_loinc_counts$LOINC_COUNT >= config$common.threshold]
# ~ 3450 r_lab_orders ordered for >= (2) unique EMPIs
print(length(common.loincs))

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
print(nrow(relevant.primary.part.cast))

# bad or old LOINCs in PDS?
# setdiff(common.loincs, relevant.primary.part.cast$LoincNumber)
# [1] "11575-5"  "41613"    "5876-7"   "103331-7" "5778-7"   "280602-2" "74217-8"  "39478-6"  "SOLOINC"  "18185-8"  "34695-0"
# [12] "6690-3"   "786-5"    "777-4"    "X29142-7" "UPE18"    "3148-9"   "4544-4"   "789-9"    "0051531"  "788-1"    "787-3"
# [23] "2890-3"   "L4985"    "718-8"    "785-7"    "L4613"

###   ###   ###

pds.with.loinc.parts <-
  left_join(x = ehr_loinc_counts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

# 3836
print(nrow(pds.with.loinc.parts))
#  OOPS went up ?

# loinc method can legitimately be NA, but if any other column is NA,
# either the PDS count should be below the threshold
# or the LOINC code should be in the bad/old list above

###   ###   ###

# make some manual decisions

# time choice:
#   Pt TIME MUCH more common than next time (24H). XXX and * more common than 24H.
# PT or 24H time
# LP6960-1
# TIME
# Pt
# Point in time (spot)
#
# LP6924-7
# TIME
# 24H
# 24 hours

# scale choice:
#   use Qn SCALE only (5x more entries than next most commons, Nom and Ord)
# LP7753-9
# SCALE
# Qn
# Quantitative

###   ###   ###

# start merging in reviewed time, system, scale and property choices
# see csv reads at top of script
# put these in config?
pds.with.loinc.parts <-
  pds.with.loinc.parts[pds.with.loinc.parts$PROPERTY %in% turbo_loinc_properties_reviewed$PartTypeVal &
                         pds.with.loinc.parts$SCALE == "LP7753-9" &
                         pds.with.loinc.parts$SYSTEM %in% turbo_loinc_systems_reviewed$PartTypeVal &
                         (pds.with.loinc.parts$TIME == "LP6960-1" |
                            pds.with.loinc.parts$TIME == "LP6924-7") , ]

# ~ 1800
print(nrow(pds.with.loinc.parts))

###   ###   ###

has.numerator <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'NUMERATOR']

has.adjustment <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'ADJUSTMENT']

# exclude assays with numerators or adjustments
pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% intersect(has.numerator, has.adjustment),]
print(nrow(pds.with.loinc.parts))

###   ###   ###

has.challenge <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'CHALLENGE']

has.challenge <-
  intersect(has.challenge, pds.with.loinc.parts$LOINC)

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% has.challenge &
                                    LoincPartLink$PartTypeName == 'CHALLENGE']))
# include peak and trough?
acceptable.challenges <- c('peak', 'trough')
acceptable.challenges <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'CHALLENGE' &
                              LoincPartLink$PartName %in% acceptable.challenges]

unacceptable.challenges <-
  setdiff(has.challenge, acceptable.challenges)

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% unacceptable.challenges ,]
print(nrow(pds.with.loinc.parts))

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC &
                                    LoincPartLink$PartTypeName == 'CHALLENGE']))

###   ###   ###

# addressed some additional ways below
# assays that have a divisior part
has.divisor  <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'DIVISOR']

# assays that have a divisior part are are metnioned in PDS
has.divisor <- intersect(has.divisor, pds.with.loinc.parts$LOINC)

# table of divisor parts mentioned by PDS
sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% has.divisor &
                                    LoincPartLink$PartTypeName == 'DIVISOR']))

acceptable.divisors <- c('100 leukocytes', 'Creatinine')

# assays that mention an acceptable divisor
acceptable.divisors <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'DIVISOR' &
                              LoincPartLink$PartName %in% acceptable.divisors]

print(table(LoincPartLink$PartName[LoincPartLink$PartTypeName == 'DIVISOR' &
                                     LoincPartLink$LoincNumber %in% acceptable.divisors]))

unacceptable.divisors <-
  setdiff(has.divisor, acceptable.divisors)

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% unacceptable.divisors , ]

print(nrow(pds.with.loinc.parts))

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC &
                                    LoincPartLink$PartTypeName == 'DIVISOR']))

###   ###   ###

has.suffix  <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'SUFFIX']

has.suffix <- intersect(has.suffix, pds.with.loinc.parts$LOINC)

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% has.suffix &
                                    LoincPartLink$PartTypeName == 'SUFFIX']))

# include Ab.IgG and Ab.IgE?
acceptable.suffix <- c('Ab.IgG', 'Ab.IgE')
acceptable.suffix <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'SUFFIX' &
                              LoincPartLink$PartName %in% acceptable.suffix]

unacceptable.suffix <-
  setdiff(has.suffix, acceptable.suffix)

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% unacceptable.suffix , ]
print(nrow(pds.with.loinc.parts))

###   ###   ###

has.method  <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'METHOD']

has.method <- intersect(has.method, pds.with.loinc.parts$LOINC)

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% has.method &
                                    LoincPartLink$PartTypeName == 'METHOD']))

# include  Automated, Calculated, Manual, Electrophoresis, Automated count, Confirm, Manual count, IA ?
acceptable.method <-
  c(
    'Automated',
    'Calculated',
    'Manual',
    'Electrophoresis',
    'Automated count',
    'Confirm',
    'Manual count',
    'IA'
  )
acceptable.method <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'METHOD' &
                              LoincPartLink$PartName %in% acceptable.method]

unacceptable.method <-
  setdiff(has.method, acceptable.method)

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% unacceptable.method , ]
print(nrow(pds.with.loinc.parts))

### check merge lossses below
sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC &
                                    LoincPartLink$PartTypeName == 'METHOD']))

# all methods still there

###   ###   ###

multi.details <- c()

challenge.frame <-
  unique(LoincPartLink[LoincPartLink$PartTypeName == 'CHALLENGE'  &
                         LoincPartLink$LoincNumber %in%  pds.with.loinc.parts$LOINC , c('LoincNumber', "PartName")])
term.table <- table(challenge.frame$LoincNumber)

term.table <-
  cbind.data.frame(names(term.table), as.numeric(term.table))
names(term.table) <- c('LoincNumber', 'count')
term.table$LoincNumber <- as.character(term.table$LoincNumber)
multi.details <-
  union(multi.details, term.table$LoincNumber[term.table$count > 1])

###   ###   ###

divisor.frame <-
  unique(LoincPartLink[LoincPartLink$PartTypeName == 'DIVISOR'  &
                         LoincPartLink$LoincNumber %in%  pds.with.loinc.parts$LOINC , c('LoincNumber', "PartName")])
term.table <- table(divisor.frame$LoincNumber)

term.table <-
  cbind.data.frame(names(term.table), as.numeric(term.table))
names(term.table) <- c('LoincNumber', 'count')
term.table$LoincNumber <- as.character(term.table$LoincNumber)
multi.details <-
  union(multi.details, term.table$LoincNumber[term.table$count > 1])

###   ###   ###

suffix.frame <-
  unique(LoincPartLink[LoincPartLink$PartTypeName == 'SUFFIX'  &
                         LoincPartLink$LoincNumber %in%  pds.with.loinc.parts$LOINC , c('LoincNumber', "PartName")])
term.table <- table(suffix.frame$LoincNumber)

term.table <-
  cbind.data.frame(names(term.table), as.numeric(term.table))
names(term.table) <- c('LoincNumber', 'count')
term.table$LoincNumber <- as.character(term.table$LoincNumber)
multi.details <-
  union(multi.details, term.table$LoincNumber[term.table$count > 1])

###   ###   ###

method.frame <-
  unique(LoincPartLink[LoincPartLink$PartTypeName == 'METHOD'  &
                         LoincPartLink$LoincNumber %in%  pds.with.loinc.parts$LOINC , c('LoincNumber', "PartName")])

method.agg <-
  aggregate(
    method.frame$PartName,
    by = list(method.frame$LoincNumber),
    FUN = function(current.chunks) {
      return(paste0(sort(unique(current.chunks)), collapse = "|"))
    }
  )
names(method.agg) <- c('LoincNumber', 'agg')

term.table <- table(method.frame$LoincNumber)

term.table <-
  cbind.data.frame(names(term.table), as.numeric(term.table))
names(term.table) <- c('LoincNumber', 'count')
term.table$LoincNumber <- as.character(term.table$LoincNumber)

method.frame <- base::merge(x = method.agg, y = term.table)

single.collapsed.methods <-
  method.frame[method.frame$count == 1 |
                 method.frame$agg %in% c('Automated|Automated count', 'Manual|Manual count') ,]

ditchers <-
  setdiff(method.frame$LoincNumber,
          single.collapsed.methods$LoincNumber)

multi.details <-
  union(multi.details, ditchers)

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% multi.details ,]

###   ###   ###

# start merging for robot

cumulative.merges <- pds.with.loinc.parts

cumulative.merges <-
  base::merge(
    x = cumulative.merges,
    y = challenge.frame,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE,
    suffixes = c('.method', '.challenge')
  )

cumulative.merges <-
  base::merge(
    x = cumulative.merges,
    y = divisor.frame,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE,
    suffixes = c('.challenge', '.divisor')
  )

cumulative.merges <-
  base::merge(
    x = cumulative.merges,
    y = suffix.frame,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE,
    suffixes = c('.divisor', '.suffix')
  )

###   ###   ###

sort(table(LoincPartLink$PartName[LoincPartLink$LoincNumber %in% cumulative.merges$LOINC &
                                    LoincPartLink$PartTypeName == 'METHOD']))

cumulative.merges <-
  base::merge(
    x = cumulative.merges,
    y = single.collapsed.methods,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE,
    suffixes = c('', '.method')
  )

print(nrow(cumulative.merges))

any.dupes.Q <- table(cumulative.merges$LOINC)
any.dupes.Q <-
  cbind.data.frame(names(any.dupes.Q), as.numeric(any.dupes.Q))
colnames(any.dupes.Q) <- c('LOINC', 'count')
any.dupes.Q$LOINC <- as.character(any.dupes.Q$LOINC)

print(any.dupes.Q$LOINC[any.dupes.Q$count > 1])

pds.with.loinc.parts <- cumulative.merges


###   ###   ###

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

###   ###   ###

# components, esp. mapping

# ~ 1350

original.needs.component.mapping <-
  unique(pds.with.loinc.parts$COMPONENT)

current.needs.component.mapping <- original.needs.component.mapping
current.component.mapping.complete <- c()
current.component.mapping.frame <- matrix(ncol = 4)
colnames(current.component.mapping.frame) <-
  c('source.id', 'source.term', 'target.term', 'authority')

###   ###   ###

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

###   ###   ###

#### begin BioPortal mappings retrieval
# ~ 5 minutes vs remote BioPortal

retreived.and.parsed <-
  bp.map.retreive.and.parse(sort(unique(current.needs.component.mapping)))

retreived.and.parsed <-
  do.call(rbind.data.frame, retreived.and.parsed)
retreived.and.parsed[] <-
  lapply(retreived.and.parsed[], as.character)

# who did we get mapped to?
rap.tab <- table(retreived.and.parsed$target.ontology)
rap.tab <- cbind.data.frame(names(rap.tab), as.numeric(rap.tab))
names(rap.tab) <- c("target.ontology", "map.count")

# acceptable.ontologies <- config$acceptable.map.onts.from.loinc

harmonized_ontology_rankings$bp.onto.iri <-
  paste0('http://data.bioontology.org/ontologies/',
         harmonized_ontology_rankings$rhs)

# acceptable.retreived.and.parsed <-
#   unique(retreived.and.parsed[retreived.and.parsed$target.ontology %in%
#                                 acceptable.ontologies , c("source.term", "target.term", "target.ontology")])


acceptable.retreived.and.parsed <-
  base::merge(x = retreived.and.parsed ,
              y = harmonized_ontology_rankings ,
              by.x = 'target.ontology',
              by.y = 'bp.onto.iri')

agg <-
  aggregate(
    acceptable.retreived.and.parsed$priority,
    by = list(acceptable.retreived.and.parsed$source.term),
    FUN = min
  )
colnames(agg) <- c('source.term', 'priority')

acceptable.retreived.and.parsed <-
  base::merge(x = acceptable.retreived.and.parsed , y = agg)

acceptable.retreived.and.parsed <-
  unique(acceptable.retreived.and.parsed[, c("source.term", "target.term")])

acceptable.retreived.and.parsed$target.term <-
  as.character(acceptable.retreived.and.parsed$target.term)
acceptable.retreived.and.parsed$source.id <-
  gsub(pattern = "http://purl.bioontology.org/ontology/LNC/",
       replacement = "",
       x = acceptable.retreived.and.parsed$source.term)
acceptable.retreived.and.parsed$source.term <-
  as.character(acceptable.retreived.and.parsed$source.term)

update.accounting(
  acceptable.retreived.and.parsed,
  "source.id",
  "target.term",
  "BioPortal IRI mappings"
)

###   ###   ###

ptv.to.iri <-
  unique(turbo_loinc_systems_reviewed[, c('PartTypeVal', 'iri')])

pds.with.loinc.parts <-
  base::merge(
    x = pds.with.loinc.parts,
    y = ptv.to.iri,
    by.x  = "SYSTEM",
    by.y = "PartTypeVal",
    all.x  = TRUE
  )

# > table(pds.with.loinc.parts$TIME)
#
# LP6924-7 LP6960-1
# 53     1392

pds.with.loinc.parts$time.iri <- NA
pds.with.loinc.parts$time.iri[pds.with.loinc.parts$TIME == 'LP6924-7'] <-
  'obo:TURBO_0010724'

###   ###   ###

# starting string searches now
# use core analytes and searchable components (separately)
# with and without trailing s

# strip match labels of 'atom$' ?
# turns out not to be necessary when comparing query to synonym

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

current.needs.component.mapping.details <-
  LoincPartLink[LoincPartLink$PartTypeName == 'COMPONENT' &
                  LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC , ]

# > table( current.needs.component.mapping.details$Property , current.needs.component.mapping.details$LinkTypeName)
#
# DetailedModel Primary Search SyntaxEnhancement
# http://loinc.org/property/analyte                656       0      0                 0
# http://loinc.org/property/analyte-core             0       0      0               289
# http://loinc.org/property/COMPONENT                0     675      0                 0
# http://loinc.org/property/search                   0       0    116                 0

component.component.counts <-
  reshape2::dcast(
    data = current.needs.component.mapping.details,
    formula = LoincNumber ~ Property,
    fun.aggregate = length,
    value.var = 'Property'
  )
colnames(component.component.counts) <-
  sub(
    pattern = 'http://loinc.org/property/',
    replacement = '',
    x = colnames(component.component.counts)
  )
LoincNumber <- component.component.counts$LoincNumber
component.component.counts <-
  component.component.counts[, c("analyte", "analyte-core", "search")]
colnames(component.component.counts) <-
  paste0(colnames(component.component.counts), ".counts")

component.component.vals <-
  reshape2::dcast(
    data = current.needs.component.mapping.details,
    formula = LoincNumber ~ Property,
    value.var = 'PartName',
    fun.aggregate = function(current.chunks) {
      temp <- sort(unique(current.chunks))
      temp <- paste0(temp, collapse = '|')
      return(temp)
    }
  )
colnames(component.component.vals) <-
  sub(
    pattern = 'http://loinc.org/property/',
    replacement = '',
    x = colnames(component.component.vals)
  )
component.component.vals <-
  component.component.vals[, c("analyte", "analyte-core", "COMPONENT", "search"), ]
colnames(component.component.vals) <-
  paste0(colnames(component.component.vals), ".vals")

part.details <-
  cbind.data.frame(LoincNumber,
                   component.component.counts,
                   component.component.vals)

colnames(part.details) <- make.names(colnames(part.details))

part.details$analyte.vals[part.details$analyte.counts == 0] <- NA
part.details$analyte.core.vals[part.details$analyte.core.counts == 0] <-
  NA
part.details$search.vals[part.details$search.counts == 0] <-
  NA

part.details <-
  part.details[, setdiff(colnames(part.details),
                         c('analyte.counts', 'analyte.core.counts')),]

###   ###   ###

colnames(pds.with.loinc.parts) <-
  c(
    "SYSTEM",
    "LOINC",
    "LOINC_COUNT",
    "COMPONENT",
    "METHOD",
    "PROPERTY",
    "SCALE",
    "TIME",
    "PartName.challenge",
    "PartName.divisor",
    "PartName.suffix",
    "PartName.method",
    "method.count",
    "system.iri",
    "time.iri"
  )

colnames(part.details) <-
  c(
    "LoincNumber",
    "search.counts",
    "analyte.vals",
    "analyte.core.vals",
    "COMPONENT.vals",
    "search.vals"
  )

pds.with.loinc.parts <-
  base::merge(
    x = pds.with.loinc.parts ,
    y = part.details,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE
  )

current.component.mapping.frame <-
  current.component.mapping.frame[complete.cases(current.component.mapping.frame),]

num.to.name <- unique(Part[, c('PartNumber', 'PartName')])

current.component.mapping.frame <-
  base::merge(
    x = current.component.mapping.frame,
    y = num.to.name,
    by.x = 'source.id',
    by.y = 'PartNumber',
    all.x = TRUE
  )

num.to.name <- table(current.component.mapping.frame$source.id)
num.to.name <-
  cbind.data.frame(names(num.to.name), as.numeric(num.to.name))
names(num.to.name) <- c('PartNumber', 'map.count')
num.to.name$PartNumber <- as.character(num.to.name$PartNumber)

singles <- num.to.name$PartNumber[num.to.name$map.count == 1]
singles <-
  current.component.mapping.frame[current.component.mapping.frame$source.id %in% singles ,]

multi.mappers <- num.to.name$PartNumber[num.to.name$map.count > 1]

# multi.mappers <-
#   current.component.mapping.frame[current.component.mapping.frame$source.id %in% multi.mappers ,]
#
# multi.mappers$ontology <- multi.mappers$target.term
#
# multi.mappers$ontology <-
#   sub(pattern = '^http://purl.bioontology.org/ontology/',
#       replacement = '',
#       x = multi.mappers$ontology)
#
# multi.mappers$ontology <-
#   sub(pattern = '^http://purl.obolibrary.org/obo/',
#       replacement = '',
#       x = multi.mappers$ontology)
#
# multi.mappers$ontology <-
#   sub(pattern = '^http://purl.org/sig/ont/fma/fma',
#       replacement = 'fma_',
#       x = multi.mappers$ontology)
#
# multi.mappers$ontology <-
#   sub(pattern = '/.*$',
#       replacement = '',
#       x = multi.mappers$ontology)
#
# multi.mappers$ontology <-
#   sub(pattern = '_.*$',
#       replacement = '',
#       x = multi.mappers$ontology)
#
# # mm.onto.tab <- table(multi.mappers$ontology)
# # mm.onto.tab <-
# #   cbind.data.frame(names(mm.onto.tab), as.numeric(mm.onto.tab))
# # names(mm.onto.tab) <- c('ontology', 'count')
#
# #### STOP
#
# # this has already been taken care of above
# # if needed again, can use harmonized_ontology_rankings
#
# # ontology_rankings <- read_csv("ontology_rankings.csv")
# # ontology_rankings <- ontology_rankings[, c('ontology', 'priority')]
#
# multi.mappers <-
#   base::merge(
#     x = multi.mappers,
#     y = ontology_rankings,
#     by.x = 'ontology',
#     by.y = 'ontology' ,
#     all.x = TRUE
#   )
#
# mm.agg <-
#   aggregate(multi.mappers$priority, by = list(multi.mappers$source.id), min)
# colnames(mm.agg) <- c('source.id', 'priority')
#
# multi.mappers <-
#   unique(base::merge(x = multi.mappers,
#                      y = mm.agg))
#
#
# mmsi.table <- table(multi.mappers$source.id)
# mmsi.table <-
#   cbind.data.frame(names(mmsi.table), as.numeric(mmsi.table))
# names(mmsi.table) <- c('PartNumber', 'map.count')
# mmsi.table$PartNumber <- as.character(mmsi.table$PartNumber)
#
#
# # constrain here (or earlier on) to just human proteins?
# lost.cause <- mmsi.table$PartNumber[mmsi.table$map.count > 1]
# lost.cause <-
#   multi.mappers[multi.mappers$source.id %in% lost.cause ,]
#
# multi.mappers <-
#   multi.mappers[!(multi.mappers$source.id %in% lost.cause) , ]
# 
# multi.mappers <-
#   multi.mappers[, intersect(colnames(multi.mappers), colnames(singles))]

# singles <-
#   singles[, intersect(colnames(multi.mappers), colnames(singles))]

current.component.mapping.frame <-
  rbind.data.frame(singles, multi.mappers)

ccmf.name.term <-
  unique(current.component.mapping.frame[, c("PartName", "target.term")])

pds.with.loinc.parts <-
  base::merge(
    x = pds.with.loinc.parts,
    y = ccmf.name.term,
    by.x = 'analyte.core.vals',
    by.y = 'PartName' ,
    all.x = TRUE
  )

text.search.input <-
  unique(pds.with.loinc.parts[is.na(pds.with.loinc.parts$target.term), c('COMPONENT', 'analyte.core.vals')])

colnames(text.search.input) <- c('PartNumber', 'original.name')

###   ###   ###

source.relevant.links <-
  LoincPartLink[LoincPartLink$LoincNumber %in% ehr_loinc_counts$LOINC , ]

source.relevant.links <-
  source.relevant.links[source.relevant.links$PartNumber %in% original.needs.component.mapping ,]

###   ###   ###

current.needs.component.mapping.details <-
  source.relevant.links[source.relevant.links$Property == 'http://loinc.org/property/analyte-core' ,]

#### had been including these... might want to go back to that pattern
# system names and display names
# non-core names, display names, search terms
# divisor names


###   ###   ###

# also loook for searchable components

# # COMPONENT
# # Search
# # http://loinc.org/property/search

###   ###   ###

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


# below: not explicitly filtering out carent, plus, period etc
# not including searchable compoent tokens

text.search.input$singularized <- sub(pattern = 's$',
                                      replacement = '',
                                      x = text.search.input$original.name)

original <- text.search.input[, c('PartNumber', 'original.name')]
singularized <- text.search.input[, c('PartNumber', 'singularized')]
colnames(singularized) <- colnames(original)
text.search.input <-
  unique(rbind.data.frame(original, singularized))
colnames(text.search.input) <- c('PartNumber', 'any.name')

text.search.input <-
  text.search.input[order(text.search.input$any.name),]

###  DRY
table(text.search.input$PartNumber %in% current.component.mapping.complete)
table(text.search.input$any.name %in% Part$PartName[Part$PartNumber %in% current.component.mapping.complete])

# substitute '.' with ' '  for neutrophils.immature etc.
# substitute 'spp' with ''  for genus-level NCBI taxon entities


# search deeper than 99 (bicarbonate)
# or finally, use ontology filter

onto.list.for.ols <-
  paste0(sort(tolower(harmonized_ontology_rankings$rhs)), collapse = ',')

pre.ols.attempts <-
  apply(
    X = text.search.input,
    MARGIN = 1,
    FUN = function(current.row) {
      temp <-
        ols.serch.term.labels.universal(
          current.string = current.row[['any.name']],
          current.id = current.row[['PartNumber']],
          strip.final.s = FALSE,
          ontology.filter = paste0("&ontology=", onto.list.for.ols),
          kept.row.count = 600,
          req.exact = 'false'
        )
      if (is.data.frame(temp)) {
        if (nrow(temp) > 0) {
          # other options?
          # colum-mismatch tolerating rbind?
          temp <- temp[, setdiff(colnames(temp), "description")]
          return(temp)
        }
      }
    }
  )

col.counts <- sapply(pre.ols.attempts, length)
table(col.counts)

niners <- pre.ols.attempts[col.counts == 9]
niners <- do.call(rbind.data.frame, niners)
niners$synonym <- ""

tenners <- pre.ols.attempts[col.counts == 10]
tenners <- do.call(rbind.data.frame, tenners)

ols.attempts <- rbind.data.frame(tenners, niners)

# ols.attempts.saved <- ols.attempts

# bicarb LP15441-6
# I submitted some issues to the OLS github

# do ontology filtering HERE?
# don't bother calculating sting distnace is we're not going to use the term?

# search.res.dists <-
#   apply(
#     X = ols.attempts,
#     MARGIN = 1,
#     FUN = function(current.row) {
#       print(current.row$query)
#       temp <- unique(c(current.row$label, current.row$synonym))
#       inner <- lapply(
#         X = temp,
#         FUN = function(current.single) {
#           current.single <- as.character(current.single)
#           the.dist <-
#             stringdist::stringdist(
#               a = tolower(current.row$query),
#               b = tolower(current.single),
#               method = 'cosine'
#             )
#           returnable <- list(
#             candidated.id = current.row$short_form,
#             candidate.str = current.single,
#             cosine = the.dist
#           )
#           return(returnable)
#         }
#       )
#       inner <- do.call(rbind.data.frame, inner)
#       min.row <- which.min(inner$cosine)
#       best.row <- inner[min.row, ]
#       best.row$query <- current.row$query
#       best.row$part.number <- current.row$loinc.part
#       best.row$ont.name <- current.row$ontology_name
#       return(best.row)
#     }
#   )

# ols.attempts <- ols.attempts[ols.attempts$loinc.part == 'LP15441-6' , ]

# ols.attempts <- ols.attempts.saved

search.res.dists <-
  apply(
    X = ols.attempts,
    MARGIN = 1,
    FUN = function(current.row) {
      print(current.row$query)
      temp <- list(current.row$label, current.row$synonym)
      # print(temp)
      temp <- lapply(temp, tolower)
      # print(temp)
      if (tolower(current.row$query) %in% temp) {
        return(current.row[c(
          "iri",
          "short_form",
          "obo_id",
          "label",
          "ontology_name",
          "ontology_prefix",
          "query",
          "loinc.part",
          "rank"
        )])
      }
    }
  )

col.counts <- sapply(search.res.dists, length)
table(col.counts)

niners <- search.res.dists[col.counts == 9]

col.counts <- sapply(niners, length)
table(col.counts)

search.res.dists <- do.call(rbind.data.frame, niners)
# search.res.dists <- data.table::rbindlist(niners)
# search.res.dists[] <- do.call(as.character, search.res.dists[])

placeholder <- lapply(c(
  "iri",
  "short_form",
  "obo_id",
  "label",
  "ontology_name",
  "ontology_prefix",
  "query",
  "loinc.part"
), function(current.col) {
  search.res.dists[, current.col] <<-
    as.character(search.res.dists[, current.col])
  return()
})

kept.best <- base::merge(x = search.res.dists,
                         y = harmonized_ontology_rankings,
                         by.x = 'ontology_prefix',
                         by.y = 'rhs')

# need to merge priorities in, not ranks

agg <-
  aggregate(
    kept.best$priority,
    by = list(kept.best$loinc.part),
    FUN = min
  )
names(agg) <- c('loinc.part', 'priority')

#### STOP

kept.best <- base::merge(
  x = kept.best,
  y = agg,
  by.x = c('loinc.part', 'priority'),
  by.y = c('loinc.part', 'priority')
)


# # search.res.dists <- data.table::rbindlist(tenners)
# 
# # # super slow
# # search.res.dists <- do.call(rbind.data.frame, search.res.dists)
# 
# search.res.dists$candidated.id <-
#   as.character(search.res.dists$candidated.id)
# search.res.dists$candidate.str <-
#   as.character(search.res.dists$candidate.str)
# 
# ###   ###   ###
# 
# # check for acceptability (cosine and source ontology)
# 
# # there's probably lots of good ones with slightly decreased consine distances,
# # but lets' start conservatively
# 
# # trailing 'atom' had been a problem in ChEBI
# # not any more thanks to applying cosine distance to preffered labels and synonyms
# 
# # broke away from using update mechanism?
# 
# # also prioritize hits by ontology
# 
# ###   ###   ###
# 
# # put overall best AFTER ontology filtering
# # prob not because we're only taking cosine == 0 anyway
# 
# overall.best <-
#   aggregate(search.res.dists$cosine,
#             by = list(search.res.dists$part.number),
#             min)
# 
# colnames(overall.best) <- c('part.number', 'cosine')
# 
# overall.best <- base::merge(search.res.dists, overall.best)
# 
# overall.best <-
#   overall.best[, c("part.number",
#                    "query",
#                    "cosine",
#                    "candidate.str",
#                    "candidated.id",
#                    "ont.name")]
# 
# hist(overall.best$cosine, breaks = 99)
# 
# ob.on.tab <- sort(table(overall.best$ont.name))
# ob.on.tab <-
#   cbind.data.frame(names(ob.on.tab), as.numeric(ob.on.tab))
# names(ob.on.tab) <- c('ontology', 'count')
# 
# # > sort(table(overall.best$ont.name))
# #
# # aeo      agro      bspo      caro      cmpo    co_323    co_324    co_325    co_331    co_338    co_341    co_343    co_347    co_365     cteno
# # 1         1         1         1         1         1         1         1         1         1         1         1         1         1         1
# # ddanat       duo   ecocore       fao       ido       lbo       mco       mop       ogi      ogms       omp       opl      peco     planp      poro
# # 1         1         1         1         1         1         1         1         1         1         1         1         1         1         1
# # pride        pw       rex      symp      tads   taxrank        to    co_320    co_322    co_330    co_333    co_339    co_340    co_350    co_357
# # 1         1         1         1         1         1         1         2         2         2         2         2         2         2         2
# # co_360       fix     hpath     micro      ncro       oae     oarcs      obib      srao     stato        vt      zeco    co_337    co_356    co_358
# # 2         2         2         2         2         2         2         2         2         2         2         2         3         3         3
# # co_366       enm     flopo      miro       mmo       mod        po        zp    co_321    co_335        eo     mondo        ms    ehdaa2        mi
# # 3         3         3         3         3         3         3         3         4         4         4         4         4         5         5
# # pato       spd    unimod       hao     nomen      fobi       oba        rs      envo      fypo      sibo       clo       sio       tto      doid
# # 5         5         5         6         6         7         7         7         8         8         9        10        10        10        11
# # vto       aro    idomal        uo        ma       obi       xao       ero     plana      tgma      chmo        vo     emapa        mp       ogg
# # 11        12        12        12        13        13        14        15        15        15        16        17        18        19        19
# # edam        hp       afo        cl       cmo    co_334     nmrcv        so      scdo    uberon       bto       bao       zfa      ordo       gaz
# # 22        22        26        27        28        29        29        29        36        37        44        52        52        59        60
# # foodon       xco       efo        go       fma      dron        pr ncbitaxon     chebi      omit      ncit
# # 62        70        71       106       115       180       202       422       546       680      1340
# 
# ####    ####    ####    ####
# 
# # get keepers from XXX
# 
# keepers <- harmonized_ontology_rankings$rhs
# 
# kept.best <-
#   overall.best[tolower(overall.best$ont.name) %in% tolower(keepers) ,]
# 
# # # create prioritized_ols_ontology_names.csv
# # kb.tab <- sort(table(overall.best$ont.name))
# # kb.tab <- cbind.data.frame(names(kb.tab), as.numeric(kb.tab))
# # names(kb.tab) <- c('ontology', 'count')
# 
# # prioritized_ols_ontology_names <-
# #   read_csv("prioritized_ols_ontology_names.csv")
# # prioritized_ols_ontology_names <-
# #   prioritized_ols_ontology_names[, c('ontology', 'priority')]
# 
# kept.best$ont.name <- tolower(kept.best$ont.name)
# 
# harmonized_ontology_rankings$lc.rhs <-
#   tolower(harmonized_ontology_rankings$rhs)
# 
# # kept.best <-
# #   base::merge(
# #     x = kept.best,
# #     y = prioritized_ols_ontology_names,
# #     by.x = 'ont.name',
# #     by.y = 'ontology' ,
# #     all.x = TRUE
# #   )
# 
# kept.best <-
#   base::merge(
#     x = kept.best,
#     y = harmonized_ontology_rankings,
#     by.x = 'ont.name',
#     by.y = 'lc.rhs' ,
#     all.x = TRUE
#   )
# 
# kb.agg <-
#   aggregate(kept.best$priority, by = list(kept.best$part.number), min)
# colnames(kb.agg) <- c('part.number', 'priority')
# 
# kept.best <-
#   unique(base::merge(
#     x = kept.best,
#     y = kb.agg,
#     by.x = c('part.number', 'priority'),
#     by.y =  c('part.number', 'priority'),
#     all.x = FALSE
#   ))
# 
# kept.best <- kept.best[kept.best$cosine == 0 , ]

###   ###   ###

# deal with these

# excluded becasue no one definitive component mapping
# like lost causes above?
overmapped <- table(pds.with.loinc.parts$LOINC)
overmapped <-
  cbind.data.frame(names(overmapped), as.numeric(overmapped))
colnames(overmapped) <- c('LOINC', 'count')
overmapped <- overmapped$LOINC[overmapped$count > 1]
overmapped <-
  pds.with.loinc.parts[pds.with.loinc.parts$LOINC %in% overmapped,]

pds.with.loinc.parts <-
  pds.with.loinc.parts[!pds.with.loinc.parts$LOINC %in% overmapped$LOINC,]

# backmerge <-
#   unique(kept.best[, c('part.number', 'candidated.id'), ])
# backmerge$iri.from.search <-
#   paste0('obo:', backmerge$candidated.id)
# backmerge <- backmerge[, c('part.number', 'iri.from.search')]

backmerge <-
  unique(kept.best[, c('loinc.part', 'short_form'), ])
backmerge$iri.from.search <-
  paste0('obo:', backmerge$short_form)
backmerge <- backmerge[, c('loinc.part', 'iri.from.search')]


# pds.with.loinc.parts <-
#   base::merge(
#     x = pds.with.loinc.parts,
#     y = backmerge,
#     by.x = 'COMPONENT',
#     by.y = 'part.number',
#     all.x = TRUE
#   )

pds.with.loinc.parts <-
  base::merge(
    x = pds.with.loinc.parts,
    y = backmerge,
    by.x = 'COMPONENT',
    by.y = 'loinc.part',
    all.x = TRUE
  )

pds.with.loinc.parts$final.comp <- pds.with.loinc.parts$target.term
pds.with.loinc.parts$final.comp[is.na(pds.with.loinc.parts$final.comp)] <-
  pds.with.loinc.parts$iri.from.search[is.na(pds.with.loinc.parts$final.comp)]

####

ready.for.robot <-
  pds.with.loinc.parts[!is.na(pds.with.loinc.parts$final.comp), ]

# excluded because NO automatically-accepted mapping for component
drawing.board <-
  pds.with.loinc.parts[is.na(pds.with.loinc.parts$final.comp), ]

robot.row.count <- nrow(ready.for.robot)
all.blanks <- rep(x = '', robot.row.count)

term.id <- paste0('lobo:LOBO_', 1:robot.row.count)
term.label <- paste0('LOBO term ', 1:robot.row.count)
has.curation.status <- all.blanks

# alternative term ... LOINC assay label
at.frame <-
  unique(LoincPartLink[, c('LoincNumber', 'LongCommonName')])
at.index <- match(ready.for.robot$LOINC, at.frame$LoincNumber)
alternative.term <- at.frame$LongCommonName[at.index]

definition <- all.blanks

# definition.source <-
#   rep(x = 'Mark Andrew Miller, ORCID:0000-0001-9076-6066|Chris Stoeckert', robot.row.count)

# maybe a link to the script that populated the template?
definition.source <- all.blanks

example.of.usage <- all.blanks
editor.note <- all.blanks
term.editor <- definition.source
ontology.term.requester <- all.blanks
term.tracker.item <-
  rep(x = 'https://github.com/obi-ontology/obi/issues/1153', robot.row.count)

# subclass or equivalent
logical.type <-
  rep(x = 'subclass', robot.row.count)

logical.note <- all.blanks

#	clinical assay from TURBO. analyte could be inferred. annotate this with ...
#   see recent notes from CJS.
#   Also quantitative by default and PT by default
parent.class <- rep(x = "'clinical assay'", robot.row.count)

# 24 collection axiom as above if neccessary (LP6924-7 TURBO_0010724)
# switch and create 24 specimen collection objective specification
material.processing.technique <- all.blanks
material.processing.technique[!is.na(ready.for.robot$time.iri) & ready.for.robot$time.iri == 'obo:TURBO_0010724'] <-
  "'specimen collection process' and ( achieves_planned_objective some '24 hour sample collection objective' )"

# skip for now. would be based on method
detection.technique <- all.blanks

# specimens with | split. union or intersection? subclass or eq class?
evaluant <- ready.for.robot$system.iri

# put all components in here and wait for errors from things that aren't molecular entities? V
analyte <- ready.for.robot$final.comp
device <- all.blanks
reagent <- all.blanks
molecular.label <- all.blanks
# skip in favor of analyte/target entity for now	some specific kind of datum based on LOINC property.
input <- all.blanks

# datums may need to be created and defined in terms of potentiatlly new units.	 V
output <- rep(x = 'obo:IAO_0000032', robot.row.count)

# move non-molecular entity analytes into here V
target.entity <- all.blanks

# will probably be using this soon
objective <- all.blanks

associated.axioms <- all.blanks
associated.axioms[!is.na(ready.for.robot$PartName.divisor) &
                    ready.for.robot$PartName.divisor == 'Creatinine'] <- "'has part' some 'division by creatinine concentration normalization'"
associated.axioms[!is.na(ready.for.robot$PartName.divisor) &
                    ready.for.robot$PartName.divisor == '100 leukocytes'] <-  "'has part' some 'division by 100 leukocytes normalization'"

# or LOINC  https://loinc.org/2339-0/ style ?
database.cross.reference.IRI <-
  paste0('http://purl.bioontology.org/ontology/LNC/',
         ready.for.robot$LOINC)


ready.for.robot <- cbind.data.frame(
  term.id,
  term.label,
  has.curation.status,
  alternative.term,
  definition,
  definition.source,
  example.of.usage,
  editor.note,
  term.editor,
  ontology.term.requester,
  term.tracker.item,
  logical.type,
  logical.note,
  parent.class,
  material.processing.technique,
  detection.technique,
  evaluant,
  analyte,
  device,
  reagent,
  molecular.label,
  input,
  output,
  target.entity,
  objective,
  associated.axioms,
  database.cross.reference.IRI
)

write.table(ready.for.robot, file = 'lobo_template_headerless.csv', row.names = FALSE, col.names = FALSE, sep = ',')

write.table(
  lobo_template_headers_only,
  file = 'lobo_template_headers_only.csv',
  row.names = FALSE,
  col.names = FALSE,
  sep = ','
)

###   ###   ###

# reworked sting inputs... no longer contains anything but unmapped core analyte components

bulk.systems <-
  kept.best[kept.best$ont.name %in% c('cl', 'uberon') &
              kept.best$part.number %in% Part$PartNumber[Part$PartTypeName == 'SYSTEM'] &
              kept.best$cosine == 0 ,]

# see turbo_loinc_systems_reviewed above

###   ###   ###

# some challenges

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

## some X+Y terms
# vitamin d2
# ChEBI
# CHEBI
# class
# TRUE
# A vitamin D supplement and has been isolated from alfalfa.
# calciferol
#
# isovaleryl-l-carnitine
# ChEBI
# CHEBI
# class
# TRUE
# An O-isovalerylcarnitine that is the 3-methylbutanoyl (isovaleryl) derivative of L-carnitine.
# isovalerylcarnitine
#
# desmethyldoxepin
# ChEBI
# CHEBI
# class
# TRUE
# A dibenzooxepine resulting from the demethylation of the antidepressant doxepin. It is the active metabolite of doxepin.
# nordoxepin
#
# 2-methylbutyrylcarnitine
# ChEBI
# CHEBI
# class
# TRUE
# A C5-acylcarnitine having 2-methylbutyryl as the acyl substituent.
# methylbutyrylcarnitine (c5)
#
# --- doesn't look promising below ---
#
# norclomipramine
#
# methsuximide
# normethsuximide

###   ###   ###

# Make sure template starts with both headers... see XXX shell script
#
# robot template --prefix "lobo: http://example.com/lobo/" -t lobo_template.csv
