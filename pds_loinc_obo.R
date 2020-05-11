source('turbo_loinc_obo_setup.R')

# true.challenges <- c('Glucose')

true.challenges <- c(
  '1H post 50 g glucose PO',
  '1H post 50 g lactose PO',
  '1H post dose corticotropin',
  '2H post 50 g lactose PO' ,
  '2H post 75 g glucose PO',
  '30M post 50 g lactose PO',
  '30M post dose corticotropin' ,
  '3H post 100 g glucose PO',
  '3H post 50 g lactose PO',
  '8H post 1 mg dexamethasone PO overnight'
)


pk.in.challenge.slot <- c('peak', 'trough')
divisors <- c('100 leukocytes', 'Creatinine')
serology.suffixes <- c("Ab",
                       "Ab.IgA",
                       "Ab.IgE",
                       "Ab.IgG",
                       "Ab.IgM")
loinc.methods <-
  c(# 'Automated',
    'Calculated',
    # 'Manual',
    'Electrophoresis',
    'Automated count',
    'Confirm',
    'Manual count',
    'IA')

#### PUT FILENAMES/PATHS IN CONFIG
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

# THIS IS HANDY FOR DIGNOSISNG/PRIORITIZING PATTERNS FOR MANUAL CURATION
# MultiAxialHierarchy <-
#   read_csv(config$MultiAxialHierarchy.file)

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

# this isn’t all LOINC codes used in the EHR, is it???
# it does go down to counts of 1!
# how was it generated?
print(nrow(ehr_loinc_counts))

common.loincs <-
  ehr_loinc_counts$LOINC[ehr_loinc_counts$LOINC_COUNT >= config$common.threshold]

# ~ 3440 r_lab_orders ordered for >= (2) unique EMPIs
print(length(common.loincs))

# get part values for the "common", LOINC-annotated ehr lab results
relevant.primary.part.cast <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs &
                  LoincPartLink$LinkTypeName == "Primary",]

relevant.primary.part.cast <-
  relevant.primary.part.cast[, c("LoincNumber", "PartNumber", "PartTypeName")]

relevant.primary.part.cast <-
  dcast(data = relevant.primary.part.cast,
        formula = LoincNumber ~ PartTypeName,
        value.var = "PartNumber")

#  ~ 3415 LOINC terms that overlap with the "common" terms present in the EHR
print(nrow(relevant.primary.part.cast))

# ~ 25 bad or old LOINCs in PDS?

# setdiff(common.loincs, relevant.primary.part.cast$LoincNumber)
# temp <-
#   ehr_loinc_counts[ehr_loinc_counts$LOINC %in% setdiff(common.loincs, relevant.primary.part.cast$LoincNumber),]
# R_LAB_RESULT_ITEM <-
#   read.csv("~/loinc_in_obo/R_LAB_RESULT_ITEM_202005050841.csv",
#            stringsAsFactors = FALSE)
# temp <- base::merge(x = temp, y = R_LAB_RESULT_ITEM)
# write.csv(temp, file =  "bad_pds_loinc_code_usage.csv")

# [1] "11575-5"  "41613"    "5876-7"   "103331-7" "5778-7"   "280602-2" "74217-8"  "39478-6"  "SOLOINC"  "18185-8"  "34695-0"
# [12] "6690-3"   "786-5"    "777-4"    "X29142-7" "UPE18"    "3148-9"   "4544-4"   "789-9"    "0051531"  "788-1"    "787-3"
# [23] "2890-3"   "L4985"    "718-8"    "785-7"    "L4613"

###   ###   ###

# first EHR relevance/LOINC content merge!


ehr.with.loinc.parts <-
  left_join(x = ehr_loinc_counts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

# 3836
print(nrow(ehr.with.loinc.parts))
# #  OOPS went up ?
##   actually got joined with all orders, even when only ordered for a single person
# could re-filter


# would we loose anything easy to implement?
# actually, even relevant.primary.part.cast has been stripped of uncommon loinc codes at this point
ehr.rare <-
  ehr.with.loinc.parts[!(ehr.with.loinc.parts$LOINC %in% common.loincs) , ]
ehr.with.loinc.parts <-
  ehr.with.loinc.parts[ehr.with.loinc.parts$LOINC %in% common.loincs , ]
# back TO 3440
print(nrow(ehr.with.loinc.parts))


####

# loinc method can legitimately be NA, but if any other column is NA,
# either the ehr count should be below the threshold
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
# put these ACCEPTED SCALE AND TIME values in config?
# we will ahve to get to additional scales and times in teh future
# may require new ontology terms or new mapping patterns

lost.due.to.time.scale <- ehr.with.loinc.parts$LOINC
ehr.with.loinc.parts <-
  ehr.with.loinc.parts[ehr.with.loinc.parts$PROPERTY %in% turbo_loinc_properties_reviewed$PartTypeVal &
                         ehr.with.loinc.parts$SCALE == "LP7753-9" &
                         ehr.with.loinc.parts$SYSTEM %in% turbo_loinc_systems_reviewed$PartTypeVal &
                         (ehr.with.loinc.parts$TIME == "LP6960-1" |
                            ehr.with.loinc.parts$TIME == "LP6924-7") , ]

lost.due.to.time.scale <-
  setdiff(lost.due.to.time.scale, ehr.with.loinc.parts$LOINC)

# ~ 1800
print(nrow(ehr.with.loinc.parts))

###   ###   ###

has.numerator <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'NUMERATOR']

has.adjustment <-
  LoincPartLink$LoincNumber[LoincPartLink$PartTypeName == 'ADJUSTMENT']

# exclude assays with numerators or adjustments
ehr.with.loinc.parts <-
  ehr.with.loinc.parts[!ehr.with.loinc.parts$LOINC %in% intersect(has.numerator, has.adjustment),]
print(nrow(ehr.with.loinc.parts))


challenge.split <-
  split.details('CHALLENGE', c(true.challenges, pk.in.challenge.slot))
divisor.split <-
  split.details('DIVISOR', divisors)
suffix.split <-
  split.details('SUFFIX', serology.suffixes)
method.split <- split.details('METHOD',
                              loinc.methods)

###   ###   ###

# start merging for robot

cumulative.merges <- ehr.with.loinc.parts

detail.structures <-
  list(challenge.split, divisor.split, method.split, suffix.split)

placeholder <-
  lapply(detail.structures, function(current.structure) {
    temp <- current.structure[[1]]
    print(head(temp))
    cumulative.merges <<-
      base::merge(
        x = cumulative.merges,
        y = temp,
        by.x = 'LOINC',
        by.y = 'LoincNumber',
        all.x = TRUE
      )
  })

print(nrow(cumulative.merges))

any.dupes.Q <- table(cumulative.merges$LOINC)
any.dupes.Q <-
  cbind.data.frame(names(any.dupes.Q), as.numeric(any.dupes.Q))
colnames(any.dupes.Q) <- c('LOINC', 'count')
any.dupes.Q$LOINC <- as.character(any.dupes.Q$LOINC)

print(any.dupes.Q$LOINC[any.dupes.Q$count > 1])

ehr.with.loinc.parts <- cumulative.merges

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

# ehr.with.loinc.parts: one row per assay/result type
# columns for the assay/result LOINC code and each LOINC parts code
# head(ehr.with.loinc.parts)
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
  unique(ehr.with.loinc.parts$COMPONENT)

current.needs.component.mapping <- original.needs.component.mapping
current.component.mapping.complete <- c()
current.component.mapping.frame <- matrix(ncol = 4)
colnames(current.component.mapping.frame) <-
  c('source.id', 'source.term', 'target.term', 'authority')

###   ###   ###

# LOINC asserts some mappings of their own

loinc.provided.component.mappings <-
  PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% current.needs.component.mapping , ]

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

print(loinc.provided.rxnorm.only)

# adalimumab (HUMIRA) http://purl.obolibrary.org/obo/DRON_00018971
# opiates [confirmed not present verbatim in DrOn or ChEBI]

# SO JUST USE ChEBI
loinc.provided.component.mappings <-
  loinc.provided.component.mappings[loinc.provided.component.mappings$ExtCodeSystem  == "https://www.ebi.ac.uk/chebi" &
                                      loinc.provided.component.mappings$Equivalence == "equivalent" , ]

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

bp.retreived.and.parsed <-
  bp.map.retreive.and.parse(sort(unique(current.needs.component.mapping)))

bp.retreived.and.parsed <-
  do.call(rbind.data.frame, bp.retreived.and.parsed)
bp.retreived.and.parsed[] <-
  lapply(bp.retreived.and.parsed[], as.character)

# who did we get mapped to?
rap.tab <- table(bp.retreived.and.parsed$target.ontology)
rap.tab <- cbind.data.frame(names(rap.tab), as.numeric(rap.tab))
names(rap.tab) <- c("target.ontology", "map.count")

# acceptable.ontologies <- config$acceptable.map.onts.from.loinc

harmonized_ontology_rankings$bp.onto.iri <-
  paste0('http://data.bioontology.org/ontologies/',
         harmonized_ontology_rankings$rhs)

acceptable.bp.retreived.and.parsed <-
  base::merge(x = bp.retreived.and.parsed ,
              y = harmonized_ontology_rankings ,
              by.x = 'target.ontology',
              by.y = 'bp.onto.iri')

lowest.priority <-
  aggregate(
    acceptable.bp.retreived.and.parsed$priority,
    by = list(acceptable.bp.retreived.and.parsed$source.term),
    FUN = min
  )
colnames(lowest.priority) <- c('source.term', 'priority')

acceptable.bp.retreived.and.parsed <-
  base::merge(x = acceptable.bp.retreived.and.parsed , y = lowest.priority)

acceptable.bp.retreived.and.parsed <-
  unique(acceptable.bp.retreived.and.parsed[, c("source.term", "target.term")])

acceptable.bp.retreived.and.parsed$target.term <-
  as.character(acceptable.bp.retreived.and.parsed$target.term)
acceptable.bp.retreived.and.parsed$source.id <-
  gsub(pattern = "http://purl.bioontology.org/ontology/LNC/",
       replacement = "",
       x = acceptable.bp.retreived.and.parsed$source.term)
acceptable.bp.retreived.and.parsed$source.term <-
  as.character(acceptable.bp.retreived.and.parsed$source.term)

update.accounting(
  acceptable.bp.retreived.and.parsed,
  "source.id",
  "target.term",
  "BioPortal IRI mappings"
)

###   ###   ###

ptv.to.iri <-
  unique(turbo_loinc_systems_reviewed[, c('PartTypeVal', 'iri')])
colnames(ptv.to.iri) <- c('PartTypeVal', 'system.iri')

ehr.with.loinc.parts <-
  base::merge(
    x = ehr.with.loinc.parts,
    y = ptv.to.iri,
    by.x  = "SYSTEM",
    by.y = "PartTypeVal",
    all.x  = TRUE
  )

# > table(ehr.with.loinc.parts$TIME)
#
# LP6924-7 LP6960-1
# 53     1392

#### 20200505 SKIP THIS? just the LOINC timing part number and work out with rules in ready for robot section?

# ehr.with.loinc.parts$time.iri <- NA
# ehr.with.loinc.parts$time.iri[ehr.with.loinc.parts$TIME == 'LP6924-7'] <-
#   'obo:TURBO_0010724'

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
                  LoincPartLink$LoincNumber %in% ehr.with.loinc.parts$LOINC ,]

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
  component.component.vals[, c("analyte", "analyte-core", "COMPONENT", "search"),]

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
                         c('analyte.counts', 'analyte.core.counts')), ]

###   ###   ###

ehr.with.loinc.parts <-
  base::merge(
    x = ehr.with.loinc.parts ,
    y = part.details,
    by.x = 'LOINC',
    by.y = 'LoincNumber',
    all.x = TRUE
  )

current.component.mapping.frame <-
  current.component.mapping.frame[complete.cases(current.component.mapping.frame), ]

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
  current.component.mapping.frame[current.component.mapping.frame$source.id %in% singles , ]

# # could theoretically require followup, prioritizaston, etc.
# multi.mappers <- num.to.name$PartNumber[num.to.name$map.count > 1]
#
#
# current.component.mapping.frame <-
#   rbind.data.frame(singles, multi.mappers)

current.component.mapping.frame <- singles

ccmf.name.term <-
  unique(current.component.mapping.frame[, c("PartName", "target.term")])

ehr.with.loinc.parts <-
  base::merge(
    x = ehr.with.loinc.parts,
    y = ccmf.name.term,
    by.x = 'analyte.core.vals',
    by.y = 'PartName' ,
    all.x = TRUE
  )

text.search.input <-
  unique(ehr.with.loinc.parts[is.na(ehr.with.loinc.parts$target.term), c('COMPONENT', 'analyte.core.vals')])

colnames(text.search.input) <- c('PartNumber', 'original.name')

source.relevant.links <-
  LoincPartLink[LoincPartLink$LoincNumber %in% ehr_loinc_counts$LOINC ,]

source.relevant.links <-
  source.relevant.links[source.relevant.links$PartNumber %in% original.needs.component.mapping , ]

current.needs.component.mapping.details <-
  source.relevant.links[source.relevant.links$Property == 'http://loinc.org/property/analyte-core' , ]

#### had been including these... might want to go back to that pattern
# system names and display names
# non-core names, display names, search terms
# divisor names

# also loook for searchable components

# # COMPONENT
# # Search
# # http://loinc.org/property/search

###   ###   ###

all.part.types.in.source <- ehr.with.loinc.parts$LOINC
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


# below: not explicitly filtering out caret, plus, period etc
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
  text.search.input[order(text.search.input$any.name), ]

###  DRY
table(text.search.input$PartNumber %in% current.component.mapping.complete)
table(text.search.input$any.name %in% Part$PartName[Part$PartNumber %in% current.component.mapping.complete])

# substitute '.' with ' ' for neutrophils.immature etc.
# substitute 'spp$' or 'sp$' with ''  for genus-level NCBI taxon entities

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
          # matching bicarbonate aginst ChEBI hydrogencarbonate has required going 500 to 550 deep
          #   because it's "just" an annotations_trimmed
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

tenners <- pre.ols.attempts[col.counts == 10]
tenners <- do.call(rbind.data.frame, tenners)
tenners$synonym <- ""

elevens <- pre.ols.attempts[col.counts == 11]
elevens <- do.call(rbind.data.frame, elevens)

ols.attempts <- rbind.data.frame(elevens, tenners)

# ols.attempts.saved <- ols.attempts

# don't bother calculating string distnace
#   if we're only going to use idnetical matches

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

ols.deep.check <-
  apply(
    X = ols.attempts,
    MARGIN = 1,
    FUN = function(current.row) {
      print(current.row$query)
      temp <-
        list(current.row$label,
             current.row$synonym,
             current.row$annotations_trimmed)
      # print(temp)
      temp <- sort(unique(unlist(lapply(temp, tolower))))
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

col.counts <- sapply(ols.deep.check, length)
table(col.counts)


niners <- ols.deep.check[col.counts == 9]

col.counts <- sapply(niners, length)
table(col.counts)

ols.deep.check <- do.call(rbind.data.frame, niners)

hist(ols.deep.check$rank, breaks = 99)
# ols.deep.check <- data.table::rbindlist(niners)

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
  ols.deep.check[, current.col] <<-
    as.character(ols.deep.check[, current.col])
  return()
})


kept.best <- base::merge(x = ols.deep.check,
                         y = harmonized_ontology_rankings,
                         by.x = 'ontology_prefix',
                         by.y = 'rhs')

# # excluded if no one definitive component mapping
# # like lost causes

lowest.priority <-
  aggregate(kept.best$priority,
            by = list(kept.best$loinc.part),
            FUN = min)

names(lowest.priority) <- c('loinc.part', 'priority')

kept.best <- base::merge(
  x = kept.best,
  y = lowest.priority,
  by.x = c('loinc.part', 'priority'),
  by.y = c('loinc.part', 'priority')
)

# oops, lost count and label match while tidying up
# coudl just retain the lowest rank

lowest.rank <-
  aggregate(kept.best$rank,
            by = list(kept.best$loinc.part),
            FUN = min)

names(lowest.rank) <- c('loinc.part', 'rank')

kept.best <- base::merge(
  x = kept.best,
  y = lowest.rank,
  by.x = c('loinc.part', 'rank'),
  by.y = c('loinc.part', 'rank')
)

# kept.best <-
#   kept.best[kept.best$count == 1 | kept.best$label.match , ]

# # ols.deep.check <- data.table::rbindlist(elevens)
#
# # # previous super slow attempt
# # ols.deep.check <- do.call(rbind.data.frame, ols.deep.check)

# # check for acceptability (cosine and source ontology)
#
# # there's probably lots of good ones with slightly decreased consine distances,
# # but lets' start conservatively
#
# # trailing 'atom' had been a problem in ChEBI
# # not any more thanks to applying cosine distance to preffered labels and synonyms

# hist(overall.best$cosine, breaks = 99)

backmerge <-
  unique(kept.best[, c('loinc.part', 'short_form'), ])
# assumes that the short forms ahve all been drawn from OBO ontologies
# that's not necessarily the case woth OLS or BioPortal
# but we have ensured that with our harmonized ontology rankings
backmerge$iri.from.search <-
  paste0('obo:', backmerge$short_form)
backmerge <- backmerge[, c('loinc.part', 'iri.from.search')]

ehr.with.loinc.parts <-
  base::merge(
    x = ehr.with.loinc.parts,
    y = backmerge,
    by.x = 'COMPONENT',
    by.y = 'loinc.part',
    all.x = TRUE
  )

####

ehr.with.loinc.parts$final.comp <- ehr.with.loinc.parts$target.term
ehr.with.loinc.parts$final.comp[is.na(ehr.with.loinc.parts$final.comp)] <-
  ehr.with.loinc.parts$iri.from.search[is.na(ehr.with.loinc.parts$final.comp)]

ehr.with.loinc.parts$final.comp <-
  sub(
    pattern = "http://purl.obolibrary.org/obo/",
    replacement = "obo:",
    x = ehr.with.loinc.parts$final.comp,
    fixed = TRUE
  )

ehr.with.loinc.parts$system.iri <-
  gsub(
    pattern = "http://purl.obolibrary.org/obo/",
    replacement = "obo:",
    x = ehr.with.loinc.parts$system.iri,
    fixed = TRUE
  )

at.frame <-
  unique(LoincPartLink[, c('LoincNumber', 'LongCommonName')])

oboflag <-
  grepl(pattern = "^obo:", x = ehr.with.loinc.parts$final.comp)

# excluded because NO automatically-accepted mapping for component
drawing.board <-
  ehr.with.loinc.parts[is.na(ehr.with.loinc.parts$final.comp) |
                         (!oboflag), ]

####    ####    ####    ####

# get labels from ontology

turbo.ontolgy.path <-
  "/Users/markampa/Turbo-Ontology/ontologies/animals_robot_transition/turbo_merged.ttl"

turbo.ontolgy.model <-
  load.rdf(turbo.ontolgy.path, format = "TURTLE")

turbo.labels <- sparql.rdf(
  turbo.ontolgy.model,
  "select * where { ?s <http://www.w3.org/2000/01/rdf-schema#label> ?o }"
)

####    ####    ####    ####

ready.for.robot <-
  unique(ehr.with.loinc.parts[(!is.na(ehr.with.loinc.parts$final.comp)) &
                                oboflag , c(
                                  "LOINC",
                                  "COMPONENT",
                                  "analyte.core.vals",
                                  "final.comp",
                                  "COMPONENT.vals",
                                  "analyte.vals",
                                  "search.counts",
                                  "search.vals",
                                  true.challenges,
                                  pk.in.challenge.slot,
                                  divisors,
                                  serology.suffixes,
                                  "METHOD",
                                  loinc.methods,
                                  "PROPERTY",
                                  "SCALE",
                                  "SYSTEM",
                                  "system.iri",
                                  "TIME"
                                )])

ready.for.robot <- ready.for.robot[order(ready.for.robot$LOINC), ]

robot.row.count <- nrow(ready.for.robot)
all.blanks <- rep(x = '', robot.row.count)

term.id <- paste0('lobo:LOBO_', 1:robot.row.count)

has.curation.status <- all.blanks

# alternative term ... LOINC assay label
at.index <- match(ready.for.robot$LOINC, at.frame$LoincNumber)
alternative.term <- at.frame$LongCommonName[at.index]

#### system display name
sys.disp.name <- match(ready.for.robot$SYSTEM , Part$PartNumber)
sys.disp.name <- Part$PartDisplayName[sys.disp.name]

table(sys.disp.name)

# PK
# ANY MODERNIZASTION REQUIRED FOR AXIOMS?

pk.pastable <- tl.augmenter(pk.in.challenge.slot, 'COMPONENT')
pk.pastable[is.na(pk.pastable)] <- ''
pk.pastable[nchar(pk.pastable) > 0] <-
  paste0(pk.pastable[nchar(pk.pastable) > 0], ' ')
table(pk.pastable)

# todo REAL CHALLENGES (labels and axioms)
abtype <- tl.augmenter(serology.suffixes, 'COMPONENT')
abtype[is.na(abtype)] <- ''
abtype[nchar(abtype) > 0] <-
  paste0(abtype[nchar(abtype) > 0], ' against ')
table(abtype)

# todo (serology) SUFFIXES  (labels and axioms)
challengetype <- tl.augmenter(true.challenges, 'COMPONENT')
challengetype[is.na(challengetype)] <- ''
challengetype[nchar(challengetype) > 0] <-
  paste0(', ', challengetype[nchar(challengetype) > 0], ' challenge ')
table(challengetype)

#### TERM LABEL
# TODO NEW METHODS FOR PK SPECIFICTY
term.label <-
  paste0(
    'Quantitative assay for ',
    abtype,
    pk.pastable,
    ready.for.robot$analyte.core.vals,
    ' in ',
    sys.disp.name,
    challengetype
  )

label.table <- table(term.label)
label.table <- cbind.data.frame(names(label.table), as.numeric(label.table))
colnames(label.table) <- c("label", "count")
table(label.table$count)

term.label[ready.for.robot$TIME == "LP6924-7"] <-
  paste0(term.label[ready.for.robot$TIME == "LP6924-7"], " with 24-hour sample collection")

#### DEAL.. swtich to tl.augmenter approach for labels

# divisors
pastable <- tl.augmenter(divisors, 'COMPONENT')
pastable[is.na(pastable)] <- ''
pastable[nchar(pastable) > 0] <-
  paste0(", normalized by ", pastable[nchar(pastable) > 0], ' ')
table(pastable)

term.label <-
  paste0(term.label, pastable)

label.table <- table(term.label)
label.table <- cbind.data.frame(names(label.table), as.numeric(label.table))
colnames(label.table) <- c("label", "count")
table(label.table$count)


associated.axioms <- all.blanks
associated.axioms[!is.na(ready.for.robot$PartName.divisor) &
                    ready.for.robot$PartName.divisor == 'Creatinine'] <-
  "'has part' some 'division by creatinine concentration normalization'"
associated.axioms[!is.na(ready.for.robot$PartName.divisor) &
                    ready.for.robot$PartName.divisor == '100 leukocytes'] <-
  "'has part' some 'division by 100 leukocytes normalization'"

#### DEAL with methods part of labels... then axiomatic approach

# assume for now that only a single method is asserted for a given test
# but start thinking about mutiple SPECIFIC method cases
# can collapse less specific methods like glucose into 20 mg glucose
# should this be folded into teh split function above?

pastable <- tl.augmenter(loinc.methods, 'METHOD')
pastable[is.na(pastable)] <- ''
pastable[nchar(pastable) > 0] <-
  paste0(' (',  pastable[nchar(pastable) > 0], ')')
table(pastable)

term.label <- paste0(term.label, pastable)

label.table <- table(term.label)
label.table <- cbind.data.frame(names(label.table), as.numeric(label.table))
colnames(label.table) <- c("label", "count")
table(label.table$count)

####

# properties/units part of labels

temp <-
  match(ready.for.robot$PROPERTY,
        turbo_loinc_properties_reviewed$PartTypeVal)
temp <- turbo_loinc_properties_reviewed$obo.iri[temp]

setdiff(temp, turbo.labels[, 1])

temp <- match(temp, turbo.labels[,1])

table(temp, useNA = 'always')

temp <- turbo.labels[temp,2]

term.label <- paste0(term.label, ' [', temp, ']')

label.table <- table(term.label)
label.table <- cbind.data.frame(names(label.table), as.numeric(label.table))
colnames(label.table) <- c("label", "count")
table(label.table$count)


####

definition <- all.blanks

# definition.source <-
#   rep(x = 'Mark Andrew Miller, ORCID:0000-0001-9076-6066|Chris Stoeckert', robot.row.count)

# maybe a link to the script that populated the template?
definition.source <- all.blanks

example.of.usage <- all.blanks
editor.note <- all.blanks
term.editor <-
  rep(x = 'Mark Andrew Miller, ORCID:0000-0001-9076-6066|Chris Stoeckert', robot.row.count)
ontology.term.requester <- all.blanks
term.tracker.item <-
  rep(x = 'ghi:1153', robot.row.count)

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
material.processing.technique[!is.na(ready.for.robot$TIME) &
                                ready.for.robot$TIME == 'LP6924-7'] <-
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

# datums may need to be created and defined in terms of potentiatlly new units.	TODO
output <- rep(x = 'obo:IAO_0000032', robot.row.count)

# move non-molecular entity analytes into here V
target.entity <- all.blanks

# will probably be using this soon
objective <- all.blanks

# or LOINC  https://loinc.org/2339-0/ style ?
database.cross.reference.IRI <-
  paste0('loinc:',
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

temp <- make.table.frame(ready.for.robot$term.label)
temp <- temp$value[temp$count > 1]
ready.for.robot <- ready.for.robot[!ready.for.robot$term.label %in% temp,]

write.table(
  ready.for.robot,
  file = 'lobo_template_headerless.csv',
  row.names = FALSE,
  col.names = FALSE,
  sep = ','
)

label.table[] <- lapply(label.table[], as.character)
ready.for.robot <- base::merge(x = ready.for.robot, y = label.table, by.x = 'term.label', by.y = 'label')

###   ###   ###

# # reworked string inputs... no longer contains anything but unmapped core analyte components
#
# bulk.systems <-
#   kept.best[kept.best$ont.name %in% c('cl', 'uberon') &
#               kept.best$part.number %in% Part$PartNumber[Part$PartTypeName == 'SYSTEM'] &
#               kept.best$cosine == 0 , ]
#
# # see turbo_loinc_systems_reviewed above

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

# now run execute_template.sh

make.table.frame <- function(my.vector) {
  temp <- table(my.vector)
  temp <- cbind.data.frame(names(temp), as.numeric(temp))
  colnames(temp) <- c('value', 'count')
  temp$value <- as.character(temp$value)
  return(temp)
}
# 
# temp <- make.table.frame(ehr.with.loinc.parts$PROPERTY)
# temp <- base::merge(x = temp, y = Part, by.x = 'value', by.y = 'PartNumber')
# write.csv(temp, "ehr_loinc_properties.csv")
