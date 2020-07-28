library(devtools)

# requires a properly formatted "turbo_R_setup.yaml" in the home directory of the user who started this script
# see https://gist.github.com/turbomam/a3915d00ee55d07510493a9944f96696 for template
devtools::source_gist(id = "https://gist.github.com/turbomam/f082295aafb95e71d109d15ca4535e46",
                      sha1 = "c587c0b4ad1ce4dea1816b81d29f08843564d50e",
                      filename = "turbo_R_setup.R")

source("tl_augmenter.R")

do.bp.search <- FALSE
save.bp.search <- FALSE
do.ols.search <- FALSE
save.ols.search <- FALSE

# look for cells by surface antigen
# look for combination analytes


# COVID-19 tests from EHR

# | PK_LAB_RESULT_ITEM_ID | LAB_RESULT_ITEM_CODE | LAB_RESULT_ITEM_DESCRIPTION | LOINC   | PDS mentions |
# | --------------------- | -------------------- | --------------------------- | ------- | -----------: |
# | 20669003              | SARS-COV-2 IGG       | SARS-COV-2 IGG              | 94563-4 |          621 |

# from LoincPartLink

# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | LoincNumber    | PartNumber    | PartName                     | PartCodeSystem      | PartTypeName    | LinkTypeName         | Property                                  |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP417915-8    | SARS coronavirus 2 Ab.IgG    | http://loinc.org    | COMPONENT       | Primary              | http://loinc.org/property/COMPONENT       |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP417915-8    | SARS coronavirus 2 Ab.IgG    | http://loinc.org    | COMPONENT       | DetailedModel        | http://loinc.org/property/analyte         |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP417540-4    | SARS coronavirus 2           | http://loinc.org    | COMPONENT       | *SyntaxEnhancement*  | http://loinc.org/property/analyte-core    |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP14672-7     | IgG                          | http://loinc.org    | COMPONENT       | Search               | http://loinc.org/property/search          |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP16680-8     | Coronavirus                  | http://loinc.org    | COMPONENT       | Search               | http://loinc.org/property/search          |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP35807-4     | SARS coronavirus             | http://loinc.org    | COMPONENT       | Search               | http://loinc.org/property/search          |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+
# | 94563-4        | LP417914-1    | SARS coronavirus 2 Ab        | http://loinc.org    | COMPONENT       | Search               | http://loinc.org/property/search          |
# +----------------+---------------+------------------------------+---------------------+-----------------+----------------------+-------------------------------------------+

# can get part details from LOINC FHIR... but isn't all of this included in the CSV files?

# sars.cov2.fhir <-
#   GET(
#     "https://fhir.loinc.org/CodeSystem/$lookup?code=LP417540-4&system=http%3A%2F%2Floinc.org&_format=json",
#     authenticate("markampa", "ref1ned&Marked", type = "basic")
#   )
#
# sars.cov2.fhir <-
#   fromJSON(rawToChar(sars.cov2.fhir$content))$parameter

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

# HAVEN'T APPLIED AXIOMS YET
loinc.methods <-
  c(
    # 'Automated',
    # 'Manual',
    'Automated count',
    'Calculated',
    'Confirm',
    'Dialysis',
    'Direct assay',
    'Electrophoresis',
    'HA',
    'High sensitivity',
    'IA',
    'IB',
    'IF',
    'ISE',
    'Manual count',
    'Oximetry',
    'Screen',
    'SSA',
    'Test strip',
    'Test strip.automated',
    'Creatinine-based formula (MDRD)',
    'Creatinine-based formula (CKD-EPI)',
    'Detection limit <= 0.05 mIU/L',
    "With P-5'-P"
  )


# loinc:14308-1 amphetamine screen

#### PUT FILENAMES/PATHS IN CONFIG
harmonized_ontology_rankings <-
  read.csv(config$harmonized_ontology_rankings, stringsAsFactors = FALSE)

### CSV reads
### some from LOINC, some manually curated

LoincPartLink_Primary <-
  read_csv(paste0(config$loinc.root, config$LoincPartLink_Primary.file))

LoincPartLink_Supplementary <-
  read_csv(paste0(config$loinc.root, config$LoincPartLink_Supplementary.file))

LoincPartLink <-
  rbind.data.frame(LoincPartLink_Primary, LoincPartLink_Supplementary)

# additional part details
Part <-
  read_csv(paste0(config$loinc.root, config$Part.file))

# get LOINC provided component mappings
PartRelatedCodeMapping <-
  read_csv(paste0(config$loinc.root, config$PartRelatedCodeMapping.file))

# # THIS IS HANDY FOR DIGNOSISNG/PRIORITIZING PATTERNS FOR MANUAL CURATION
MultiAxialHierarchy <-
  read_csv(paste0(config$loinc.root, config$MultiAxialHierarchy.file))

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
  turbo_loinc_properties_reviewed[turbo_loinc_properties_reviewed$use == 'y' ,]

###   ###   ###

# manually curate LOINC system-OBO term mappings?
# other options:
# BioPortal mappings or term lookups,
# OLS,
# local Solr
# OntoBee API?
# LOINC provided ChEBI mappings

# they are  all single UBERON terms except
# LP7576-4 Ser/Plas http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969
# LP7579-8 Ser/Plas/Bld http://purl.obolibrary.org/obo/UBERON_0001977|http://purl.obolibrary.org/obo/UBERON_0001969|http://purl.obolibrary.org/obo/UBERON_0000178
# LP7536-8 RBC http://purl.obolibrary.org/obo/CL_0000232
# LP7720-8 WBC http://purl.obolibrary.org/obo/CL_0000738

# bootstrapped from earlier iterations of this script
# but systems are not searched any more, only components
# could easily be reactivated

turbo_loinc_systems_reviewed <-
  read_csv(config$reviewed.systems.file)
turbo_loinc_systems_reviewed$use[is.na(turbo_loinc_systems_reviewed$use)] <-
  'n'

turbo_loinc_systems_reviewed <-
  turbo_loinc_systems_reviewed[turbo_loinc_systems_reviewed$use == 'y' , ]

###   ###   ###

print(nrow(ehr_loinc_counts))

common.loincs <-
  ehr_loinc_counts$LOINC[ehr_loinc_counts$LOINC_COUNT >= config$common.threshold]

print(length(common.loincs))

###   ###   ###

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

print(nrow(relevant.primary.part.cast))

# ~ 25 bad or old LOINCs in PDS?
setdiff(common.loincs, relevant.primary.part.cast$LoincNumber)

# bad.loinc.freqs <-
#   ehr_loinc_counts[ehr_loinc_counts$LOINC %in% setdiff(common.loincs, relevant.primary.part.cast$LoincNumber),]
# R_LAB_RESULT_ITEM <-
#   read.csv("~/loinc_in_obo/R_LAB_RESULT_ITEM_202005050841.csv",
#            stringsAsFactors = FALSE)
# bad.loinc.freqs <- base::merge(x = bad.loinc.freqs, y = R_LAB_RESULT_ITEM)
# write.csv(bad.loinc.freqs, file =  "bad_pds_loinc_code_usage.csv")

###   ###   ###

# first EHR relevance/LOINC content merge!

ehr.with.loinc.parts <-
  left_join(x = ehr_loinc_counts,
            y = relevant.primary.part.cast,
            by = c("LOINC" = "LoincNumber"))

print(nrow(ehr.with.loinc.parts))

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

# LP6960-1
# Pt
# Point in time (spot)
#
# LP6924-7
# 24H
# 24 hours

# scale choice:
# Qn SCALE 5x more entries than next most commons, Nom and Ord
# LP7753-9
# Qn
# Quantitative

# for 94563-4SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay
# LP7751-3
# Ord

###   ###   ###

# start merging in reviewed time, system, scale and property choices
# see csv reads at top of script
# put these ACCEPTED SCALE AND TIME values in config?
# we will have to get to additional scales and times in the future
# may require new ontology terms or new mapping patterns

lost.due.to.time.scale <- ehr.with.loinc.parts$LOINC
ehr.with.loinc.parts <-
  ehr.with.loinc.parts[ehr.with.loinc.parts$PROPERTY %in% turbo_loinc_properties_reviewed$PartTypeVal &
                         ehr.with.loinc.parts$SCALE %in% c("LP7753-9", "LP7751-3") &
                         ehr.with.loinc.parts$SYSTEM %in% turbo_loinc_systems_reviewed$PartTypeVal &
                         (ehr.with.loinc.parts$TIME == "LP6960-1" |
                            ehr.with.loinc.parts$TIME == "LP6924-7") , ]

lost.due.to.time.scale <-
  setdiff(lost.due.to.time.scale, ehr.with.loinc.parts$LOINC)

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

# component properties
# "PartCodeSystem"
# "PartTypeName"
# "LinkTypeName"
# "Property"

unique.loincs <- unique(ehr.with.loinc.parts$LOINC)
unique.loincs <-
  LoincPartLink[LoincPartLink$LoincNumber %in% unique.loincs &
                  LoincPartLink$PartTypeName == "COMPONENT" , c("PartNumber",
                                                                "PartName")]

original.needs.component.mapping <- unique(unique.loincs$PartNumber)

current.needs.component.mapping <- original.needs.component.mapping
current.component.mapping.complete <- c()
current.component.mapping.frame <- matrix(ncol = 6)
colnames(current.component.mapping.frame) <-
  c(
    'loinc.part.code',
    'loinc.part.name',
    'obo.uri',
    'obo.label',
    'rank',
    'justification'
  )

###   ###   ###

# start with component mappings that we have reviewed

previously.reviewed  <- read_csv("core_analyte_success.csv")

previously.reviewed <-
  unique(previously.reviewed[, c("PartNumber",
                                 "PartName" ,
                                 "obo.uri" ,
                                 "obo.label" ,
                                 "rank" ,
                                 "justification")])

justification.reviewed <-
  paste0(previously.reviewed$justification, ", reviewed")

previously.reviewed$justification <- justification.reviewed

check.reviewed <- make.table.frame(previously.reviewed$PartNumber)
check.reviewed <- check.reviewed$value[check.reviewed$count == 1]
previously.reviewed <-
  previously.reviewed[previously.reviewed$PartNumber %in% check.reviewed ,]

update.accounting(
  data = previously.reviewed ,
  loinc.part.code =  "PartNumber",
  loinc.part.name = "PartName" ,
  obo.uri = "obo.uri" ,
  obo.label = "obo.label" ,
  rank = "rank" ,
  justification = previously.reviewed$justification
)

current.component.mapping.frame <-
  current.component.mapping.frame[!is.na(current.component.mapping.frame$loinc.part.code),]

###   ###   ###

# LOINC provides some mappings itself

relevant.links <-
  PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% current.needs.component.mapping ,]

obo.friendly.sources <-
  c("https://www.ebi.ac.uk/chebi",
    "https://www.ncbi.nlm.nih.gov/taxonomy")

obo.friendly.links <-
  relevant.links[relevant.links$ExtCodeSystem %in% obo.friendly.sources ,]

###   ###   ###

friend.of.a.friend <-
  c(
    "http://pubchem.ncbi.nlm.nih.gov",
    "http://www.nlm.nih.gov/research/umls/rxnorm"
  )

friend.of.a.friend <-
  relevant.links[relevant.links$ExtCodeSystem %in% friend.of.a.friend ,]

friend.of.a.friend <-
  friend.of.a.friend[!(friend.of.a.friend$PartNumber %in% obo.friendly.links$PartNumber) , ]

pubchem.cids <-
  friend.of.a.friend$ExtCodeId[friend.of.a.friend$ExtCodeSystem == "http://pubchem.ncbi.nlm.nih.gov"]

pubchem.xrefs <-
  lapply(
    X = pubchem.cids,
    FUN = function(current.cid) {
      # current.cid <- 769
      print(current.cid)
      current.url <-
        paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
          current.cid ,
          "/xrefs/SBURL/JSON"
        )
      # print(current.url)
      current.xrefs <- GET(current.url)
      current.xrefs <- rawToChar(current.xrefs$content)
      current.xrefs <- fromJSON(current.xrefs)
      current.xrefs <-
        current.xrefs$InformationList$Information$SBURL
      current.xrefs <- unlist(current.xrefs)
      current.xrefs <-
        cbind.data.frame(rep(x = current.cid, n = length(current.xrefs)), current.xrefs)
      # print(current.xrefs)
      return(current.xrefs)
    }
  )

pubchem.xrefs <- do.call(rbind.data.frame, pubchem.xrefs)

colnames(pubchem.xrefs) <- c("ExtCodeId", "xref.url")

pubmed.to.chebi <-
  pubchem.xrefs[grepl(pattern = "chebi", x = pubchem.xrefs$xref.url) ,]

if (nrow(pubmed.to.chebi) > 0) {
  pubmed.to.chebi$next.step <-
    sub(
      pattern = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=",
      replacement = "",
      x = pubmed.to.chebi$xref.url,
      fixed = TRUE
    )
  
  pubmed.to.chebi$ExtCodeSystem <- "http://pubchem.ncbi.nlm.nih.gov"
  
  pubmed.to.chebi <-
    pubmed.to.chebi[c("ExtCodeSystem", "ExtCodeId", "next.step")]
  
  reunited <-
    base::merge(x = friend.of.a.friend, y = pubmed.to.chebi)
  
  reunited$ExtCodeSystem <- "https://www.ebi.ac.uk/chebi"
  reunited$ExtCodeId <- reunited$next.step
  
  reunited <- reunited[, setdiff(colnames(reunited), "next.step")]
  
  reunited$ExtCodeDisplayName <- ""
  reunited$Equivalence <- ""
  
  reunited <- rbind.data.frame(obo.friendly.links, reunited)
} else {
  reunited <- obo.friendly.links
}

reunited$ext.iri <- reunited$ExtCodeId

reunited$ext.iri[reunited$ExtCodeSystem == "https://www.ebi.ac.uk/chebi"] <-
  sub(pattern = '^CHEBI:',
      replacement = 'http://purl.obolibrary.org/obo/CHEBI_',
      x = reunited$ext.iri[reunited$ExtCodeSystem == "https://www.ebi.ac.uk/chebi"])

reunited$ext.iri[reunited$ExtCodeSystem == "https://www.ncbi.nlm.nih.gov/taxonomy"] <-
  paste0("http://purl.obolibrary.org/obo/NCBITaxon_",
         reunited$ext.iri[reunited$ExtCodeSystem == "https://www.ncbi.nlm.nih.gov/taxonomy"])

reunited$rank <- NA

reunited.check <- make.table.frame(reunited$PartNumber)
reunited.check <- reunited.check$value[reunited.check$count == 1]
reunited <- reunited[reunited$PartNumber %in% reunited.check ,]

update.accounting(
  data = reunited,
  loinc.part.code = "PartNumber",
  loinc.part.name = "PartName",
  obo.uri = "ext.iri",
  obo.label = "ExtCodeDisplayName",
  rank = "rank",
  justification = "LOINC provided mappings"
)

###   ###   ###

#### begin BioPortal mappings retrieval
# how long should we expect this to take?

if (do.bp.search) {
  bp.retreived.and.parsed <-
    bp.map.retreive.and.parse(sort(unique(current.needs.component.mapping)))
  if (save.bp.search) {
    save(bp.retreived.and.parsed, file = "bp_retreived_and_parsed.Rdata")
  }
} else {
  load("bp_retreived_and_parsed.Rdata")
}


bp.retreived.and.parsed <-
  do.call(rbind.data.frame, bp.retreived.and.parsed)

bp.retreived.and.parsed[] <-
  lapply(bp.retreived.and.parsed[], as.character)

# who did we get mapped to?
rap.tab <- table(bp.retreived.and.parsed$target.ontology)
rap.tab <- cbind.data.frame(names(rap.tab), as.numeric(rap.tab))
names(rap.tab) <- c("target.ontology", "map.count")

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

acceptable.bp.retreived.and.parsed$target.label <- ""

name.and.num <- Part[, c("PartNumber", "PartName")]

acceptable.bp.retreived.and.parsed <-
  left_join(x = acceptable.bp.retreived.and.parsed,
            y = name.and.num,
            by = c("source.id" = "PartNumber"))

acceptable.bp.retreived.and.parsed$rank <- NA

check.bp <-
  make.table.frame(acceptable.bp.retreived.and.parsed$source.term)
check.bp <- check.bp$value[check.bp$count == 1]
acceptable.bp.retreived.and.parsed <-
  acceptable.bp.retreived.and.parsed[acceptable.bp.retreived.and.parsed$source.term %in% check.bp ,]

update.accounting(
  data = acceptable.bp.retreived.and.parsed,
  loinc.part.code = "source.id",
  loinc.part.name = "PartName",
  obo.uri = "target.term",
  obo.label = "target.label",
  rank = "rank",
  justification = "BioPortal mapping"
)

###   ###   ###

text.search.input <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs &
                  LoincPartLink$PartNumber %in% current.needs.component.mapping &
                  LoincPartLink$Property == "http://loinc.org/property/analyte-core" , c("PartNumber", "PartName")]

text.search.input <- unique(text.search.input)

text.search.input$singularized <- sub(pattern = 's$',
                                      replacement = '',
                                      x = text.search.input$PartName)

original <- text.search.input[, c('PartNumber', 'PartName')]
singularized <- text.search.input[, c('PartNumber', 'singularized')]
colnames(singularized) <- colnames(original)
text.search.input <-
  unique(rbind.data.frame(original, singularized))
colnames(text.search.input) <- c('PartNumber', 'any.name')

text.search.input <-
  text.search.input[order(text.search.input$any.name), ]

# ###   ###   ###

###  have we already searched for any of these?
table(text.search.input$PartNumber %in% current.component.mapping.complete)
table(text.search.input$any.name %in% Part$PartName[Part$PartNumber %in% current.component.mapping.complete])

# substitute '.' with ' ' for neutrophils.immature etc.
# substitute 'spp$' or 'sp$' with ''  for genus-level NCBI taxon entities

onto.list.for.ols <-
  paste0(sort(tolower(harmonized_ontology_rankings$rhs)), collapse = ',')

if (do.ols.search) {
  ols.attempts <-
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
            # kept.row.count = 600,
            kept.row.count = 200,
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
  
  col.counts <- sapply(ols.attempts, length)
  table(col.counts)
  
  tenners <- ols.attempts[col.counts == 10]
  tenners <- do.call(rbind.data.frame, tenners)
  tenners$synonym <- ""
  
  elevens <- ols.attempts[col.counts == 11]
  elevens <- do.call(rbind.data.frame, elevens)
  
  ols.bound <- rbind.data.frame(elevens, tenners)
  if (save.ols.search) {
    save(ols.bound, file = "ols_bound.Rdata")
  }
} else {
  load("ols_bound.Rdata")
}

# requires that the query matches some field from OLS exactly
# could go back to string distances or something more complex

ols.deep.check <-
  apply(
    X = ols.bound,
    MARGIN = 1,
    FUN = function(current.row) {
      # current.row <- ols.attempts[1,]
      print(current.row$query)
      temp <-
        c(
          unlist(current.row$label),
          unlist(current.row$synonym),
          unlist(current.row$annotations_trimmed)
        )
      # print(temp)
      temp <- sort(unique(tolower(temp)))
      temp <- sort(unique(unlist(lapply(temp, tolower))))
      temp.cosines <-
        stringdist::stringdist(a = tolower(current.row$query),
                               b = temp,
                               method = "cosine")
      best.match <- temp[which.min(temp.cosines)]
      best.cosine <- temp.cosines[which.min(temp.cosines)]
      # print(temp)
      return(
        list(
          iri = current.row[["iri"]],
          short_form = current.row[["short_form"]],
          obo_id = current.row[["obo_id"]],
          label = current.row[["label"]],
          ontology_name = current.row[["ontology_name"]],
          ontology_prefix = current.row[["ontology_prefix"]],
          query = current.row[["query"]],
          loinc.part = current.row[["loinc.part"]],
          rank = current.row[["rank"]],
          best.match = best.match,
          best.cosine = best.cosine
        )
      )
    }
  )

col.counts <- sapply(ols.deep.check, length)
table(col.counts)

elevens <- ols.deep.check[col.counts == 11]

ols.deep.check <- do.call(rbind.data.frame, elevens)

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

kept.best$loinc.vs.label <-
  stringdist(kept.best$query, kept.best$label, method = 'cosine')

# create a combination scoer? with that variables and what operations?

####    ####    ####    ####

lowest.best <-
  aggregate(kept.best$best.cosine,
            by = list(kept.best$loinc.part),
            FUN = min)

names(lowest.best) <- c('loinc.part', 'best.cosine')

lowest.best <- base::merge(
  x = kept.best,
  y = lowest.best,
  by.x = c('loinc.part', 'best.cosine'),
  by.y = c('loinc.part', 'best.cosine')
)

####    ####    ####    ####

lowest.label <-
  aggregate(kept.best$loinc.vs.label,
            by = list(kept.best$loinc.part),
            FUN = min)

names(lowest.label) <- c('loinc.part', 'loinc.vs.label')

lowest.label <- base::merge(
  x = kept.best,
  y = lowest.label,
  by.x = c('loinc.part', 'loinc.vs.label'),
  by.y = c('loinc.part', 'loinc.vs.label')
)

lowest.both <- unique(rbind.data.frame(lowest.best, lowest.label))

codes.paths <-
  MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% lowest.both$loinc.part , c("CODE", "PATH_TO_ROOT")]

labeled.paths <-
  sapply(codes.paths$PATH_TO_ROOT, function(current.path) {
    print(current.path)
    temp <-
      unlist(strsplit(current.path, split = '.', fixed = TRUE))
    temp <- match(temp, Part$PartNumber)
    temp <- Part$PartName[temp]
    temp <- paste0(temp, collapse = " -> ")
    return(temp)
  })
codes.paths$labeled.paths <- labeled.paths

lowest.both <-
  base::merge(x = lowest.both,
              y = codes.paths,
              by.x = 'loinc.part' ,
              by.y = 'CODE')

lowest.both$plus <-
  grepl(pattern = '+',
        x = lowest.both$query,
        fixed = TRUE)

lowest.both$period <-
  grepl(pattern = '.',
        x = lowest.both$query,
        fixed = TRUE)

average.cosine <-
  apply(
    lowest.both,
    MARGIN = 1,
    FUN = function(current.row) {
      # print(current.row)
      return((as.numeric(current.row[['loinc.vs.label']]) + as.numeric(current.row[['best.cosine']])) /
               2)
    }
  )

lowest.both$average.cosine <- average.cosine

lowest.both <-
  lowest.both[, c(
    "loinc.part",
    "query",
    "labeled.paths",
    'plus',
    'period' ,
    "obo_id",
    "domain",
    "priority",
    "rank",
    "label",
    "best.match" ,
    "average.cosine",
    "iri",
    "ontology_prefix",
    "short_form",
    "ontology_name",
    "loinc.vs.label" ,
    "best.cosine"
  )]

names(lowest.both) <- c(
  "loinc.part",
  "part.name",
  "labeled.paths",
  'plus',
  'period',
  "obo_id",
  "domain",
  "ontology.priority",
  "ols.search.rank",
  "obo.label",
  "best.match" ,
  "average.cosine",
  "iri",
  "ontology_prefix",
  "short_form",
  "ontology_name",
  "loinc.vs.label.cosine" ,
  "best.cosine"
)

# write.csv(lowest.both, file = "loinc_ols_results.csv", row.names = FALSE)

####    ####    ####    ####

accepted.ols.res <-
  unique(lowest.both[lowest.both$loinc.vs.label.cosine == 0, c(
    "loinc.part",
    "part.name",
    "obo.label",
    "iri",
    "ontology.priority",
    "loinc.vs.label.cosine" ,
    "average.cosine"
  )])

if (nrow(accepted.ols.res) > 0) {
  lowest.priority <-
    aggregate(
      accepted.ols.res$ontology.priority,
      by = list(accepted.ols.res$loinc.part),
      FUN = min
    )
  colnames(lowest.priority) <- c('loinc.part', 'ontology.priority')
  
  preferred.ols <-
    base::merge(x = accepted.ols.res , y = lowest.priority)
  
  update.accounting(
    data = preferred.ols,
    loinc.part.code = "loinc.part",
    loinc.part.name = "part.name",
    obo.uri = "iri",
    obo.label = "obo.label",
    rank = "average.cosine",
    justification = "ols text search exact matches"
  )
}

# save.image("20200626.Rdata")
# write.csv(current.component.mapping.frame, "20200626.csv")
# manual review 1: review and possible reject some of the OLS mappings
# manual review 2:
#  current.component.mapping.frame will contain several kinds of component terms,
#  like core analytes and search term
#  need to merge the loinc number back in and pick the most detailed component type

current.component.mapping.frame$obo.uri <-
  sub(pattern = "http://purl.bioontology.org/ontology/NCBITAXON/",
      replacement = "	http://purl.obolibrary.org/obo/NCBITaxon_",
      x = current.component.mapping.frame$obo.uri)

preferred.component.property <-
  LoincPartLink[LoincPartLink$LoincNumber %in% common.loincs ,]

preferred.component.property <-
  inner_join(
    x = preferred.component.property,
    y =  current.component.mapping.frame,
    by = c("PartNumber" = "loinc.part.code")
  )

table(preferred.component.property$Property)

component_property_priorities <-
  structure(
    list(
      property = c(
        "http://loinc.org/property/analyte-core",
        "http://loinc.org/property/analyte",
        "http://loinc.org/property/COMPONENT",
        "http://loinc.org/property/search"
      ),
      priority = c(1, 2, 3, 4)
    ),
    class = c("data.frame"),
    row.names = c(NA,-4L),
    spec = structure(list(
      cols = list(
        property = structure(list(), class = c("collector_character",
                                               "collector")),
        priority = structure(list(), class = c("collector_double",
                                               "collector"))
      ),
      default = structure(list(), class = c("collector_guess",
                                            "collector")),
      skip = 1
    ), class = "col_spec")
  )

preferred.component.property <-
  inner_join(x = preferred.component.property,
             y =  component_property_priorities,
             by = c("Property" = "property"))

lowest.priority <-
  aggregate(
    preferred.component.property$priority,
    by = list(preferred.component.property$LoincNumber),
    FUN = min
  )

colnames(lowest.priority) <- c('LoincNumber', 'priority')

preferred.component.property <-
  base::merge(x = preferred.component.property , y = lowest.priority)

table(preferred.component.property$Property)

core.analyte.success <-
  preferred.component.property[preferred.component.property$Property == "http://loinc.org/property/analyte-core" , ]

check.core.analytes <-
  make.table.frame(core.analyte.success$LoincNumber)

#### add filtering for singles

# write.csv(core.analyte.success, file = "core_analyte_success.csv", row.names = FALSE)

no.core.analyte <-
  preferred.component.property[preferred.component.property$Property != "http://loinc.org/property/analyte-core" , ]

# search.term.promiscuity <- table(no.core.analyte$obo.uri)
search.term.promiscuity <- make.table.frame(no.core.analyte$obo.uri)
names(search.term.promiscuity) <-
  c("obo.uri", "search.term.promiscuity")
no.core.analyte <-
  base::merge(x = no.core.analyte, y = search.term.promiscuity, by = "obo.uri")

no.core.analyte$assay.vs.label <-
  stringdist(no.core.analyte$LongCommonName,
             no.core.analyte$obo.label,
             method = 'cosine')

no.core.analyte <-
  no.core.analyte[, c(
    "LoincNumber",
    "LongCommonName",
    "PartNumber" ,
    "PartName",
    "obo.uri",
    "obo.label",
    "assay.vs.label",
    "LinkTypeName",
    "justification",
    "search.term.promiscuity"
  )]

# write.csv(x = no.core.analyte,
#           file = "no_core_analyte_found.csv",
#           row.names = FALSE)

####

ehr.with.loinc.parts <-  base::merge(
  x = ehr.with.loinc.parts,
  y = core.analyte.success,
  by.x = 'LOINC',
  by.y = 'LoincNumber',
  all.x = TRUE
)

ehr.with.loinc.parts <-
  ehr.with.loinc.parts[!is.na(ehr.with.loinc.parts$obo.uri),]

conflicting.rows <- make.table.frame(ehr.with.loinc.parts$LOINC)
singleton.loincs <-
  conflicting.rows$value[conflicting.rows$count == 1]
ehr.with.loinc.parts <-
  ehr.with.loinc.parts[ehr.with.loinc.parts$LOINC %in% singleton.loincs , ]

# ####    ####    ####    ####

### todo add systems back in

system.list <-
  match(x = ehr.with.loinc.parts$SYSTEM, table = turbo_loinc_systems_reviewed$PartTypeVal)
system.list <- turbo_loinc_systems_reviewed$label[system.list]
ehr.with.loinc.parts$system.phrase <- system.list

####

# make sure we're reporting ( write.csv , save.image )
# all assays and components that haven't been included
# in addition to the more detailed exports above

####    ####    ####    ####

turbo.url <-
  "https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/robot_implementation/turbo_merged.ttl"

ont.tf <- tempfile()
download.file(turbo.url, ont.tf)

turbo.model <-
  rrdf::load.rdf(filename = ont.tf, format = "TURTLE")

turbo.label.query <-
  "select
distinct
?s ?o
where
{
?s <http://www.w3.org/2000/01/rdf-schema#label> ?prel .
# bind(lcase(str(?prel)) as ?o)
bind(str(?prel) as ?o)
filter(isiri(?s))
}"
  
  Sys.time()
  system.time(turbo.labels <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.label.query))
  
  turbo.labels <- as.data.frame(turbo.labels)
  
  check.turbo.labels <- make.table.frame(turbo.labels$s)
  check.turbo.labels <-
    check.turbo.labels$value[check.turbo.labels$count == 1]
  
  turbo.labels <-
    turbo.labels[turbo.labels$s %in% check.turbo.labels ,]
  
  ####    ####    ####    ####
  
  expected.cols <- c(
    "LOINC",
    "COMPONENT",
    "PartName",
    "obo.uri",
    true.challenges,
    pk.in.challenge.slot,
    divisors,
    serology.suffixes,
    "METHOD",
    loinc.methods,
    "PROPERTY",
    "SCALE",
    "SYSTEM",
    "system.phrase",
    "TIME"
  )
  
  available.cols <- colnames(ehr.with.loinc.parts)
  
  # this should be the empty set... meaning all required columns are avaialble
  setdiff(expected.cols, available.cols)
  
  # this dataframe is easier to examine than ehr.with.loinc.parts
  # but isn't really necessary
  # we could just get the columns directly from ehr.with.loinc.parts
  
  pre.ready.for.robot <-
    unique(ehr.with.loinc.parts[(!is.na(ehr.with.loinc.parts$obo.uri)) ,  c(
      "LOINC",
      "COMPONENT",
      "PartName",
      "obo.uri",
      true.challenges,
      pk.in.challenge.slot,
      divisors,
      serology.suffixes,
      "METHOD",
      loinc.methods,
      "PROPERTY",
      "SCALE",
      "SYSTEM",
      "system.phrase",
      "TIME"
    )])
  
  pre.ready.for.robot <-
    pre.ready.for.robot[order(pre.ready.for.robot$LOINC),]
  robot.row.count <- nrow(pre.ready.for.robot)
  all.blanks <- rep(x = '', robot.row.count)
  
  # any of these all blank columns could be overwritten below
  has.curation.status <- all.blanks
  definition <- all.blanks
  example.of.usage <- all.blanks
  editor.note <- all.blanks
  ontology.term.requester <- all.blanks
  logical.note <- all.blanks
  device <- all.blanks
  reagent <- all.blanks
  molecular.label <- all.blanks
  target.entity <- all.blanks
  
  # will probably be using this soon
  objective <- all.blanks
  
  # definition.source <-
  #   rep(x = 'Mark Andrew Miller, ORCID:0000-0001-9076-6066|Chris Stoeckert', robot.row.count)
  # maybe a link to the script that populated the template?
  definition.source <- all.blanks
  
  # skip for now. would be based on method
  detection.technique <- all.blanks
  
  # skip in favor of analyte/target entity for now	some specific kind of datum based on LOINC property.
  input <- all.blanks
  
  evaluant <- pre.ready.for.robot$system.phrase
  
  term.editor <-
    rep(x = 'Mark Andrew Miller, ORCID:0000-0001-9076-6066|Chris Stoeckert', robot.row.count)
  
  term.tracker.item <-
    rep(x = 'ghi:1153', robot.row.count)
  
  # or LOINC  https://loinc.org/2339-0/ style ?
  database.cross.reference.IRI <-
    paste0('loinc:',
           pre.ready.for.robot$LOINC)
  
  # subclass or equivalent
  logical.type <-
    rep(x = 'subclass', robot.row.count)
  
  #  quantitative by default and PT timing by default
  parent.class <-
    rep(x = "'assay of specimen from organism'", robot.row.count)
  
  # 24 collection axiom as above if necessary (LP6924-7 TURBO_0010724)
  # switch and create 24 specimen collection objective specification
  material.processing.technique <- all.blanks
  material.processing.technique[!is.na(pre.ready.for.robot$TIME) &
                                  pre.ready.for.robot$TIME == 'LP6924-7'] <-
    "'specimen collection process' and ( achieves_planned_objective some '24 hour sample collection objective' )"
  
  ####
  
  associated.axioms <- all.blanks
  associated.axioms[!is.na(pre.ready.for.robot$PartName.divisor) &
                      pre.ready.for.robot$PartName.divisor == 'Creatinine'] <-
    "'has part' some 'division by creatinine concentration normalization'"
  associated.axioms[!is.na(pre.ready.for.robot$PartName.divisor) &
                      pre.ready.for.robot$PartName.divisor == '100 leukocytes'] <-
    "'has part' some 'division by 100 leukocytes normalization'"
  
  ####
  
  # properties/units
  prop.unit <-
    match(pre.ready.for.robot$PROPERTY,
          turbo_loinc_properties_reviewed$PartTypeVal)
  prop.unit <-
    turbo_loinc_properties_reviewed[prop.unit,]
  prop.quote.stripped <-
    gsub(pattern = "'",
         replacement = "",
         x = prop.unit$label)
  table(prop.quote.stripped, useNA = 'always')
  
  output <- prop.unit$label
  ax.out <-
    paste0(
      "('scalar measurement datum' and 'has measurement unit label' some ",
      prop.unit$label,
      ")"
    )
  output[prop.unit$scalar] <- ax.out[prop.unit$scalar]
  
  
  # assays for antibodies
  # would assays with RNA suffix go in here, too?
  sero.backtrack <- pre.ready.for.robot[, c(serology.suffixes)]
  sero.backtrack <- sero.backtrack * 1
  rownames(sero.backtrack) <- ehr.with.loinc.parts$LOINC
  sb.rowsums <- rowSums(sero.backtrack, na.rm = TRUE)
  sb.colsums <- colSums(sero.backtrack, na.rm = TRUE)
  names(sb.colsums) <- colnames(sero.backtrack)
  sero.backtrack$rowsums <- sb.rowsums
  sero.backtrack$LOINC <- rownames(sero.backtrack)
  
  # alternative term ... LOINC assay description
  alt.term.frame <-
    unique(LoincPartLink[, c('LoincNumber', 'LongCommonName')])
  at.index <-
    match(pre.ready.for.robot$LOINC, alt.term.frame$LoincNumber)
  alternative.term <- alt.term.frame$LongCommonName[at.index]
  
  #### system display name
  sys.disp.name <-
    match(pre.ready.for.robot$SYSTEM , Part$PartNumber)
  sys.disp.name <- Part$PartDisplayName[sys.disp.name]
  table(sys.disp.name)
  
  # REFACTOR
  # PK pharmacokinetics
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
  
  # divisors
  div.pastable <- tl.augmenter(divisors, 'COMPONENT')
  div.pastable[is.na(div.pastable)] <- ''
  div.pastable[nchar(div.pastable) > 0] <-
    paste0(", normalized by ", div.pastable[nchar(div.pastable) > 0], ' ')
  table(div.pastable)
  
  #### DEAL with methods part of labels... then axiomatic approach
  # assume for now that only a single method is asserted for a given test
  # but start thinking about multiple SPECIFIC method cases
  # can collapse less specific methods like glucose into 20 mg glucose
  # should this be folded into teh split function above?
  
  # this tl.augmenter is failing now (too few entries)
  meth.pastable <- tl.augmenter(loinc.methods, 'METHOD')
  
  meth.pastable[is.na(meth.pastable)] <- ''
  meth.pastable[nchar(meth.pastable) > 0] <-
    paste0(' (',  meth.pastable[nchar(meth.pastable) > 0], ')')
  table(meth.pastable)
  
  ####
  
  # put all components in here and wait for errors from things that aren't molecular entities? V
  analyte <- pre.ready.for.robot$obo.uri
  table(is.na(analyte))
  
  analyte.index <- match(analyte, turbo.labels$s)
  table(is.na(analyte.index))
  
  pre.analyte <- turbo.labels$o[analyte.index]
  table(is.na(pre.analyte))
  
  temp <- cbind.data.frame(analyte, analyte.index)
  temp <- temp[is.na(temp$analyte.index),]
  
  single.quote.prob <- grepl(pattern = "'", x = pre.analyte)
  table(single.quote.prob)
  
  analyte <- paste0("'", pre.analyte, "'")
  
  # inefficient
  analyte[sero.backtrack$rowsums == 1] <-
    paste0("('immunoglobulin complex' and 'has disposition to bind' some '",
           pre.analyte[sero.backtrack$rowsums == 1],
           "')")
  
  analyte[(!(is.na(sero.backtrack$Ab.IgA)))] <-
    paste0(
      "('IgA immunoglobulin complex' and 'has disposition to bind' some '",
      pre.analyte[(!(is.na(sero.backtrack$Ab.IgA)))],
      "')"
    )
  
  analyte[(!(is.na(sero.backtrack$Ab.IgE)))] <-
    paste0(
      "('IgE immunoglobulin complex' and 'has disposition to bind' some '",
      pre.analyte[(!(is.na(sero.backtrack$Ab.IgE)))],
      "')"
    )
  
  analyte[(!(is.na(sero.backtrack$Ab.IgG)))] <-
    paste0(
      "('IgG immunoglobulin complex' and 'has disposition to bind' some '",
      pre.analyte[(!(is.na(sero.backtrack$Ab.IgG)))],
      "')"
    )
  
  analyte[(!(is.na(sero.backtrack$Ab.IgM)))] <-
    paste0(
      "('IgM immunoglobulin complex' and 'has disposition to bind' some '",
      pre.analyte[(!(is.na(sero.backtrack$Ab.IgM)))],
      "')"
    )
  
  prop.lab.component <- rep("Assay", nrow(pre.ready.for.robot))
  prop.lab.component[prop.unit$scalar] <- "Quantitative assay"
  
  ####
  
  # will eventually switch to OBI IDs
  # because these have a static start of 3000001, asn the number of successful assays keeps growing,
  # ID/assay pairings will not be constant
  id.prefix <- "turbo:TURBO_"
  id.start <- 3000001
  id.end <- id.start + (robot.row.count - 1)
  term.suffixes <- id.start:id.end
  term.id <- paste0(id.prefix , term.suffixes)
  
  
  #### TERM LABEL
  # TODO NEW METHODS FOR PK SPECIFICTY
  #### use tl.augmenter approach as much as possible (above?)
  term.label <-
    paste0(
      prop.lab.component,
      ' for ',
      abtype,
      pk.pastable,
      pre.ready.for.robot$PartName,
      ' in ',
      sys.disp.name,
      " specimen ",
      challengetype
    )
  
  label.table()
  print(length(term.label))
  
  term.label[pre.ready.for.robot$TIME == "LP6924-7"] <-
    paste0(term.label[pre.ready.for.robot$TIME == "LP6924-7"], " with 24-hour sample collection")
  
  label.table()
  
  term.label <-
    paste0(term.label, div.pastable)
  
  label.table()
  
  term.label <- paste0(term.label, meth.pastable)
  
  label.table()
  
  term.label <- paste0(term.label, ' [', prop.quote.stripped, ']')
  
  label.table()
  print(length(term.label))
  
  ####
  
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
  
  ready.for.robot <- ready.for.robot[(!(single.quote.prob)) , ]
  
  na.prob <- grepl(pattern = "'NA'", x = ready.for.robot$analyte)
  table(na.prob, useNA = 'always')
  na.prob.rows  <-
    ready.for.robot[na.prob,]
  ready.for.robot <-
    ready.for.robot[(!(na.prob)),]
  
  unique.labels.only <- make.table.frame(ready.for.robot$term.label)
  unique.labels.only <-
    unique.labels.only$value[unique.labels.only$count == 1]
  shared.labels <-
    ready.for.robot[(!(ready.for.robot$term.label %in% unique.labels.only)),]
  shared.label.loinc.nums <-
    sub(
      pattern = "loinc:",
      replacement = "",
      x = shared.labels$database.cross.reference.IRI
    )
  shared.label.vs.ehr <-
    ehr.with.loinc.parts[ehr.with.loinc.parts$LOINC %in% shared.label.loinc.nums ,]
  ready.for.robot <-
    ready.for.robot[ready.for.robot$term.label %in% unique.labels.only,]
  
  write.table(
    ready.for.robot,
    file = 'lobo_template_headerless.csv',
    row.names = FALSE,
    col.names = FALSE,
    sep = ','
  )
  
  ###   ###   ###
  
  # some "failures"
  
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
  # methsuximide
  # normethsuximide
  
  # now run execute_template.sh
  