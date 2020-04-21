###   ###   ###

# the next several block are all just alternative houskeeping or discovery tricks

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

# PS how did I make the turbo_loinc_*_reviewed files above?
# write.csv(pds.prominent.loinc.parts[pds.priminent.loinc.parts$PartTypeName == "PROPERTY",], "turbo_loinc_properties.csv")
# etc.


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

###   ###   ###

# find common ancestors of unmapped LOINC components
# rdflib, rrdf, (SPARQL... against what endpoint), igraph (from LOINC tabular files)

# some more perspective/prioritization regarding
# what kind of patterns to look for

# why not start with all component parts or even all part terms
# as opposed to "salvage mode" below

MultiAxialHierarchy.relevant <-
  MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% pds.prominent.loinc.parts$PartTypeVal , ]

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
  node.appearances[node.appearances$appearances >= 2 , ]

common.nodes <-
  left_join(common.nodes, Part, by = c("part" = "PartNumber"))

common.nodes <- common.nodes[complete.cases(common.nodes),]

###   ###   ###

# previously looked at common super classes
# now look at common tokens

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
  tabulated.tokens[tabulated.tokens$count > 3 , ]

###   ###   ###

# # moot, now that I fixed core analyte column
# no.suffix <-
#   apply(
#     X = pds.with.loinc.parts,
#     MARGIN = 1,
#     FUN = function(current.row) {
#       # print(names(current.row))
#       # print(current.row[['PartName.suffix']])
#       # print(current.row[['analyte.vals']])
#       if ('analyte.vals' %in% names(current.row) &
#           'PartName.suffix' %in% names(current.row)) {
#         if (!is.na(current.row[['PartName.suffix']]) > 0 &
#             !is.na(current.row[['analyte.vals']]) > 0) {
#           inner <- sub(pattern = current.row[['PartName.suffix']],
#                        replacement = '',
#                        x = current.row[['analyte.vals']])
#           print(inner)
#           return(inner)
#         } else {
#           return(NA)
#         }
#       }
#     }
#   )
#
# pds.with.loinc.parts$no.suffix <- no.suffix

###   ###   ###

# first pass after text searching... but just use XXX instead

keepers <-
  c(
    'aeo',
    'aro',
    'bto',
    'caro',
    'chebi',
    'chmo',
    'cl',
    'clo',
    'cmo',
    'doid',
    'dron',
    'envo',
    'fma',
    'foodon',
    'go',
    'hp',
    'mi',
    'micro',
    'mmo',
    'mod',
    'mp',
    'ncbitaxon',
    'oae',
    'oba',
    'obi',
    'obib',
    'ogg',
    'ogms',
    'pato',
    'pr',
    'so',
    'uberon',
    'vo',
    'vto',
    'xco'
  )

###   ###   ###

names.and.numbers <-
  unique(source.relevant.links[, c("PartNumber", "PartName")])

common.core.analyte.numbers <- table(source.relevant.links$PartNumber)
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

###   ###   ###

common.divisors <-
  LoincPartLink[LoincPartLink$LoincNumber %in% ehr_loinc_counts$LOINC &
                  LoincPartLink$Property == 'http://loinc.org/property/analyte-divisor' ,]

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

###   ###   ###

# look for low-hanging fruit patterns
followup <-
  LoincPartLink[LoincPartLink$LoincNumber %in% pds.with.loinc.parts$LOINC , ]


