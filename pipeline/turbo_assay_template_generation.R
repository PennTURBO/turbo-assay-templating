library(devtools)

# todo: move this gist into the PennTURBO organization

# requires a properly formatted "turbo_R_setup.yaml" in the home directory of the user who started this script
# see https://gist.github.com/turbomam/a3915d00ee55d07510493a9944f96696 for template
source_gist(id = "https://gist.github.com/turbomam/f082295aafb95e71d109d15ca4535e46",
            sha1 = "dbc656aaf63b23dfdd35d875f6772e7c468170a4",
            filename = "turbo_R_setup.R")

# library(tidyverse)
library(tidyr)
library(tibble)

source("turbo_assay_template_generation_functions.R")

# clarify requirement for Loinc_2 folder in location pipeline/data

####

# keep monitoring these for "lost assays"
#   ehr.bad.loincs
#   ehr.na.loincs

#   deal with multi-gene assays
#   see gene.check

#   assays that are dropped because they include a numerator or a normalization

#   LOINC codes whose usage in the EHR is less than the threshold
#     we're not currently counting unique EMPIS per LOINC

#   parts that are silently (?) dropped or less-annotated because they have multiple paths to root?

####

# make sure remote mapping requests are being submitted in a usefully sorted order

# for picking best OLS search result:
#   could calculate  cosine or jaccard between query and ols label, synonyms and truncatees, all concatenated together

# expect 'reviewed part frame' from staged.loinc.mapping to be empty unless passed as an argument
# but start collecting mappings for all parts into one file
# bp mappings and ols search results should have
#  loinc code and label
#  external (obo?) id (IRI?) and label
#  hierarchy if available
#  mapping quality metadata

# staged.loinc.mapping variable shouldn't contain the word 'system'... they are part agnostic

# > library(tidyverse)
# ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
# ✓ tibble  3.0.1     ✓ purrr   0.3.4
# ✓ tidyr   1.1.0     ✓ forcats 0.5.0
# ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# x NLP::annotate()         masks ggplot2::annotate()
# x tibble::as_data_frame() masks igraph::as_data_frame(), dplyr::as_data_frame()
# x randomForest::combine() masks dplyr::combine()
# x purrr::compose()        masks igraph::compose()
# x NLP::content()          masks httr::content()
# x tidyr::crossing()       masks igraph::crossing()
# x dplyr::filter()         masks stats::filter()
# x purrr::flatten()        masks jsonlite::flatten()
# x igraph::groups()        masks dplyr::groups()
# x dplyr::lag()            masks stats::lag()
# x purrr::lift()           masks caret::lift()
# x randomForest::margin()  masks ggplot2::margin()
# x caret::progress()       masks httr::progress()
# x purrr::simplify()       masks igraph::simplify()

# setwd("~/GitHub/turbo-assay-templating/")

# VPN and tunnel may be required
# set that up outside of this script

repeat.ehr.loinc.query <- FALSE

min.results <- 20

id.prefix <- "http://purl.obolibrary.org/obo/OBI_"
id.start <- 3000001

Part.min.cols <-
  c("PartNumber", "PartName", "PartDisplayName", "Status")

excluded.Properties <- c("http://loinc.org/property/category")

bioontologies.base <- "http://data.bioontology.org/ontologies/"

obo.base <- "http://purl.obolibrary.org/obo/"

bp.loinc.prefix <- "http://purl.bioontology.org/ontology/LNC/"

chebi.sys.iri <- "https://www.ebi.ac.uk/chebi"

ncbitaxon.sys.iri <- "https://www.ncbi.nlm.nih.gov/taxonomy"

turbo.url <-
  "https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/robot_implementation/turbo_merged.ttl"

turbo.format <- "TURTLE"

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
  
  turbo.molent.query <- "
  select
  *
  where
  {?s <http://www.w3.org/2000/01/rdf-schema#subClassOf>* <http://purl.obolibrary.org/obo/CHEBI_23367> ;
    <http://www.w3.org/2000/01/rdf-schema#label> ?l .
  }
  "
  
  
  turbo.atom.query <- "
  select
  *
  where
  {?s <http://www.w3.org/2000/01/rdf-schema#subClassOf>* <http://purl.obolibrary.org/obo/CHEBI_33250> ;
    <http://www.w3.org/2000/01/rdf-schema#label> ?l .
  }
  "
  
  hor.path <- "data/"
  
  needs.review.fn <-  "core_analyte_mapping_needs_review.csv"
  
  cd_markersmapping_fn <- "data/cd_markers.csv"
  
  loinc_combo_analytes_fn <- "data/loinc_combo_analytes.csv"
  
  ehr_reference_labs_result_fn <-
    "data/ehr_reference_labs_result.Rdata"
  
  loinc_to_obo_mapping_reviewed <-
    read_csv("data/loinc_to_obo_mapping_reviewed.csv")
  
  ####
  
  MultiAxialHierarchy <-
    read_csv("data/Loinc_2/AccessoryFiles/MultiAxialHierarchy/MultiAxialHierarchy.csv")
  LoincPartLink_Primary <-
    read_csv("data/Loinc_2/AccessoryFiles/PartFile/LoincPartLink_Primary.csv")
  LoincPartLink_Supplementary <-
    read_csv("data/Loinc_2/AccessoryFiles/PartFile/LoincPartLink_Supplementary.csv")
  LoincPartLink <-
    rbind.data.frame(LoincPartLink_Primary, LoincPartLink_Supplementary)
  Part <- read_csv("data/Loinc_2/AccessoryFiles/PartFile/Part.csv")
  Loinc <- read.csv("data/Loinc_2/LoincTable/Loinc.csv")
  PartRelatedCodeMapping <-
    read_csv("data/Loinc_2/AccessoryFiles/PartFile/PartRelatedCodeMapping.csv")
  
  Part.min <- Part[, Part.min.cols]
  
  harmonized_ontology_rankings <-
    read.csv(paste0(hor.path, config$harmonized_ontology_rankings),
             stringsAsFactors = FALSE)
  
  # say rank, not priority
  ontology.rankings.min <-
    cbind.data.frame(
      ontology.iri = paste(bioontologies.base,
                           harmonized_ontology_rankings$rhs,
                           sep = ""),
      ontology.rank  = harmonized_ontology_rankings$priority
    )
  
  onto.list.for.ols <-
    paste0(sort(tolower(harmonized_ontology_rankings$rhs)), collapse = ',')
  
  ehr.driver <-
    JDBC(driverClass = "oracle.jdbc.OracleDriver",
         classPath = config$oracle.jdbc.path)
  
  ehr.con.string <- paste0(
    "jdbc:oracle:thin:@//",
    config$pds.host,
    ":",
    config$pds.port,
    "/",
    config$pds.database
  )
  
  ehr.connection <- NULL
  
  ####
  
  if (repeat.ehr.loinc.query) {
    ehr.test.and.refresh()
    
    ehr.reference.labs.query <- "
SELECT
	rlri.PK_LAB_RESULT_ITEM_ID ,
	RLRI.LAB_RESULT_ITEM_CODE ,
	RLRI.LAB_RESULT_ITEM_DESCRIPTION ,
	RLRI.LOINC ,
	count(DISTINCT pe.EMPI ) AS UNIQUE_EMPI_COUNT
FROM
	MDM.R_LAB_RESULT_ITEM rlri
JOIN MDM.LAB_RESULTS lr ON
	rlri.PK_LAB_RESULT_ITEM_ID = lr.FK_LAB_RESULT_ITEM_ID
JOIN mdm.PATIENT_ENCOUNTER pe ON
	lr.FK_PATIENT_ENCOUNTER_ID = pe.PK_PATIENT_ENCOUNTER_ID
GROUP BY
	rlri.PK_LAB_RESULT_ITEM_ID ,
	RLRI.LAB_RESULT_ITEM_CODE ,
	RLRI.LAB_RESULT_ITEM_DESCRIPTION ,
	RLRI.LOINC"
    
    print(Sys.time())
    timed.system <- system.time(
      ehr.reference.labs.result <-
        dbGetQuery(conn = ehr.connection, statement = ehr.reference.labs.query)
    )
    print(Sys.time())
    print(timed.system)
    
    # "
    # SELECT
    # 	lr1.PK_LAB_RESULT_ID
    # FROM
    # 	mdm.LAB_RESULTS lr1 MINUS
    # SELECT
    # 	lr2.PK_LAB_RESULT_ID
    # FROM
    # 	mdm.LAB_RESULTS lr2,
    # 	mdm.R_LAB_RESULT_ITEM rlri2
    # WHERE
    # 	lr2.FK_LAB_RESULT_ITEM_ID = rlri2.PK_LAB_RESULT_ITEM_ID"
    # # takes 15 minutes
    # # so far, have not any observed lab results that don't map to a R lab result item
    
    # # Close connection
    # dbDisconnect(ehr.connection)
    save(ehr.reference.labs.result, file = ehr_reference_labs_result_fn)
  } else {
    # # deprecated?
    # load(ehr_reference_labs_result_fn)
    
    # save filename in a config file
    # there's no SQL statements to regenerate this in the script yet
    ehr_loinc_utilizastion <-
      read_csv("data/ehr_loinc_utilizastion.csv")
  }
  
  loinc.utilization <-
    make.table.frame(ehr_loinc_utilizastion$RESULT_COUNT)
  names(loinc.utilization) <-
    c("loinc.utilization", "number.of.loincs")
  loinc.utilization$loinc.utilization <-
    as.numeric(loinc.utilization$loinc.utilization)
  loinc.utilization$number.of.loincs <-
    as.numeric(loinc.utilization$number.of.loincs)
  
  ####
  
  # ehr.na.loincs <-
  #   ehr.reference.labs.result[is.na(ehr.reference.labs.result$LOINC) ,]
  
  # some of these may have new/replacement codes
  ehr.bad.loincs <-
    ehr_loinc_utilizastion[(!(is.na(ehr_loinc_utilizastion$LOINC))) &
                             (!(
                               ehr_loinc_utilizastion$LOINC %in% LoincPartLink$LoincNumber
<<<<<<< HEAD
                             )) ,]
=======
                             )) , ]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  # ####
  
  lpl.min <-
    LoincPartLink[, c("PartTypeName", "LinkTypeName", "Property")]
  
  lpl.wide <-
    lpl.min %>% pivot_wider(names_from = PartTypeName,
                            values_from = LinkTypeName,
                            values_fn = length)
  
  lpl.long <-
    lpl.wide %>% pivot_longer(-Property, names_to = "PartTypeName", values_to = "count")
<<<<<<< HEAD
  lpl.long <- lpl.long[complete.cases(lpl.long),]
=======
  lpl.long <- lpl.long[complete.cases(lpl.long), ]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  lpl.long$search.flag <-
    lpl.long$Property == "http://loinc.org/property/search"
  lpl.long$doc.flag <-
    grepl(pattern = "http://loinc.org/property/document-",
          x = lpl.long$Property,
          fixed = TRUE)
  lpl.long$rad.flag <-
    grepl(pattern = "http://loinc.org/property/rad-",
          x = lpl.long$Property,
          fixed = TRUE)
  
  flag.flag <- grepl(pattern = "flag", x = colnames(lpl.long))
  flag.flag <- lpl.long[, flag.flag]
  flag.flag[] <- lapply(flag.flag[], as.numeric)
  flag.flag$rowSums <- rowSums(flag.flag)
  flag.flag <- flag.flag$rowSums == 0
<<<<<<< HEAD
  lpl.long <- lpl.long[flag.flag, ]
=======
  lpl.long <- lpl.long[flag.flag,]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  flag.flag <- grepl(pattern = "flag", x = colnames(lpl.long))
  flag.flag <- colnames(lpl.long)[flag.flag]
  
  lpl.long <- lpl.long[, setdiff(colnames(lpl.long), flag.flag)]
  
  lpl.wide <-
    lpl.long %>% pivot_wider(names_from = PartTypeName ,
                             values_from = count,
                             values_fn = length)
  
  lpl.wide$rowSums <-
    rowSums(lpl.wide[, 2:ncol(lpl.wide)], na.rm = TRUE)
  print(sort(lpl.wide$rowSums))
  # no Property is associated with more than one PartTypeName
  
  lpl.wide.colSums <-
    colSums(lpl.wide[, 2:ncol(lpl.wide)], na.rm = TRUE)
  print(sort(lpl.wide.colSums))
  # some PartTypeNames are associated with more than one Property
  
  # skip category... almost guaranteed many to many
  # some columns may be redundant... check jaccards?
  
  gene.check <-
    unique(LoincPartLink[LoincPartLink$Property == "http://loinc.org/property/analyte-gene" , c("LoincNumber", "PartName")])
  gene.check <- make.table.frame(gene.check$LoincNumber)
  gene.check <- gene.check$value[gene.check$count > 1]
  
  accepted.assays <- setdiff(LoincPartLink$LoincNumber, gene.check)
  accepted.properties <-
    setdiff(lpl.wide$Property, excluded.Properties)
  lpl <-
    unique(LoincPartLink[LoincPartLink$Property %in% accepted.properties &
                           LoincPartLink$LoincNumber %in% accepted.assays,
                         c("LoincNumber", "PartName", "Property")])
  
  max.one.gene.wide <-
    lpl.wide.builder(names.from = 'Property',
                     values.from = 'PartNumber',
                     to.strip = "http://loinc.org/property/")
  
  # not really using this any more
  # had been concerned that some columns were very similar but not identical, so not suitabel for deleting
  # like TIME and time.core
  codes.jaccard.calcs <- do.jaccard.calcs(max.one.gene.wide, 2)
  
  Loinc.min <-
    Loinc[Loinc$LOINC_NUM %in% max.one.gene.wide$LoincNumber, c(
      "LOINC_NUM",
      "VersionLastChanged",
      "CHNG_TYPE",
      "DefinitionDescription",
      "STATUS",
      "CONSUMER_NAME",
      "CLASSTYPE",
      "FORMULA",
      "UNITSREQUIRED",
      "SUBMITTED_UNITS",
      "ORDER_OBS",
      "CDISC_COMMON_TESTS",
      "HL7_FIELD_SUBFIELD_ID",
      "EXTERNAL_COPYRIGHT_NOTICE",
      "EXAMPLE_UNITS",
      "LONG_COMMON_NAME",
      "UnitsAndRange",
      "EXAMPLE_UCUM_UNITS",
      "EXAMPLE_SI_UCUM_UNITS",
      "STATUS_REASON",
      "STATUS_TEXT",
      "CHANGE_REASON_PUBLIC",
      "COMMON_TEST_RANK",
      "COMMON_ORDER_RANK",
      "COMMON_SI_TEST_RANK",
      "HL7_ATTACHMENT_STRUCTURE",
      "EXTERNAL_COPYRIGHT_LINK",
      "PanelType",
      "AskAtOrderEntry",
      "AssociatedObservations",
      "VersionFirstReleased",
      "ValidHL7AttachmentRequest"
    )]
  
  ####
  
  part.names.sorted <-
    sort(colnames(max.one.gene.wide[2:ncol(max.one.gene.wide)]))
  Part.min.cols.merge.prep <- setdiff(Part.min.cols, "PartNumber")
  placeholder <-
    lapply(Part.min.cols.merge.prep, function(current_col) {
      print(current_col)
      max.one.gene.wide <<-
        max.one.gene.wide %>% add_column(!!(current_col) := NA)
    })
  
  for (i in 1:length(part.names.sorted)) {
    joinee = "base"
    if (i > 1) {
      joinee <- part.names.sorted[i - 1]
    }
    joiner <- part.names.sorted[i]
    by.trick <- "PartNumber"
    names(by.trick) <- joiner
    print(by.trick)
    max.one.gene.wide <<-
      left_join(
        x = max.one.gene.wide,
        y = Part.min,
        by = by.trick,
        all.x = TRUE,
        suffix = c(
          paste0(".", joinee, collapse = ""),
          paste0(".", joiner, collapse = "")
        )
      )
  }
  
  max.one.gene.wide <-
    left_join(x = max.one.gene.wide,
              y = Loinc.min,
              by = c("LoincNumber" = "LOINC_NUM"))
  
  ####    ####    ####    ####
  
  commonly.ordered.loinc.codes <-
<<<<<<< HEAD
    ehr_loinc_utilizastion[ehr_loinc_utilizastion$RESULT_COUNT > min.results ,]
=======
    ehr_loinc_utilizastion[ehr_loinc_utilizastion$RESULT_COUNT > min.results , ]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  max.one.gene.wide <-
    inner_join(x = max.one.gene.wide,
               y = commonly.ordered.loinc.codes,
               by = c("LoincNumber" = "LOINC"))
  
  
  dot.base.cols <-
    grepl(pattern = ".base$", x = colnames(max.one.gene.wide))
  not.dot.base.cols <-
    colnames(max.one.gene.wide)[(!(dot.base.cols))]
  max.one.gene.wide <- max.one.gene.wide[, not.dot.base.cols]
  
  # may not really need this any more... similar constrains are applied by the merges below
  common.constrained.assays <-
    get.common.constrained.assays()
  
  ####
  
  # TURBO labels aren't really being used for anything any more
  # they had been handy for checking whether the OBO terms necessary for the assays had been imported into the TURBO ontology
  # but many of the OBO targets of LOINC term mappings are complex class expression now
  
  # now leaning towards checking for eligible analytes externally to this R script
  #  SAPRQL query via ROBOT?
  
  ont.tf <- tempfile()
  download.file(turbo.url, ont.tf)
  
  turbo.model <-
    rrdf::load.rdf(filename = ont.tf, format = turbo.format)
  
  Sys.time()
  system.time(turbo.labels <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.label.query))
  
  
  Sys.time()
  system.time(turbo.molents <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.molent.query))
  
  Sys.time()
  system.time(turbo.atoms <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.atom.query))
  
  # what about double quotation?
  # what about combo analytes?
  # turbo.analyte.eligible <- sort(union(turbo.atoms[,2], turbo.molents[,2]))
  # make them all evaluants and promote with reasoning or explicit SPARQL
  
  turbo.labels <- as.data.frame(turbo.labels)
  
  check.turbo.labels <- make.table.frame(turbo.labels$s)
  check.turbo.labels <-
    check.turbo.labels$value[check.turbo.labels$count == 1]
  
  turbo.labels <-
<<<<<<< HEAD
    turbo.labels[turbo.labels$s %in% check.turbo.labels ,]
=======
    turbo.labels[turbo.labels$s %in% check.turbo.labels , ]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  ####
  
  common.constrained.assays <-
    get.common.constrained.assays(component.code.list = loinc_to_obo_mapping_reviewed$PartNumber)
  
  mergable <-
<<<<<<< HEAD
    loinc_to_obo_mapping_reviewed[loinc_to_obo_mapping_reviewed$pure.obo &
                                    is.na(loinc_to_obo_mapping_reviewed$problematic), c("PartNumber",
                                                                                        "label.axiom", "text.for.label")]
=======
    loinc_to_obo_mapping_reviewed[loinc_to_obo_mapping_reviewed$pure.obo, c("PartNumber",
                                                                            "label.axiom", "text.for.label")]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  ####
  
  next.merge <- max.one.gene.wide
  
  # max.one.gene.wide constrained only by minimum result count
  # common.constrained.assays
  #  already constrained (but not merged) by
  #  component mapped to OBO label based axioms
  #  fairly loose default constraints on component, method.code.list, property, scale, system, time
  #  "min.min.patients" parameter name may be confusing now that we're tracking LOINC code utilization, not patient count
  
  # throwaway ensuring that next merges get helpful suffixes
  next.merge <-
    inner_join(
      x = next.merge,
      y = mergable,
      by = c("analyte.core" = "PartNumber")
    )
  
  # for component
  next.merge <-
    inner_join(
      x = next.merge,
      y = mergable,
      by = c("analyte.core" = "PartNumber"),
      suffix = c(".placeholder", ".analyte.core")
    )
  
  next.merge <-
    inner_join(
      x = next.merge,
      y = mergable,
      by = c("PROPERTY" = "PartNumber"),
      suffix = c(".analyte.core", ".PROPERTY")
    )
  
<<<<<<< HEAD
  next.merge <- next.merge[next.merge$TIME_ASPCT == "LP6960-1", ]
=======
  next.merge <- next.merge[next.merge$TIME_ASPCT == "LP6960-1",]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  next.merge <-
    inner_join(
      x = next.merge,
      y = mergable,
      by = c("SYSTEM" = "PartNumber"),
      suffix = c(".PROPERTY", ".SYSTEM")
    )
  
<<<<<<< HEAD
  next.merge <- next.merge[next.merge$SCALE_TYP == "LP7753-9", ]
  
  # leave as is
  next.merge <- next.merge[is.na(next.merge$METHOD_TYP), ]
  next.merge <- next.merge[is.na(next.merge$analyte.divisor), ]
  next.merge <- next.merge[is.na(next.merge$challenge), ]
  next.merge <- next.merge[is.na(next.merge$analyte.gene), ]
  
  # get these back in there by adding to OBI
  next.merge <- next.merge[is.na(next.merge$analyte.numerator), ]
  next.merge <- next.merge[is.na(next.merge$adjustment), ]
  next.merge <- next.merge[is.na(next.merge$count), ]
  
  next.merge.nosuff <-
    next.merge[is.na(next.merge$analyte.suffix),]
=======
  next.merge <- next.merge[next.merge$SCALE_TYP == "LP7753-9",]
  
  # leave as is
  next.merge <- next.merge[is.na(next.merge$METHOD_TYP),]
  next.merge <- next.merge[is.na(next.merge$analyte.divisor),]
  next.merge <- next.merge[is.na(next.merge$challenge),]
  next.merge <- next.merge[is.na(next.merge$analyte.gene),]
  
  # get these back in there by adding to OBI
  next.merge <- next.merge[is.na(next.merge$analyte.numerator),]
  next.merge <- next.merge[is.na(next.merge$adjustment),]
  next.merge <- next.merge[is.na(next.merge$count),]
  
  next.merge.nosuff <-
    next.merge[is.na(next.merge$analyte.suffix), ]
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  next.merge.suffixes <-
    inner_join(
      x = next.merge,
      y = mergable,
      by = c("analyte.suffix" = "PartNumber"),
      suffix = c(".SYSTEM", ".analyte.suffix"),
      na_matches = "never"
    )
  
  next.merge <-
    plyr::rbind.fill(next.merge.suffixes, next.merge.nosuff)
  
  next.merge.row.count <- nrow(next.merge)
  all.blanks <- rep("", next.merge.row.count)
  
  id.end <- id.start + (next.merge.row.count - 1)
  term.suffixes <- id.start:id.end
  
  next.merge$robot.01.term.id <- paste0(id.prefix , term.suffixes)
  
  next.merge$robot.03.curation.status <- all.blanks
  
  # skip alternative term except for use by a private user, or with approval by loinc (in which case, use their display label)
  
  next.merge$robot.04.alt.term <- all.blanks
  
<<<<<<< HEAD
  next.merge$robot.06.def.source <-
    rep(
      "Inspired by LOINC. Templated by https://github.com/PennTURBO/turbo-assay-templating",
      next.merge.row.count
    )
=======
  next.merge$robot.06.def.source <- all.blanks
>>>>>>> 8ec5e99c1b2407cb25087eaf8e4474030e7ed6c5
  
  next.merge$robot.07.usage.example <- all.blanks
  
  next.merge$robot.08.ed.note <- all.blanks
  
  next.merge$robot.09.term.ed <-
    "Mark Andrew Miller|ORCID:0000-0001-9076-6066|Chris Stoeckert|ORCID:0000-0002-5714-991X"
  
  next.merge$robot.10.requestor <- all.blanks
  
  next.merge$robot.11.issue.tracker.id <-
    "https://github.com/obi-ontology/obi/issues/1153"
  
  next.merge$robot.12.logical.type <- "subclass"
  
  next.merge$robot.13.local.note <- all.blanks
  
  # 'assay of specimen from organism' is still a TURBO term
  next.merge$robot.14.parent.class <- 'assay'
  
  next.merge$robot.15.mat.proc <- all.blanks
  
  next.merge$robot.16.det.tech <- all.blanks
  
  next.merge$robot.17.evaluant <- next.merge$label.axiom.SYSTEM
  
  # analyte? only if it's a molecular entity
  # otherwise, target entity column
  
  next.merge$robot.18.analyte <- all.blanks
  
  next.merge$robot.19.device <- all.blanks
  
  next.merge$robot.20.reagent <- all.blanks
  
  next.merge$robot.21.mol.lab <- all.blanks
  
  next.merge$robot.22.input <- all.blanks
  
  next.merge$robot.23.output <- next.merge$label.axiom.PROPERTY
  
  next.merge$robot.24.target.ent <-
    next.merge$label.axiom.analyte.core
  
  next.merge$robot.25.obj <- all.blanks
  
  next.merge$robot.26.assoc.ax <- all.blanks
  
  # with approval, or for private use, LOINC codes will go here
  # next.merge$robot.27.dbxr <- all.blanks
  # http://purl.bioontology.org/ontology/LNC/4100-4
  # https://loinc.org/4100-4/
  # next.merge$LoincNumber
  next.merge$robot.27.dbxr <-
    paste0("http://purl.bioontology.org/ontology/LNC/",
           next.merge$LoincNumber)
  
  # there will different patterns for other suffixes, like RNA, DNA and Ag
  next.merge$robot.24.target.ent[!is.na(next.merge$analyte.suffix)] <-
    paste0(
      next.merge$label.axiom[!is.na(next.merge$analyte.suffix)],
      " and 'has disposition to bind' some ",
      next.merge$robot.24.target.ent[!is.na(next.merge$analyte.suffix)]
    )
  
  tidy.targent <- next.merge$text.for.label.analyte.core
  
  tidy.targent[!is.na(next.merge$analyte.suffix)] <-
    paste0(next.merge$text.for.label[!is.na(next.merge$analyte.suffix)],
           " against ",
           tidy.targent[!is.na(next.merge$analyte.suffix)])
  
  
  tidy.evaluant <- next.merge$text.for.label.SYSTEM
  
  tidy.units <- next.merge$text.for.label.PROPERTY
  
  next.merge$robot.02.label <-
    paste0("Assay for ",
           tidy.targent,
           " in ",
           tidy.evaluant,
           " [",
           tidy.units,
           "]")
  
  label.redundancy <- make.table.frame(next.merge$robot.02.label)
  next.merge$label.redundancy <- label.redundancy$count
  
  next.merge$robot.05.definition <-
    paste0(
      "An assay for ",
      tidy.targent,
      " in ",
      tidy.evaluant,
      " that generates a scalar datum with ",
      tidy.units,
      " units"
    )
  
  robot.cols <- grepl(pattern = "robot", x = colnames(next.merge))
  robot.cols <- sort(colnames(next.merge)[robot.cols])
  
  robot.frame <-
    next.merge[next.merge$label.redundancy == 1 , robot.cols]
  
  unofficial.dbxr <-
    next.merge[, c("LoincNumber", "robot.01.term.id")]
  
  write.table(
    robot.frame,
    file = 'build/turbo_assay_template_headerless.csv',
    row.names = FALSE,
    col.names = FALSE,
    sep = ','
  )
  
  write.table(
    unofficial.dbxr,
    file = 'build/unofficial_dbxr.csv',
    row.names = FALSE,
    col.names = FALSE,
    sep = ','
  )
  
  save.image('build/turbo_assays_latest_image.Rdata')
  
  head(next.merge[, c("LoincNumber", "robot.01.term.id")])
  