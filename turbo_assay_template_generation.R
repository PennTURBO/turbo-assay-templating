library(devtools)

# requires a properly formatted "turbo_R_setup.yaml" in the home directory of the user who started this script
# see https://gist.github.com/turbomam/a3915d00ee55d07510493a9944f96696 for template
source_gist(id = "https://gist.github.com/turbomam/f082295aafb95e71d109d15ca4535e46",
            sha1 = "dbc656aaf63b23dfdd35d875f6772e7c468170a4",
            filename = "turbo_R_setup.R")

library(tidyverse)

source("turbo_assay_template_generation_functions.R")

# clarify requirement for LOINC folder in location XXX

# # monitor
# ehr.bad.loincs
# ehr.na.loincs

# deal with multi-gene assays
# gene.check

# assays that are dropped because they include a numerator or a normalization

# ehr results or loincs with too few unique EMPIs
# we're not currently counting unique EMPIS per LOINC

# parts that are silently (?) dropped or less-annotated because they have multiple paths to root

####

# looks like some of the mapping work isn't going in a usefully sorted order
# could do cosine or jaccard between query and ols label, synonyms and truncatees all concatenated together

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
min.empis <- 2

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
    load(ehr_reference_labs_result_fn)
    # save filename in a config file
    # there's no SQL statements to regenrrate this in the script yet
    ehr_loinc_utilizastion <- read_csv("data/ehr_loinc_utilizastion.csv")
  }
  
  temp <- make.table.frame(ehr_loinc_utilizastion$RESULT_COUNT)
  temp$value <- as.numeric(temp$value)
  temp$count <- as.numeric(temp$count)
  
  ####
  
  # some of these may have new/replacement codes
  ehr.bad.loincs <-
    ehr.reference.labs.result[(!(is.na(ehr.reference.labs.result$LOINC))) &
                                (!(
                                  ehr.reference.labs.result$LOINC %in% LoincPartLink$LoincNumber
                                )) , ]
  
  ehr.na.loincs <-
    ehr.reference.labs.result[is.na(ehr.reference.labs.result$LOINC) , ]
  
  # ####
  
  lpl.min <-
    LoincPartLink[, c("PartTypeName", "LinkTypeName", "Property")]
  
  lpl.wide <-
    lpl.min %>% pivot_wider(names_from = PartTypeName,
                            values_from = LinkTypeName,
                            values_fn = length)
  
  lpl.long <-
    lpl.wide %>% pivot_longer(-Property, names_to = "PartTypeName", values_to = "count")
  lpl.long <- lpl.long[complete.cases(lpl.long),]
  
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
  lpl.long <- lpl.long[flag.flag, ]
  
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
  # some columns may be redundant... check jaccards
  # do names and codes???
  
  gene.check <-
    unique(LoincPartLink[LoincPartLink$Property == "http://loinc.org/property/analyte-gene" , c("LoincNumber", "PartName")])
  gene.check <- make.table.frame(gene.check$LoincNumber)
  gene.check <- gene.check$value[gene.check$count > 1]
  
  # (!(LoincPartLink %in% gene.check)) &
  
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
  
  codes.jaccard.calcs <- do.jaccard.calcs(max.one.gene.wide, 2)
  
  ehr.reference.labs.result.min.empi.count <-
    ehr.reference.labs.result[(!(is.na(ehr.reference.labs.result$LOINC))), c("LOINC", "UNIQUE_EMPI_COUNT")]
  ehr.reference.labs.result.min.empi.count <-
    aggregate(
      ehr.reference.labs.result.min.empi.count$UNIQUE_EMPI_COUNT,
      by = list(ehr.reference.labs.result.min.empi.count$LOINC),
      FUN = max
    )
  
  colnames(ehr.reference.labs.result.min.empi.count) <-
    c("LOINC", "UNIQUE_EMPI_COUNT")
  
  ehr.reference.labs.result.min.empi.count <-
    ehr.reference.labs.result.min.empi.count[ehr.reference.labs.result.min.empi.count$UNIQUE_EMPI_COUNT >= min.empis ,]
  
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
  
  max.one.gene.wide <-
    inner_join(x = max.one.gene.wide,
               y = ehr.reference.labs.result.min.empi.count,
               by = c("LoincNumber" = "LOINC"))
  
  dot.base.cols <-
    grepl(pattern = ".base$", x = colnames(max.one.gene.wide))
  not.dot.base.cols <-
    colnames(max.one.gene.wide)[(!(dot.base.cols))]
  max.one.gene.wide <- max.one.gene.wide[, not.dot.base.cols]
  
  common.constrained.assays <-
    get.common.constrained.assays()
  
  ####
  
  
  ont.tf <- tempfile()
  download.file(turbo.url, ont.tf)
  
  turbo.model <-
    rrdf::load.rdf(filename = ont.tf, format = turbo.format)
  
  Sys.time()
  system.time(turbo.labels <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.label.query))
  
  turbo.labels <- as.data.frame(turbo.labels)
  
  check.turbo.labels <- make.table.frame(turbo.labels$s)
  check.turbo.labels <-
    check.turbo.labels$value[check.turbo.labels$count == 1]
  
  turbo.labels <-
    turbo.labels[turbo.labels$s %in% check.turbo.labels ,]
  
  ####
  
  common.constrained.assays <-
    get.common.constrained.assays(component.code.list = loinc_to_obo_mapping_reviewed$PartNumber)
  
  analytes.from.reviewed.mapping <-
    intersect(loinc_to_obo_mapping_reviewed$PartNumber,
              Part$PartNumber[Part$PartTypeName == "COMPONENT"])
  
  ####    ####    ####    ####
  
  component.analyte.core.any.map.meth <-
    staged.loinc.mapping(
      needs.mapping = sort(unique(max.one.gene.wide$analyte.core)),
      already.reviewed = analytes.from.reviewed.mapping,
      ontology.rankings = ontology.rankings.min,
      current.code.col.name = "analyte.core",
      current.name.col.name = "PartName.analyte.core"
    )
  
####    ####    ####    ####
  
  temp <- component.analyte.core.any.map.meth$loinc.provided.secondary
  
  # bicarbonate
  
  # only single mappers are returned
  temp <- component.analyte.core.any.map.meth$paths.from.leaves
  temp <- make.table.frame(temp$CODE)
  
  
mappings.with.hints <- component.analyte.core.any.map.meth$any.method
  
mappings.with.hints$plus <-
  grepl(pattern = "+",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$period <-
  grepl(pattern = ".",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$gene <-
  grepl(pattern = "gene",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$cell <-
  grepl(pattern = "cell",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$free <-
  grepl(pattern = "free",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$total <-
  grepl(pattern = "total",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints$bound <-
  grepl(pattern = "bound",
        x = mappings.with.hints$PartDisplayName,
        fixed = TRUE)
mappings.with.hints <-
  unique(mappings.with.hints[, c(
    "LOINC_PartNumber",
    "PartDisplayName",
    "cell",
    "gene",
    "period",
    "plus",
    "free",
    "total",
    "bound",
    "target.term",
    "target.label",
    "onto.abbrev",
    "path.to.LOINC.root",
    "mapping.method",
    "rank",
    "q.vs.label.cosine",
    "other.cosine",
    "best.cosine",
    "best.cosine.scaled",
    "priority",
    "target.domain"
  )])

write_excel_csv(mappings.with.hints, needs.review.fn)

####

  # component (core analyte + ?)
  # method (no hierarchy, few good mappings)
  # property
  
  # scale
  # do Qn, (Ord & Nom both "categorical"?... need to create instances)
  # > sort(table(max.one.gene.wide$PartDisplayName.SCALE_TYP))
  # - OrdQn   Doc   Nar   Nom   Ord    Qn 
  # 17    21    25    60   251   721  2313 
  
  # system: GOOD!
  
  # time: make axioms but don't bother searching?
  # do PT, 24 (one EHR-relevant assay each for 72, 5, 48, 12)
  # > sort(table(max.one.gene.wide$PartDisplayName.TIME_ASPCT))
  # 12 hours             48 hours              5 hours             72 hours                    *          Unspecified 
  # 1                    1                    1                    1                    2                   12 
  # 24 hours Point in time (spot) 
  # 120                 3270                 12 
  
  # map
  # constrain
  # add more component subparts
  # no hierarchy found for systems or methods
  # very little useful method mapping
  
  # property
  # scale... using Qn, Ord and Nom only now
  # map property and scale to (output) datum types?
  
  
  # time... using Pt and 24 only right now
  
  # units.brainstorm <-
  #   as.data.frame.matrix(
  #     table(
  #       common.constrained.assays$PartName.PROPERTY,
  #       common.constrained.assays$PartName.SCALE_TYP
  #     )
  #   )
  # units.brainstorm$PROPERTY <- rownames(units.brainstorm)
  # units.brainstorm <-
  #   pivot_longer(data = units.brainstorm,
  #                -PROPERTY,
  #                names_to = "SCALE" ,
  #                values_to = "count")
  # units.brainstorm <- units.brainstorm[units.brainstorm$count > 0 , ]
  # 
  # prop.sum <-
  #   aggregate(units.brainstorm$count,
  #             by = list(units.brainstorm$PROPERTY),
  #             FUN = sum)
  # 
  # names(prop.sum) <- c("PROPERTY", "count.sum")
  # 
  # prop.sum$prop.lc <- tolower(prop.sum$PROPERTY)
  # 
  # absent_needed_props <- read_excel("absent_needed_props.xlsx",
  #                                   col_names = FALSE)
  # colnames(absent_needed_props) <- c("absent_needed_props")
  # absent_needed_props <-
  #   pull(absent_needed_props, absent_needed_props)
  # 
  # prop.sum <- prop.sum[prop.sum$prop.lc %in% absent_needed_props ,]
  # 
  # turbo_loinc_properties_reviewed <- read_csv("turbo_loinc_properties_reviewed.csv")
  # 
  # understood.properties <-
  #   turbo_loinc_properties_reviewed$PartName[(!(is.na(turbo_loinc_properties_reviewed$use))) &
  #                                              turbo_loinc_properties_reviewed$use == "y"]
  
  ####
  
  # # # saves columns, not tied together in any way
  # prop.sum.json <- jsonlite::toJSON(prop.sum)
  # attr(prop.sum.json, "tag") <- "tagattr"
  # yaml::write_yaml(prop.sum.json, "properties_followup.yaml")
  
  # goal here: get the reviewed mappings, the LOINC provided mappings, 
  #   the (highly scrutinized) bp mappings and the minimially scrutinized ols search results 
  #   into one table that's easy to review
  
  # need better column ordering
  # query (loinc) ID, query (loinc) name, OBO:ID, OBO name, quality metrics...
  
  # these don't show ow important the parts are in terms of number of loinc assays,  
  #   number of EHR assays, or (minimum?) number of patients
  
  ####
  
  # column.usefulness <- get.column.usefulness(common.constrained.assays)
  # write.csv(column.usefulness, "column_usefulness.csv", row.names = FALSE)

####

cd_markers <- read_csv(cd_markersmapping_fn)

temp <-
  apply(
    X = cd_markers,
    MARGIN = 1,
    FUN = function(current_marker) {
      print(current_marker[["PartNumber"]])
      print(current_marker[["PartDisplayName"]])
      
      inner_temp <- ols.serch.term.labels.universal(
        current.string = current_marker[["PartDisplayName"]],
        current.id = current_marker[["PartNumber"]],
        strip.final.s = FALSE,
        ontology.filter = paste0("&ontology=", onto.list.for.ols),
        kept.row.count = 30,
        req.exact = 'false'
      )
      return(inner_temp)
      
    }
  )

temp <- do.call(rbind.data.frame, temp)


temp <-
  temp[, c(
    "loinc.part",
    "query",
    "iri",
    "label",
    "rank",
    "short_form",
    "obo_id",
    "ontology_name",
    "ontology_prefix"
  )]

write.csv(temp, file = "cd_markers_needs_review.csv", row.names = FALSE)

####

loinc_combo_analytes <- read_csv(loinc_combo_analytes_fn)
loinc_combo_analytes <-
  left_join(loinc_combo_analytes, turbo.labels, by = c("obo1" = "s"))
loinc_combo_analytes <-
  left_join(loinc_combo_analytes, turbo.labels, by = c("obo2" = "s"))

loinc_combo_analytes$axiom.labels <-
  paste0("'(",
         loinc_combo_analytes$o.x,
         "' and '",
         loinc_combo_analytes$o.y,
         "')")

write.csv(loinc_combo_analytes, "loinc_combo_analytes_new.csv", row.names = FALSE)

#### many "cells by surface marker" and "combo analytes" are already integrated into 
####   loinc_to_obo_mapping_reviewed.csv

