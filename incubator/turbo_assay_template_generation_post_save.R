
# 
# # this was the last merge, so it didn't get a suffix
# # should be prepred that the column name may change in the future
# robot.with.suffix.frame <-
#   next.merge[next.merge$label.redundancy == 1 &
#                !is.na(next.merge$analyte.suffix), c(robot.cols, "label.axiom")]
# 
# robot.with.suffix.frame$robot.02.label <-
#   paste0(
#     robot.with.suffix.frame$label.axiom,
#     " and 'has disposition to bind' some ",
#     robot.with.suffix.frame$robot.02.label
#   )

# reorder columns
# handle suffixes
# but save analyte/target entity for later

# robot.frame$robot.label <-
#   paste0(
#     "Quantitiative assay for ",
#     robot.frame$robot.target.ent,
#     " in ",
#     robot.frame$robot.evaluant, " [", robot.frame$robot.output,"]"
#   )
# 
# label.redundancy <- make.table.frame(robot.frame$robot.label)
# 
# robot.frame$robot.definition <-
#   paste0(
#     "An assay for ",
#     robot.frame$robot.target.ent,
#     " in ",
#     robot.frame$robot.evaluant,
#     " that generates a ",
#     robot.frame$robot.output,
#     " on a quantitative scale"
#   )


####    ####    ####    ####

analytes.from.reviewed.mapping <-
  intersect(loinc_to_obo_mapping_reviewed$PartNumber,
            Part$PartNumber[Part$PartTypeName == "COMPONENT"])

####    ####    ####    ####

if (FALSE) {
  component.analyte.core.any.map.meth <-
    staged.loinc.mapping(
      needs.mapping = sort(unique(max.one.gene.wide$analyte.core)),
      already.reviewed = analytes.from.reviewed.mapping,
      ontology.rankings = ontology.rankings.min,
      current.code.col.name = "analyte.core",
      current.name.col.name = "PartName.analyte.core"
    )
}


####    ####    ####    ####

method.utilizaston <-
  aggregate(
    max.one.gene.wide$RESULT_COUNT,
    by = list(max.one.gene.wide$METHOD_TYP),
    FUN = sum
  )

colnames(method.utilizaston) <-
  c("METHOD_TYP", "RESULT_COUNT.sum")

method.utilizaston <-
  inner_join(x = method.utilizaston,
             y = Part,
             by = c("METHOD_TYP" = "PartNumber"))

# write.csv(method.utilizaston, "ehr_common_loinc_methods.csv")

####    ####    ####    ####

challenge.utilizaston <-
  aggregate(
    max.one.gene.wide$RESULT_COUNT,
    by = list(max.one.gene.wide$challenge),
    FUN = sum
  )

colnames(challenge.utilizaston) <-
  c("challenge", "RESULT_COUNT.sum")

challenge.utilizaston <-
  inner_join(x = challenge.utilizaston,
             y = Part,
             by = c("challenge" = "PartNumber"))

# write.csv(challenge.utilizaston, "ehr_common_loinc_challenges.csv")

####



####    ####    ####    ####

# for.meth.chal.counting <-
#   inner_join(x = max.one.gene.wide,
#              y = loinc_to_obo_mapping_reviewed,
#              by = c("analyte.core" = "PartNumber"))
#
# ehr.common.loinc.methods <- make.table.frame(for.meth.chal.counting$METHOD_TYP)
# colnames(ehr.common.loinc.methods) <- c("PartNumber","count")
# ehr.common.loinc.methods <-
#   left_join(x = ehr.common.loinc.methods,
#             y = Part,
#             by = c("PartNumber" = "PartNumber"))
# write.csv(ehr.common.loinc.methods, "ehr_common_loinc_methods.csv")

####    ####    ####    ####

# temp <- component.analyte.core.any.map.meth$loinc.provided.secondary
#
#
# # only single mappers are returned
# temp <- component.analyte.core.any.map.meth$paths.from.leaves
# temp <- make.table.frame(temp$CODE)


####

# bicarbonate

####

# will want to add methods and challenges in "soon"?
# adjustments and numerators seem like lower priorities
# add ordinal soon
# nominal later
# "Nom", "Ord"
# "24H"/24 hours sample collection objective is a turbo term
# divison by XXX tranformation process terms are from turbo
sure.things <-
  max.one.gene.wide[is.na(max.one.gene.wide$METHOD_TYP) &
                      is.na(max.one.gene.wide$challenge) &
                      is.na(max.one.gene.wide$adjustment) &
                      is.na(max.one.gene.wide$analyte.numerator) &
                      max.one.gene.wide$analyte.core  %in% loinc_to_obo_mapping_reviewed$PartNumber &
                      max.one.gene.wide$PartDisplayName.SCALE_TYP %in% c("Qn") &
                      max.one.gene.wide$PartDisplayName.analyte.suffix %in% c(NA, "IgA", "IgM", "Ab", "IgG", "IgE") &
                      # max.one.gene.wide$PartDisplayName.analyte.divisor %in% c(NA, "100 cells", "100 leukocytes", "Creatinine") &
                      is.na(max.one.gene.wide$PartDisplayName.analyte.divisor) &
                      max.one.gene.wide$PartName.TIME_ASPCT %in% c(NA, "Pt") ,]

# temp <-
#   as.data.frame.matrix(
#     table(
#       sure.things$PartDisplayName.PROPERTY,
#       sure.things$PartDisplayName.SCALE_TYP
#     )
#   )

# prop.scale.utilizastion <-
#   aggregate(
#     sure.things$RESULT_COUNT,
#     by = list(
#       sure.things$PartDisplayName.PROPERTY,
#       sure.things$PartDisplayName.SCALE_TYP
#     ),
#     FUN = sum
#   )
#
# colnames(prop.scale.utilizastion) <-
#   c("PartDisplayName.PROPERTY",
#     "PartDisplayName.SCALE_TYP",
#     "RESULT_COUNT.sum")

prop.utilizastion <-
  aggregate(sure.things$RESULT_COUNT,
            by = list(sure.things$PROPERTY),
            FUN = sum)

colnames(prop.utilizastion) <-
  c("PROPERTY",
    "RESULT_COUNT.sum")

prop.utilizastion <-
  inner_join(x = prop.utilizastion,
             y = Part,
             by = c("PROPERTY" = "PartNumber"))

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

cd_markers <-
  read_csv(cd_markersmapping_fn)

temp <-
  apply(
    X = cd_markers,
    MARGIN = 1,
    FUN = function(current_marker) {
      print(current_marker[["PartNumber"]])
      print(current_marker[["PartDisplayName"]])
      
      inner_temp <-
        ols.serch.term.labels.universal(
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

loinc_combo_analytes <-
  read_csv(loinc_combo_analytes_fn)
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

write.csv(loinc_combo_analytes,
          "loinc_combo_analytes_new.csv",
          row.names = FALSE)

#### many "cells by surface marker" and "combo analytes" are already integrated into
####   loinc_to_obo_mapping_reviewed.csv

mappings.with.hints <-
  component.analyte.core.any.map.meth$any.method

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

# several of the hinted analytes have already been integrated into loinc_to_obo_mapping_reviewed.csv
