# make all file imports unique and complete

# in addition to tracking the number of result items, order items and specimens,
#  we should also track the total number of results that can be mapped to classes in the resulting ontology

## what do initial or terminal characters like -, _, ., z mean in orders?

# could change all of the counts to by-patient instead of by-order
# that would be more similar to the medication mapping approach,  but probably much slower

source(
  "https://raw.githubusercontent.com/PennTURBO/turbo-globals/master/turbo_R_setup.R"
)

library(skimr)
library(googlesheets4)

unwanted.cols <-
  c("ACTIVE_YN",
    "MDM_LAST_UPDATE_DATE",
    "MDM_INSERT_UPDATE_FLAG")

# rlri_obo.fp <- "rlri_obo_loincsafe.csv"
rlri_obo.fp <- "rlri_res_20201105.csv"

roi_normalized.fp <- "roi_label_normalization.csv"

obi_assay_headers.fp <-  "obi_assay_headers.txt"

specimen_with_method_table_annotated.fp <-
  "specimen_with_method_table_annotated.csv"

dont.trust.manual.spreadsheets <- function(spreadsheet) {
  spreadsheet <- unique(spreadsheet[complete.cases(spreadsheet), ])
}

# all lowercase with no leading or trailing space
#   could also remove multi whitespace
tidy.for.merge <- function(input.vector) {
  temp <-
    tolower(input.vector)
  temp <-
    sub(pattern = "^\\s+",
        replacement = "",
        x = temp)
  temp <-
    sub(pattern = "\\s+$",
        replacement = "",
        x = temp)
  temp <-
    sub(pattern = "\\s+",
        replacement = " ",
        x = temp)
  return(temp)
}

# sometimes the ids.labels datfarame ends up empty
#  when this script is run beginning to end
# I added 10 second sleep breaks between the steps last time I ran it
#  and it turned out OK
# 5 seconds was ok at least once too
sleepsecs <- 2


####

empi.q <- "
PREFIX pmbb: <http://www.itmat.upenn.edu/biobank/>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX ontologies: <http://transformunify.org/ontologies/>
select
#distinct
#?d ?l
#*
?p ?r
where {
    graph pmbb:expanded {
        ?s a <http://purl.obolibrary.org/obo/IAO_0000028> ;
           obo:BFO_0000050 ?c ;
           <http://transformunify.org/ontologies/TURBO_0010094> ?r .
        ?c a <http://purl.obolibrary.org/obo/IAO_0000578> ;
           obo:BFO_0000051 ontologies:TURBO_0010295 ;
           obo:IAO_0000219 ?p .
        ?p a obo:NCBITaxon_9606 .
    }
}
"

empi.res <-
  q2j2df(query = empi.q, auth = saved.authentication, repo = "mam_drivetrain_prd")

empi.cat <- paste0(empi.res$r, collapse = ", ")

# SQL setup

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

temp.q <- "select * from MDM.R_LOCATION"

if (exists("pdsConnection")) {
  result <- tryCatch({
    temp <- dbGetQuery(pdsConnection, temp.q)
  }, error = function(e) {
    pdsConnection <-
      dbConnect(pdsDriver,
                pds.con.string,
                config$pds.user,
                config$pds.pw)
  })
} else {
  pdsConnection <-
    dbConnect(pdsDriver,
              pds.con.string,
              config$pds.user,
              config$pds.pw)
}

####

# what order primary types are observed for teh ~ 50 ish patients
# > table(order.items$ORDER_ITEM_PRIMARY_TYPE)
#
# CARDIOLOGY PROCEDURES            DIAGNOSTIC      ECHOCARDIOGRAPHY               IMAGING            LABORATORY            PROCEDURES   VASCULAR ULTRASOUND
#                    13                    10                     3                     5                   693                     7                     3

####

# what are ALL of the order primary types

# ORDER_ITEM_PRIMARY_TYPE frequencies in the source data

# MDM.R_LAB_RESULT_ITEM 47k rows... small and relevant enough that we can just pull the whole thing
# MDM.R_ORDER_ITEM 340k rows, including meds, microbiology, procedures

# X DISCHARGE MEDICATION	205186
# X MEDICATIONS	84225

# LABORATORY	20200
# DIAGNOSTIC	10469

# ? OTHER	10191

# PROCEDURES	9505
# IMAGING	4176

# X SUPPLIES	1438
# X IMMUNIZATION	810
# X CONSULTS	485
# X NURSING ORDERS	466

# CARDIOLOGY PROCEDURES	451
# ...
# VASCULAR ULTRASOUND	34
# ...
# ECHOCARDIOGRAPHY	18
# ...more...


# SQL queries

rlri.q <- "select * from MDM.R_LAB_RESULT_ITEM"

# seconds
print(Sys.time())
timed.system <- system.time(rlri.res <-
                              dbGetQuery(pdsConnection, rlri.q))
print(Sys.time())
print(timed.system)


rlri.res <- rlri.res[, setdiff(colnames(rlri.res), unwanted.cols)]

####

rlri.lr.count.q <- "
SELECT
	rlri.PK_LAB_RESULT_ITEM_ID,
	count(1) AS lab_result_count
FROM
	mdm.LAB_RESULTS lr
JOIN MDM.R_LAB_RESULT_ITEM rlri ON
	lr.FK_LAB_RESULT_ITEM_ID = rlri.PK_LAB_RESULT_ITEM_ID
GROUP BY
	rlri.PK_LAB_RESULT_ITEM_ID"

# 5 minutes
print(Sys.time())
timed.system <- system.time(rlri.lr.count.res <-
                              dbGetQuery(pdsConnection, rlri.lr.count.q))
print(Sys.time())
print(timed.system)


rlri.res <-
  full_join(x = rlri.res, y = rlri.lr.count.res, by = "PK_LAB_RESULT_ITEM_ID")

rlri.res <- rlri.res[!is.na(rlri.res$LAB_RESULT_COUNT), ]

####

rlri_obo <-
  read_csv(rlri_obo.fp)

rlri_obo <- dont.trust.manual.spreadsheets(rlri_obo)

#

rlri.res <-
  full_join(x = rlri.res, y = rlri_obo, by = "PK_LAB_RESULT_ITEM_ID")

####

propigate.by.loinc <-
  rlri.res[(!(is.na(rlri.res$obo))) &
             (!(is.na(rlri.res$LOINC))) ,]

propigate.by.loinc <- unique(propigate.by.loinc$LOINC)

propigate.by.loinc <-
  unique(rlri.res[rlri.res$LOINC %in% propigate.by.loinc, c("LOINC", "obo")])

colnames(propigate.by.loinc) <- c("LOINC", "loinc.obo")

propigate.by.loinc <-
  unique(propigate.by.loinc[complete.cases(propigate.by.loinc),])

pbl.loinc.tab <- make.table.frame(propigate.by.loinc$LOINC)
noisy.loincs <- pbl.loinc.tab$value[pbl.loinc.tab$count > 1]
pbl.obo.tab <- make.table.frame(propigate.by.loinc$loinc.obo)

propigate.by.loinc <-
  propigate.by.loinc[!(propigate.by.loinc$LOINC %in% noisy.loincs) , ]


rlri.res <-
  full_join(x = rlri.res,
            y = propigate.by.loinc,
            by = c("LOINC"))

rlri.res$obo[is.na(rlri.res$obo)] <-
  rlri.res$loinc.obo[is.na(rlri.res$obo)]

# rlri.res <-
#   rlri.res[(!(is.na(rlri.res$LAB_RESULT_COUNT))),]

rlri.res <-
  rlri.res[(!(is.na(rlri.res$obo))),]

# rlri.res <-
#   rlri.res[(!(is.na(rlri.res$LAB_RESULT_COUNT))) &
#              rlri.res$LAB_RESULT_COUNT > 1 ,]

####

# CHECK FOR CASES WHERE MORE THAN ONE ANALYTE/TARGET ENTITY HAS BEEN MAPPED TO A LOINC CODE
# WITHOUT EXAMINING WHAT THE LOINC CODE MEANS

# make a roi wanted cols vector

roi.q <- "
select
*
from MDM.R_ORDER_ITEM
where ORDER_ITEM_PRIMARY_TYPE in
( 'CARDIOLOGY PROCEDURES', 'DIAGNOSTIC', 'ECHOCARDIOGRAPHY', 'IMAGING', 'LABORATORY',  'PROCEDURES', 'VASCULAR ULTRASOUND')"

# seconds
print(Sys.time())
timed.system <- system.time(roi.res <-
                              dbGetQuery(pdsConnection, roi.q))
print(Sys.time())
print(timed.system)

roi.res <- roi.res[, setdiff(colnames(roi.res), unwanted.cols)]

####

roi.lr.count.q <- "
SELECT
	roi.PK_ORDER_ITEM_ID,
	count(1) AS lab_result_count
FROM
	mdm.LAB_RESULTS lr
JOIN mdm.R_ORDER_ITEM roi ON
	lr.FK_ORDER_ITEM_ID = roi.PK_ORDER_ITEM_ID
GROUP BY
	roi.PK_ORDER_ITEM_ID"

# 5 minutes
print(Sys.time())
timed.system <- system.time(roi.lr.count.res <-
                              dbGetQuery(pdsConnection, roi.lr.count.q))
print(Sys.time())
print(timed.system)

roi.res <-
  left_join(x = roi.res, y = roi.lr.count.res, by = "PK_ORDER_ITEM_ID")

roi.res <- roi.res[!is.na(roi.res$LAB_RESULT_COUNT) , ]

# roi.res <-
#   roi.res[(!(is.na(roi.res$LAB_RESULT_COUNT))) &
#             roi.res$LAB_RESULT_COUNT > 1 , ]

####

# looks like I created the mappings based on data that had already been stripped of 1-count rows

roi_normalized <-
  read_csv(roi_normalized.fp)

roi_normalized <- dont.trust.manual.spreadsheets(roi_normalized)

roi.res <-
  full_join(x = roi.res, y = roi_normalized, by = "PK_ORDER_ITEM_ID")

roi.res <- roi.res[!is.na(roi.res$tidied) , ]

####

# lab.result.instances.q <- "
# SELECT
# 	pe.EMPI,
# 	ROI.PK_ORDER_ITEM_ID,
# 	RLRI.PK_LAB_RESULT_ITEM_ID,
# -- MDM.LAB_RESULTS columns for class building, not instantiation
#  LR.PK_LAB_RESULT_ID,
# 	LR.ORDER_NAME,
# -- order text dates, draw instructions boilerplate, BUT ALSO occasionally some diagnosis/reason information
# --	LR.ORDER_TEXT,
# -- DIC called out
# --	LR.ORDER_SET_NAME,
# -- (coag?) times called out
# --  LR.ORDER_GROUP,
# 	LR.ORDER_TYPE,
# 	LR.SPECIMEN_COLLECTION_CLASS,
# 	LR.SPECIMEN,
# 	LR.SPECIMEN_COLLECTION_METHOD,
# 	LR.SPECIMEN_CONTAINER,
# --	LR.SPECIMEN_CONTAINER_VOLUME,
# --   FREETEXT       ALPHA CALCULATION     NUMERIC
# --        220        1487        5529       98799
# 	LR.RESULT_TYPE,
# 	LR.UNIT_OF_MEASURE,
# -- result text can clarify order/result in a small number of cases
# --	LR.RESULT_TEXT,
# -- instruments with location etc
# 	LR.RESULT_RESOURCE
# FROM
# 	mdm.patient_encounter pe
# -- inner join: don't report encounters with no lab results
# inner JOIN mdm.LAB_RESULTS lr ON
# 	pe.pK_Patient_Encounter_ID = lr.FK_Patient_Encounter_ID
# -- left join: dont skip any lab results
# full outer JOIN mdm.R_LAB_RESULT_ITEM rlri ON
# 	lr.FK_LAB_RESULT_ITEM_ID = RLRI.PK_LAB_RESULT_ITEM_ID
# -- left join: dont skip any lab results
# full outer JOIN mdm.R_ORDER_ITEM roi ON
# 	lr.FK_ORDER_ITEM_ID = roi.PK_ORDER_ITEM_ID
# WHERE
# 	pe.empi IN (
# 	1000107956, 1000278527, 1000281816, 1000367534, 1000405272, 1000546847, 1000854351, 1000897214, 1000960588, 1000960825,
# 	1001049004, 1001254964, 1001360256, 1001582758, 1001602831, 1001708405, 1001844191, 1002672027, 1002807554, 1003050271,
# 	1003624988, 1003806685, 1003848085, 1004053205, 1004152048, 1004284348, 1004485686, 1004530378, 8000125255, 8000417910,
# 	8440439230, 8441315009, 8442075693, 8444590624, 8444918080, 8446649113, 8449483403, 8449558691, 8450303284, 8451818660,
# 	8453511073 )
# 	AND lr.RESULT_VALUE_NUM IS NOT NULL
#   AND UNIT_OF_MEASURE IS NOT NULL
# 	AND RESULT_TYPE IN ('CALCULATION', 'NUMERIC')
# -- AND SPECIMEN IS NOT NULL
# "

lab.result.instances.q <- "
SELECT
PE.EMPI,
RLRI.PK_LAB_RESULT_ITEM_ID,
ROI.PK_ORDER_ITEM_ID,
LR.ORDER_NAME,
LR.ORDER_TYPE,
LR.PK_LAB_RESULT_ID,
LR.RESULT_RESOURCE,
LR.RESULT_TYPE,
LR.SPECIMEN,
LR.SPECIMEN_COLLECTION_CLASS,
LR.SPECIMEN_COLLECTION_METHOD,
LR.SPECIMEN_CONTAINER,
LR.UNIT_OF_MEASURE
FROM
	MDM.PATIENT_ENCOUNTER PE
INNER JOIN MDM.LAB_RESULTS LR ON
PE.PK_PATIENT_ENCOUNTER_ID = LR.FK_PATIENT_ENCOUNTER_ID
INNER JOIN MDM.R_LAB_RESULT_ITEM RLRI ON
LR.FK_LAB_RESULT_ITEM_ID = RLRI.PK_LAB_RESULT_ITEM_ID
INNER JOIN MDM.R_ORDER_ITEM ROI ON
LR.FK_ORDER_ITEM_ID = ROI.PK_ORDER_ITEM_ID
WHERE
	PE.EMPI IN (
	1000107956, 1000278527, 1000281816, 1000367534, 1000405272, 1000546847, 1000854351, 1000897214, 1000960588, 1000960825,
	1001049004, 1001254964, 1001360256, 1001582758, 1001602831, 1001708405, 1001844191, 1002672027, 1002807554, 1003050271,
	1003624988, 1003806685, 1003848085, 1004053205, 1004152048, 1004284348, 1004485686, 1004530378, 8000125255, 8000417910,
	8440439230, 8441315009, 8442075693, 8444590624, 8444918080, 8446649113, 8449483403, 8449558691, 8450303284, 8451818660,
	8453511073 )
AND LR.RESULT_VALUE_NUM IS NOT NULL
AND UNIT_OF_MEASURE IS NOT NULL
--	AND RESULT_TYPE IN ('CALCULATION',
--	'NUMERIC')
--	-- AND SPECIMEN IS NOT NULL"

# seconds for a group of 41 EMPIs
print(Sys.time())
timed.system <- system.time(lab.result.instances.res <-
                              dbGetQuery(pdsConnection, lab.result.instances.q))
print(Sys.time())
print(timed.system)

skim(lab.result.instances.res)

table(lab.result.instances.res$RESULT_TYPE, useNA = 'a')

lab.result.instances.res <-
  lab.result.instances.res[(!(
    lab.result.instances.res$RESULT_TYPE %in% c("ALPHA", "FREETEXT")
  )) , ]


# # 15 minutes ... only relevant for deciding which to map?
# specimen.method.counts.q <- "
# SELECT
# 	SPECIMEN ,
# 	SPECIMEN_COLLECTION_METHOD,
# 	COUNT(1) AS lab_result_count
# FROM
# 	mdm.LAB_RESULTS lr
# WHERE
# 	SPECIMEN IS NOT NULL
# GROUP BY
# 	SPECIMEN ,
# 	SPECIMEN_COLLECTION_METHOD
# "
#
# print(Sys.time())
# timed.system <- system.time(specimen.method.counts.res <-
#                               dbGetQuery(pdsConnection, specimen.method.counts.q))
# print(Sys.time())
# print(timed.system)

####

lab.result.instances.res$spec.tidy <-
  tidy.for.merge(lab.result.instances.res$SPECIMEN)

lab.result.instances.res$coll.meth.tidy <-
  tidy.for.merge(lab.result.instances.res$SPECIMEN_COLLECTION_METHOD)

####

specimen_with_method_table_annotated <-
  read_csv(specimen_with_method_table_annotated.fp)

# specimen_with_method_table_annotated <-
#   dont.trust.manual.spreadsheets(specimen_with_method_table_annotated)

specimen_with_method_table_annotated <-
  unique(specimen_with_method_table_annotated[!is.na(specimen_with_method_table_annotated$obo),
                                              c("spec.tidy", "coll.meth.tidy", "obo")])

specimen_with_method_table_annotated$coll.meth.tidy[is.na(specimen_with_method_table_annotated$coll.meth.tidy)] <-
  ""

####

lab.result.instances.working <- lab.result.instances.res
lab.result.instances.working$coll.meth.tidy[is.na(lab.result.instances.working$coll.meth.tidy)] <-
  ""
lab.result.instances.working$spec.tidy[is.na(lab.result.instances.working$spec.tidy)] <-
  ""


lab.result.instances.working <-
  left_join(
    x = lab.result.instances.working,
    y = specimen_with_method_table_annotated,
    by = c("spec.tidy", "coll.meth.tidy"),
    suffix = c(".instance", ".specimen")
  )

lab.result.instances.working <-
  left_join(
    x = lab.result.instances.working,
    y = roi.res,
    by = c("PK_ORDER_ITEM_ID"),
    suffix = c(".specimen", ".order")
  )

lab.result.instances.working <-
  left_join(
    x = lab.result.instances.working,
    y = rlri.res,
    by = c("PK_LAB_RESULT_ITEM_ID"),
    suffix = c(".order", ".result")
  )

# require that a term has been found for the target entity
lab.result.instances.working <-
  lab.result.instances.working[!is.na(lab.result.instances.working$EMPI) &
                                 !is.na(lab.result.instances.working$obo.result) , ]

lab.result.instances.working <-
  unique(lab.result.instances.working[, c("tidied", "obo.result", "obo.order")])

# something got out of whack with the suffices
colnames(lab.result.instances.working) <-
  c("tidied", "obo.result", "obo.specimen")

####

temp <-
  strsplit(lab.result.instances.working$obo.result,
           split = ':',
           fixed = TRUE)
temp <- do.call(rbind.data.frame, temp)
lab.result.instances.working$result.ontology <- temp[, 1]

temp <-
  strsplit(lab.result.instances.working$obo.specimen,
           split = ':',
           fixed = TRUE)
temp <- do.call(rbind.data.frame, temp)
lab.result.instances.working$specimen.ontology <- temp[, 1]

####

temp <- colnames(lab.result.instances.working)
temp <- temp[grepl(pattern = "\\.ontology$", x = temp)]

lab.result.instances.working.ontologies <-
  make.table.frame(unlist(as.matrix(lab.result.instances.working[, temp])))

####

# making lots of OBO assumptions about the following term labeling patterns

placeholder <-
  lapply(
    sort(unique(
      lab.result.instances.working.ontologies$value
    )),
    FUN = function(current.ontology) {
      print(current.ontology)
      q <- paste0("select
*
where
{
    graph <http://purl.obolibrary.org/obo/",
  tolower(current.ontology) , ".owl> {
        ?term <http://www.w3.org/2000/01/rdf-schema#label> ?label .
    }
    bind(lang(?label) as ?lablang)
    filter(isiri(?term))
}")
      res <-
        q2j2df(endpoint = "http://localhost:7200",
               query = q,
               repo = "bottom_up")
      res$onto.abbrev <- current.ontology
      return(res)
    })

names(placeholder) <-
  sort(unique(lab.result.instances.working.ontologies$value))

Sys.sleep(sleepsecs)

####

ids.labels <- do.call(rbind.data.frame, placeholder)

Sys.sleep(sleepsecs)

ids.labels$term.owner <-
  mapply(grepl, pattern = ids.labels$onto.abbrev, x = ids.labels$term)

Sys.sleep(sleepsecs)

table(ids.labels$term.owner, useNA = 'a')

ids.labels <- ids.labels[ids.labels$term.owner , ]

Sys.sleep(sleepsecs)

ids.labels$termnum <-
  sub(pattern = "^.*_",
      replacement = "",
      x = ids.labels$term)

Sys.sleep(sleepsecs)

ids.labels$id.recon <-
  paste(ids.labels$onto.abbrev, ids.labels$termnum, sep = ":")

Sys.sleep(sleepsecs)

ids.labels$sq.flag <- grepl(pattern = "'", x = ids.labels$label)

Sys.sleep(sleepsecs)

ids.labels$dq.flag <- grepl(pattern = '"', x = ids.labels$label)

Sys.sleep(sleepsecs)

ids.labels <-
  ids.labels[(!(ids.labels$sq.flag)) &
               (!(ids.labels$dq.flag))  ,]

####

# take hard coding out
lab.result.instances.working <-
  lab.result.instances.working[lab.result.instances.working$result.ontology %in% c("CHEBI", "PR") , ]

lab.result.instances.working <-
  lab.result.instances.working[(is.na(lab.result.instances.working$specimen.ontology)) |
                                 lab.result.instances.working$specimen.ontology == "UBERON" ,]

####

lab.result.instances.working <-
  left_join(
    x = lab.result.instances.working,
    y = ids.labels,
    by = c("obo.result" = "id.recon") ,
    suffix = c(".instances", ".result")
  )


lab.result.instances.working <-
  left_join(
    x = lab.result.instances.working,
    y = ids.labels,
    by = c("obo.specimen" = "id.recon") ,
    suffix = c(".result", ".specimen")
  )

# lab.result.instances.working <-
#   lab.result.instances.working[(!(is.na(
#     lab.result.instances.working$label.specimen
#   ))), ]

####

lab.result.instances.working <-
  unique(lab.result.instances.working[, c("tidied", "label.specimen", "label.result")])

lab.result.instances.working <-
  lab.result.instances.working[(!(is.na(
    lab.result.instances.working$label.result
  ))), ]

#### so far so good

# take a matrix approach from the beginning

lriw.len <- nrow(lab.result.instances.working)
blanks <- rep("", lriw.len)

curation_status <-
  alt_term <-
  usage_example <-
  editor_note <-
  term_requestor <-
  term_tracker_id <-
  logical_note <-
  mat_proc_tech <-
  det_tech <-
  device <-
  reagent <-
  mol_lab <-
  input <-
  output <- targ_ent <- objective <- assoc_axes <- dbxref <- blanks


# paramaterize the following strings
obi.prefix <- "http://purl.obolibrary.org/obo/OBI_"
obi.start <- 3100001
obi.stop <- obi.start + lriw.len - 1
iri.values <- seq(obi.start, obi.stop)


ID <- paste0(obi.prefix, iri.values)

label <- rep("assay for", lriw.len)
current.flag <- (!(is.na(lab.result.instances.working$tidied)))
label[current.flag] <-
  paste(lab.result.instances.working$tidied[current.flag], label[current.flag])
label <- paste(label, lab.result.instances.working$label.result)
label.suffix <- rep("", lriw.len)
current.flag <-
  (!(is.na(
    lab.result.instances.working$label.specimen
  )))
label.suffix[current.flag] <-
  paste(" in", lab.result.instances.working$label.specimen[current.flag])
label <- paste0(label, label.suffix)

def <-
  paste(
    "An analyte assay that measures the abundance of",
    lab.result.instances.working$label.result
  )
current.flag <-
  (!(is.na(
    lab.result.instances.working$label.specimen
  )))
def[current.flag] <-
  paste(def[current.flag], "in", lab.result.instances.working$label.specimen[current.flag])
current.flag <- (!(is.na(lab.result.instances.working$tidied)))
def[current.flag] <-
  paste(def[current.flag],
        "as part of a",
        lab.result.instances.working$tidied[current.flag],
        "order")

def_source <- rep("PennTURBO team", lriw.len)

term_editor <-
  rep(
    "Mark A. Miller, ORCID:0000-0001-9076-6066 & Chris Stoeckert, ORCID ORCID:0000-0002-5714-991X",
    lriw.len
  )

logical_type <- rep("subclass", lriw.len)
parent_class <- rep("'analyte assay'", lriw.len)

evaluant <- blanks
current.flag <-
  (!(is.na(
    lab.result.instances.working$label.specimen
  )))
evaluant[current.flag] <-
  paste0("'", lab.result.instances.working$label.specimen[current.flag], "'")


evaluant <- blanks
current.flag <-
  (!(is.na(
    lab.result.instances.working$label.specimen
  )))
evaluant[current.flag] <-
  paste0(
    "(specimen and ('derives from' some ('part of' some '",
    lab.result.instances.working$label.specimen[current.flag],
    "')))"
  )

analyte <-
  paste0("'", lab.result.instances.working$label.result, "'")

template.for.robot <-
  cbind.data.frame(
    ID,
    label,
    curation_status,
    alt_term,
    def,
    def_source,
    usage_example,
    editor_note,
    term_editor,
    term_requestor,
    term_tracker_id,
    logical_type,
    logical_note,
    parent_class,
    mat_proc_tech,
    det_tech,
    evaluant,
    analyte,
    device,
    reagent,
    mol_lab,
    input,
    output,
    targ_ent,
    objective,
    assoc_axes,
    dbxref
  )

template.for.robot <- as.matrix(template.for.robot)
colnames(template.for.robot) <- NULL

obi_assay_headers <-
  read_delim(
    obi_assay_headers.fp,
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE
  )

obi_assay_headers <- as.matrix(obi_assay_headers)
obi_assay_headers[is.na(obi_assay_headers)] <- ""
colnames(obi_assay_headers) <- NULL

template.for.robot <-
  as.data.frame(rbind(obi_assay_headers, template.for.robot))

template.for.robot <-
  template.for.robot[template.for.robot$V18 != "'NA'" ,]

unique.label.check <- make.table.frame(template.for.robot$V2)

print(max(unique.label.check$count))

analyte.tab <- make.table.frame(template.for.robot$V18)

print(max(analyte.tab$count))

####

write_delim(
  x = template.for.robot,
  path = "bottom_up_with_defs.tsv",
  delim = "\t",
  col_names = FALSE
)

####

# googlewb.name <- "TURBO assay template - Clinically relevant analytes"
# sheet.name <- "bottom-up-template"
# googlewb <-
#   gs4_create(googlewb.name)

googlewb <-
  gs4_create()

temp <- template.for.robot
colnames(temp) <- temp[1, ]
temp <- temp[-1, ]

# first write times out?!
# prefers to name sheet and workbook iteself?
# sheet_write(data = template.for.robot, ss = googlewb, sheet = sheet.name)
sheet_write(data = temp, ss = googlewb)
# then get rid of "Sheet1"? or just write into it in the first place?

temp <- read_sheet(ss = googlewb)

#
# tidy_results <- text_df %>%
#   tidytext::unnest_tokens(word, text)
#
# tidy_results <- tidy_results %>%
#   anti_join(stop_words)
#
# tidy_counts <- tidy_results %>%
#   count(word, sort = TRUE)
#
# write.csv(x = rlri.res.annotated,
#           file = "rlri_annotated.csv",
#           row.names = FALSE)
