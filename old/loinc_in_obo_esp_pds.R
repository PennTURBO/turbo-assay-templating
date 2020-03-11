library(SPARQL)
library(config)
library(httr)
library(jsonlite)
library(rJava)
library(RJDBC)

# source("/turbo_graphdb_setup.R")

source("~/../loinc_in_obo/turbo_graphdb_setup.R")

config.yaml.file <- "~/labs_config.yaml"

sparql.namespaces <- c(
  'dc',
  '<http://purl.org/dc/elements/1.1/>',
  'rdfs',
  '<http://www.w3.org/2000/01/rdf-schema#>',
  'mydata',
  '<http://example.com/resource/>',
  'owl',
  '<http://www.w3.org/2002/07/owl#>',
  'skos',
  '<http://www.w3.org/2004/02/skos/core#>'
)

jdbcDriver <-
  JDBC("oracle.jdbc.OracleDriver", classPath = "C:/ojdbc8.jar")

config.bootstrap <-
  config::get(file = config.yaml.file)

umls.semantic.types.uri <- config.bootstrap$umls.semantic.types.uri
loinc.uri <- config.bootstrap$loinc.uri

selected.relational.ehr.configuration <-
  config::get(config = "pds",
              file = config.yaml.file)

selected.gdb.configuration <-
  config::get(config = config.bootstrap$selected.gdb.configuration,
              file = config.yaml.file)

graphdb.address.port <-
  selected.gdb.configuration$graphdb.address.port
selected.repo <- selected.gdb.configuration$selected.repo
api.user <- selected.gdb.configuration$api.user
api.pass <- selected.gdb.configuration$api.pass

saved.authentication <-
  authenticate(api.user, api.pass, type = "basic")

url.post.endpoint <-
  paste0(graphdb.address.port,
         "/rest/data/import/upload/",
         selected.repo,
         "/url")

update.endpoint <-
  paste0(graphdb.address.port,
         "/repositories/",
         selected.repo,
         "/statements")

select.endpoint <-
  paste0(graphdb.address.port, "/repositories/", selected.repo)

# does this ever get used?
bp.api.key <- config.bootstrap$bp.api.key

debug.flag <- config.bootstrap$debug.flag

host.name <- selected.relational.ehr.configuration$host
database.name <- selected.relational.ehr.configuration$dbname
db.user <- selected.relational.ehr.configuration$dbuser
# postgres.passwd <- rstudioapi::askForPassword("Database password")
db.passwd <- selected.relational.ehr.configuration$dbpass
db.port <- selected.relational.ehr.configuration$port

jdbcConnection  <- dbConnect(
  jdbcDriver,
  paste0(
    "jdbc:oracle:thin:@",
    host.name,
    ":",
    db.port,
    "/",
    database.name
  ),
  db.user,
  db.passwd
)

# ~ 1 hour?

# ehr.lab.freqs <- dbGetQuery(
#   jdbcConnection,
#   'SELECT
# 	lr.FK_Lab_Result_Item_ID,
# 	rlri.LAB_RESULT_ITEM_CODE,
# 	rlri.LAB_RESULT_ITEM_DESCRIPTION,
# 	rlri.LOINC,
# 	COUNT(1)
# FROM
# 	mdm.LAB_RESULTS lr
# JOIN mdm.R_LAB_RESULT_ITEM rlri ON
# 	lr.FK_LAB_RESULT_ITEM_ID = rlri.PK_LAB_RESULT_ITEM_ID
# GROUP BY
# 	lr.FK_Lab_Result_Item_ID,
# 	rlri.LAB_RESULT_ITEM_CODE,
# 	rlri.LAB_RESULT_ITEM_DESCRIPTION,
# 	rlri.LOINC
# ORDER BY
# 	COUNT(1) DESC'
# )

# save.image("~/ehr_loinc_usage_df_only.Rdata")
load("~/ehr_loinc_usage_df_only.Rdata")

###   ###   ###

# CLEAR REPO

last.post.time <- Sys.time()

# post.res <- POST(update.endpoint,
#                  body = list(update = "clear all"),
#                  saved.authentication)
#
# # empty for sparql statement
# last.post.status <- rawToChar(post.res$content)


# or just do it as a sparql update

last.post.status <- update.statement <- "clear all"

sparql.result <-
  SPARQL(
    url =  update.endpoint,
    update = update.statement,
    curl_args = list(
      userpwd = paste0(api.user, ":", api.pass),
      httpauth = 1
    )
  )

# Warning message:
#   In testCurlOptionsInFormParameters(.params) :
#   Found possible curl options in form parameters: userpwd, httpauth

expectation <- NULL

monitor.named.graphs()

### LOINC

update.body <- paste0(
  '{
  "type":"url",
  "format":"text/turtle",
  "context": "http://purl.bioontology.org/ontology/LNC/",
  "data": "',
  loinc.uri,
  '"
}'
)

post.res <- POST(
  url.post.endpoint,
  body = update.body,
  content_type("application/json"),
  accept("application/json"),
  saved.authentication
)

last.post.status <- rawToChar(post.res$content)
last.post.time <- Sys.time()


### semantic types
# I used to include these in each bioportal export
#   and in a separate file
# now I'm doing separate file only
# no ontology name is asserted in the file
# graph/context name?
# https://bioportal.bioontology.org/ontologies/STY ?
# http://purl.bioontology.org/ontology/STY/ ?
# https://www.nlm.nih.gov/research/umls/META3_current_semantic_types.html ?


update.body <- paste0(
  '{
  "type":"url",
  "format":"text/turtle",
  "context": "http://purl.bioontology.org/ontology/STY/",
  "data": "',
  umls.semantic.types.uri,
  '"
}'
)


post.res <- POST(
  url.post.endpoint,
  body = update.body,
  content_type("application/json"),
  accept("application/json"),
  saved.authentication
)


last.post.status <- rawToChar(post.res$content)
last.post.time <- Sys.time()

###

expectation <-
  c(
    "http://purl.bioontology.org/ontology/LNC/",
    "http://purl.bioontology.org/ontology/STY/"
  )

monitor.named.graphs()

###

# PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
# PREFIX LNC: <http://purl.bioontology.org/ontology/LNC/>
# #PREFIX LNC: <http://purl.bioontology.org/ontology/LOINC/>
# PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

#        "2345-7"

prepre <- "
PREFIX LNC: <http://purl.bioontology.org/ontology/LNC/>
select
#(count(   ?notationlab) as ?count)
(str(?notationlab) as ?notation)
(str(?componentlab) as ?comp)
(str(?propertylab) as ?prop)
#(str(?timelab) as ?time)
(str(?systemlab) as ?sys)
#(str(?scalelab) as ?scale)
(str(?methodlab) as ?meth)
#(str(?classlab) as ?class)
where {
    graph LNC: {
        ?LoincTerm skos:notation ?notationlab ;
                   #                   LNC:has_class ?LoincClass ;
                   LNC:has_component ?LoincComponent ;
                   LNC:has_property ?LoincProperty ;
                   LNC:has_scale <http://purl.bioontology.org/ontology/LNC/LP7753-9> ;
                   LNC:has_system ?LoincSystem ;
                   LNC:has_time_aspect <http://purl.bioontology.org/ontology/LNC/LP6960-1> .
        optional {
            ?LoincTerm LNC:has_method ?LoincMethod .
            ?LoincMethod skos:prefLabel ?methodlab .
        }
        ?LoincComponent skos:prefLabel ?componentlab .
        ?LoincProperty skos:prefLabel ?propertylab .
        #        ?LoincTime skos:prefLabel ?timelab .
        ?LoincSystem skos:prefLabel ?systemlab .
        #        ?LoincScale skos:prefLabel ?scalelab .
        #        ?LoincClass skos:prefLabel ?classlab .
    }
}
# order by ?notationlab
# limit  100000
# offset 100000
"

#     99999
# 9 minutes
# 8 GB

prepre <-
  gsub(pattern = " +",
       replacement = " ",
       x = prepre)

prefixed <- paste0(sparql.prefixes, "\n", prepre, "\n")

Sys.time()
sparql.res.list <- SPARQL(
  url =  select.endpoint,
  query = prefixed,
  # ns = sparql.namespaces,
  curl_args = list(
    userpwd = paste0(api.user, ":", api.pass),
    httpauth = 1
    # 'Accept-Encoding' = 'gzip, deflate'
  )
)

loinc.part.frame <- sparql.res.list$results
Sys.time()

# http://purl.bioontology.org/ontology/LNC/LP6960-1 skos:prefLabel	"Pt"@en
# http://purl.bioontology.org/ontology/LNC/LP7753-9	skos:prefLabel	"Qn"@en
