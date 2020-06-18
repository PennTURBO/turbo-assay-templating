# get tissue/specimen relationships from TURBO

library(rrdf)
# library(rdflib)
# library(SPARQL)

# # get the IRIs for entities that have been asserted to bear an evalaunt role
# lobo_template_headerless.file <-
#   "~/loinc_in_obo/lobo_template_headerless.csv"
#
# lobo_template_headerless <-
#   read.csv(lobo_template_headerless.file, header = FALSE)
# evaluant.list <- unique(lobo_template_headerless$V17)

# better yet, get them form the reviewed systems table
turbo_loinc_systems_reviewed <-
  read.csv("turbo_loinc_systems_reviewed.csv")

evaluant.list <-
  turbo_loinc_systems_reviewed$iri[turbo_loinc_systems_reviewed$use == "y"]

evaluant.list <-
  sort(unique(unlist(strsplit(
    evaluant.list, "\\|"
  ))))

http.start <- grepl(pattern = "^http?://", x = evaluant.list)
evaluant.list <- cbind.data.frame(evaluant.list, http.start)

evaluants <- evaluant.list$evaluant.list

# # print(length(evaluants))
# # cat(evaluants)
#
# # assuming they're prefixed. could check whether they begin with http or not
# evaluants.split <- strsplit(evaluants, ":")
# evaluants.split <- do.call(rbind.data.frame, evaluants.split)
# colnames(evaluants.split) <- c("prefix", "rhs")
# evaluants.prefixes <- unique(evaluants.split$prefix)
#
# # hwo to process multiple prefixes?
# current.prefix <- evaluants.prefixes[[1]]
# forward.resolved <-
#   httr::GET(paste0("https://prefix.cc/", current.prefix, ".file.json"))
# forward.resolved <- rawToChar(forward.resolved$content)
# forward.resolved <- jsonlite::fromJSON(forward.resolved)
# # as.data.frame(forward.resolved)
# current.base <- forward.resolved[[current.prefix]]
#
# evaluants.fq <-
#   sub(
#     pattern = paste0(current.prefix, ":"),
#     replacement = current.base,
#     x = evaluants,
#     fixed = TRUE
#   )

# also reverse search
# also http://prefix.cc/popular/all.file.csv

# owl.pref <-

# https://raw.githubusercontent.com/obi-ontology/obi/master/obi.owl
# http://purl.obolibrary.org/obo/obi.owl"

# faster, but only takes filename arguments
# rrdf::load.rdf()
# rdflib will tak a URL

anat.spec.url <-
  "https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/robot_implementation/turbo_merged.ttl"

# anat.spec.ont <-
#   rdflib::rdf_parse(doc =
#                       "https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/robot_implementation/turbo_merged.ttl",
#                     format = "turtle")

ont.tf <- tempfile()
download.file(anat.spec.url, ont.tf)

anat.spec.ont <-
  rrdf::load.rdf(filename = ont.tf, format = "TURTLE")

# tissue.specimen.query <- "
# PREFIX owl: <http://www.w3.org/2002/07/owl#>
# PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
# PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
# select
# ?sfo ?sfolab ?tissue ?tisslab
# where {
#     ?sfo <http://www.w3.org/2002/07/owl#equivalentClass> ?x5 ;
#          rdfs:label ?sfolab .
#     ?x5 <http://www.w3.org/2002/07/owl#intersectionOf> ?x4  .
#     ?x4 rdf:first ?x3  .
#     ?x3 <http://www.w3.org/2002/07/owl#intersectionOf> ?x2 .
#     ?x2 rdf:first <http://purl.obolibrary.org/obo/OBI_0100051> ;
#         rdf:rest ?x1  .
#     ?x1 rdf:first ?x1f ;
#         rdf:rest ?x1r .
#     ?x1f  a owl:Restriction ;
#           owl:onProperty <http://purl.obolibrary.org/obo/OBI_0000312> ;
#           owl:someValuesFrom <http://purl.obolibrary.org/obo/OBI_0600005> .
#     ?x1r rdf:first ?x1rf .
#     ?x1rf a  owl:Restriction ;
#           owl:onProperty <http://purl.obolibrary.org/obo/RO_0001000> ;
#           owl:someValuesFrom ?tissue .
#     ?tissue rdfs:label ?tisslab .
# }"

tissue.specimen.query <- "
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
select
distinct
?specimen ?speclab ?anatterm ?anatlab
where {
  # ?s (rdfs:subClassOf|rdf:first|rdf:rest|owl:equivalentClass)* <http://purl.obolibrary.org/obo/OBI_0100051> .
  ?specimen (rdfs:subClassOf|rdf:first|rdf:rest|
               owl:equivalentClass|owl:intersectionOf)* <http://purl.obolibrary.org/obo/OBI_0100051> ;
  (rdfs:subClassOf|rdf:first|rdf:rest|
     owl:equivalentClass|owl:intersectionOf)* ?r ;
  rdfs:label ?l .
  ?r a owl:Restriction ;
  owl:onProperty <http://purl.obolibrary.org/obo/RO_0001000> ;
  owl:someValuesFrom ?anatterm .
  ?anatterm rdfs:label ?atl .
  bind(lcase(str(?l)) as ?speclab)
  bind(lcase(str(?atl)) as ?anatlab)
}"


# Sys.time()
# system.time(tissue.specimen.result <-
#               rdflib::rdf_query(rdf = anat.spec.ont, query = tissue.specimen.query))

Sys.time()
system.time(
  tissue.specimen.result <-
    rrdf::sparql.rdf(model = anat.spec.ont, sparql = tissue.specimen.query)
)
tissue.specimen.result <- as.data.frame(tissue.specimen.result)



# user  system elapsed
# 52.728   0.079  52.978

evaluants.only <-
  sort(setdiff(evaluants, tissue.specimen.result$anatterm))

print(evaluants.only)

write.csv(x = tissue.specimen.result, file = "tissue_specimen_result.csv")

# #### is this necessary?
#
# Sys.time()
# system.time(
#   uber.ext <-
#     rdflib::rdf_parse(doc = "http://ontologies.berkeleybop.org/uberon/ext.owl", format = "rdfxml")
# )
#
#
#
# user  system elapsed
# # 6.604   0.463  10.076
#
# # # actually, its better to submit prefixed terms
# # #   because they take up less space in the CURL payload!
# # evaluant.valblock <-
# #   paste0("{ <", paste0(evaluants.only, collapse = "> <"), "> }")
# #
# # uber.lab.query <-
# #   paste0(
# #     "
# #     select ?s ?l
# # where {
# #   values ?s ",
# #     evaluant.valblock,
# #     " ?s <http://www.w3.org/2000/01/rdf-schema#label> ?l .
# # } "
# #   )
# #
# # uber.lab.query <-
# #   paste0(
# #     "select * where { values ?s {
# #     <http://purl.obolibrary.org/obo/UBERON_0010974> <http://purl.obolibrary.org/obo/UBERON_0000178>
# #     } . ?s <http://www.w3.org/2000/01/rdf-schema#label> ?l . } limit 3 "
# #   )
#
# uber.lab.query <-
#   paste0("select ?s ?l  where {  ?s <http://www.w3.org/2000/01/rdf-schema#label> ?l  } ")
#
# Sys.time()
# system.time(uber.lab.result <-
#               rdflib::rdf_query(rdf = uber.ext, query = uber.lab.query))
#
# evaluants.only.labels <-
#   uber.lab.result[uber.lab.result$s %in% evaluants.only , ]
#
# setdiff(evaluants.only, evaluants.only.labels$s)
