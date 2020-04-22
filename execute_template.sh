cat lobo_template_headers_only.csv lobo_template_headerless.csv > lobo_template.csv

# --input https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/animals_robot_transition/turbo_remerged.ttl \
# --input ~/Turbo-Ontology/ontologies/animals_robot_transition/turbo_remerged.ttl \
# --output-iri 'https://raw.githubusercontent.com/PennTURBO/loinc_in_obo/master/lobo.ttl' \

robot \
template \
--template lobo_template.csv \
--add-prefix "lobo: http://example.com/lobo/"  \
--add-prefix "loinc: http://purl.bioontology.org/ontology/LNC/" \
--add-prefix "ghi: https://github.com/obi-ontology/obi/issues/" \
--add-prefix "oboInOwl: http://www.geneontology.org/formats/oboInOwl#" \
--ontology-iri "https://raw.githubusercontent.com/PennTURBO/loinc_in_obo/master/lobo.ttl" \
--input ~/Turbo-Ontology/ontologies/animals_robot_transition/turbo_remerged.ttl \
--ancestors \
annotate \
--remove-annotations \
--annotation-file lobo_annotations.ttl \
--output lobo.ttl 

