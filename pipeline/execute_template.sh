cat config/turbo_assay_template_headers_only.csv build/turbo_assay_template_headerless.csv > build/turbo_assay_template.csv

rm -rf turbo_merged.ttl*
rm -rf build/turbo_merged.ttl*

wget https://raw.githubusercontent.com/PennTURBO/turbo-ontology/master/ontologies/robot_implementation/turbo_merged.ttl

mv turbo_merged.ttl build/turbo_merged.ttl 

# files only, no ontology URLs?
# --input https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/animals_robot_transition/turbo_remerged.ttl \

robot \
template \
--template build/turbo_assay_template.csv \
--add-prefix "lobo: http://example.com/lobo/"  \
--add-prefix "obo: http://purl.obolibrary.org/obo/" \
--add-prefix "loinc: http://purl.bioontology.org/ontology/LNC/" \
--add-prefix "ghi: https://github.com/obi-ontology/obi/issues/" \
--add-prefix "oboInOwl: http://www.geneontology.org/formats/oboInOwl#" \
--ontology-iri "https://raw.githubusercontent.com/PennTURBO/turbo-assay-templating/master/turbo_assays.ttl" \
--input build/turbo_merged.ttl \
--ancestors \
annotate \
--remove-annotations \
--annotation-file config/turbo_assay_annotations.ttl \
--output build/turbo_assays.ttl 

