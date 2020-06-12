robot \
template \
--template more-specimens.tsv \
--ontology-iri "https://raw.githubusercontent.com/PennTURBO/loinc_in_obo/master/more-specimens.ttl" \
--input ~/Turbo-Ontology/ontologies/robot_implementation/turbo_merged.ttl \
--ancestors \
annotate \
--remove-annotations \
--annotation-file  more-specimens_annotations.ttl \
--output more-specimens.ttl 

