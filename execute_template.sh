cat lobo_template_headers_only.csv lobo_template_headerless.csv > lobo_template.csv
#robot template --input ~/obi/obi.owl --prefix "lobo: http://example.com/lobo/" -t lobo_template.csv --output lobo.ttl
robot template --input ~/Turbo-Ontology/ontologies/animals_robot_transition/turbo_remerged.ttl --prefix "lobo: http://example.com/lobo/" -t lobo_template.csv --output lobo.ttl

