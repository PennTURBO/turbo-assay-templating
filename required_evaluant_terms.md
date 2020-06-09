```R
lobo_template_headerless <-
  read.csv("~/loinc_in_obo/lobo_template_headerless.csv", header = FALSE)

evaluants <- unique(lobo_template_headerless$V17)
evaluants <- sort(unique(unlist(strsplit(evaluants, "\\|"))))
print(length(evaluants))
cat(evaluants)
```

----

```SPARQL
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
select * where {
    values ?s {
        obo:CL_0000232 
        obo:UBERON_0000173 obo:UBERON_0000178 obo:UBERON_0001087 obo:UBERON_0001088 obo:UBERON_0001090 
        obo:UBERON_0001268 obo:UBERON_0001359 obo:UBERON_0001836 obo:UBERON_0001969 obo:UBERON_0001977 
        obo:UBERON_0006314 obo:UBERON_0012168 obo:UBERON_0013755 obo:UBERON_0013756 obo:UBERON_0013757
    }
    ?s rdfs:label ?l
}
```

----

| s                                             | l                    |
| --------------------------------------------- | -------------------- |
| http://purl.obolibrary.org/obo/UBERON_0000173 | amniotic fluid       |
| http://purl.obolibrary.org/obo/UBERON_0013755 | arterial blood       |
| http://purl.obolibrary.org/obo/UBERON_0000178 | blood                |
| http://purl.obolibrary.org/obo/UBERON_0001969 | blood plasma         |
| http://purl.obolibrary.org/obo/UBERON_0001977 | blood serum          |
| http://purl.obolibrary.org/obo/UBERON_0006314 | bodily fluid         |
| http://purl.obolibrary.org/obo/UBERON_0013757 | capillary blood      |
| http://purl.obolibrary.org/obo/UBERON_0001359 | cerebrospinal fluid  |
| http://purl.obolibrary.org/obo/UBERON_0001268 | peritoneal fluid     |
| http://purl.obolibrary.org/obo/UBERON_0001087 | pleural fluid        |
| http://purl.obolibrary.org/obo/UBERON_0001836 | saliva               |
| http://purl.obolibrary.org/obo/UBERON_0001090 | synovial fluid       |
| http://purl.obolibrary.org/obo/UBERON_0012168 | umbilical cord blood |
| http://purl.obolibrary.org/obo/UBERON_0001088 | urine                |
| http://purl.obolibrary.org/obo/UBERON_0013756 | venous blood         |
