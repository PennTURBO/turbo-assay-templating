The TURBO clinical assays module currently defines assays in terms of tissues etc., but the standard pattern in OBI uses specimens-from-organism.

### Which tissues etc. are used in the TURBO clinical assays module?

```R
lobo_template_headerless <-
  read.csv("~/loinc_in_obo/lobo_template_headerless.csv", header = FALSE)

evaluants <- unique(lobo_template_headerless$V17)
evaluants <- sort(unique(unlist(strsplit(evaluants, "\\|"))))
print(length(evaluants))
cat(evaluants)
```



Look up those terms in Uberon 



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



Compare to [specimens available in OBI](https://github.com/PennTURBO/loinc_in_obo/blob/master/specimens_in_ontobee.csv) and highlight absent or problematic specimen terms:



- Specimen term ready to connect to clinical assays module
- **Tissue term exists but specimen term absent from OBI**
- *Specimen term present in OBI, but uses subClassOf logical definition instead of equivalentCass logical definition*

| s                                                 | l                        |
| ------------------------------------------------- | ------------------------ |
| http://purl.obolibrary.org/obo/UBERON_0000173     | amniotic fluid           |
| **http://purl.obolibrary.org/obo/UBERON_0013755** | **arterial blood**       |
| *http://purl.obolibrary.org/obo/UBERON_0000178*   | *blood*                  |
| **http://purl.obolibrary.org/obo/UBERON_0001969** | **blood plasma**         |
| **http://purl.obolibrary.org/obo/UBERON_0001977** | **blood serum**          |
| **http://purl.obolibrary.org/obo/UBERON_0006314** | **bodily fluid**         |
| **http://purl.obolibrary.org/obo/UBERON_0013757** | **capillary blood**      |
| **http://purl.obolibrary.org/obo/UBERON_0001359** | **cerebrospinal fluid**  |
| http://purl.obolibrary.org/obo/UBERON_0001268     | peritoneal fluid         |
| http://purl.obolibrary.org/obo/UBERON_0001087     | pleural fluid            |
| http://purl.obolibrary.org/obo/UBERON_0001836     | saliva                   |
| http://purl.obolibrary.org/obo/UBERON_0001090     | synovial fluid           |
| **http://purl.obolibrary.org/obo/UBERON_0012168** | **umbilical cord blood** |
| *http://purl.obolibrary.org/obo/UBERON_0001088*   | *urine*                  |
| **http://purl.obolibrary.org/obo/UBERON_0013756** | **venous blood**         |



Also need to create a specimen term for erythrocyte, http://purl.obolibrary.org/obo/CL_0000232