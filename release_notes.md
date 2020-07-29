## v2020-07-28 pre-release

- Only ~ 400 rows, due to removal of all TURBO terms for any part
    - using the very generic parent class 'assay', as 'assay of specimen from organism' is still just a TURBO term
    - TURBO terms that would need to be imported into an OBO ontology like OBI can be found by viewing pipeline/data/loinc_to_obo_mapping_reviewed.csv, filtering on `pure.obo == FALSE`
- Free from mentions of non-OBO authorities 

# TODO

- add versioning to 
  - **upstream** TURBO ontology artifact (version IRI) and python script (repo tag)
  - assay template generating code including execute_robot.sh (repo tag)
  - turbo assay turtle file... version IRI from annotations file?
- improve readability of labels and definitions
- the LOINC component part has been monolithically instantiated as a target entity. Checking if it is a molecular entity or atom after the fact will allow moving the OBO version of the component into the analyte column in the template
- still using TURBO IRIs... look for an acceptable range in the OBI namespace

# Additional Challenges

- 

# Causes of failures to build assay ontology from template:

- Mentioned terms not in ontology
- Single quotes in term labels, especially in combination with parentheses
- Core analytes that have single quotes in their labels, like 5' nucleotidase

# The TURBO ontology with OBO imports was found inconsistent by Hermit

Solved by removing several terms about birth


- [DOB about day of birth](http://transformunify.org/ontologies/TURBO_0010198)
- [DOB TMD 946](http://transformunify.org/ontologies/TURBO_0010238)
- [birth instant TPO birth day](http://transformunify.org/ontologies/TURBO_0010202)
- [birth instant starts NTR](http://transformunify.org/ontologies/TURBO_0010245)
- [neonate stage OTR NTR](http://transformunify.org/ontologies/TURBO_0010242)
- [SNS OTR birth instant](http://transformunify.org/ontologies/TURBO_0010203)
- [Patient 946](http://transformunify.org/ontologies/TURBO_0010210)
- [Day of birth 946](http://transformunify.org/ontologies/TURBO_0010239)
- [Birth instant 946](http://transformunify.org/ontologies/TURBO_0010240)
- [Neonate stage (process) 946](http://transformunify.org/ontologies/TURBO_0010212)
- [Neonate stage temporal region 946](http://transformunify.org/ontologies/TURBO_0010243)
- [birth instant TPO NTR](http://transformunify.org/ontologies/TURBO_0010244)
- [birth instant](http://transformunify.org/ontologies/TURBO_0010201)
- [day of birth](http://transformunify.org/ontologies/TURBO_0010199)
- [neonate temporal region](http://transformunify.org/ontologies/TURBO_0010200)
- [Neonate stage temporal region 946](http://transformunify.org/ontologies/TURBO_0010243)
- [Start of neonate stage/birth 946](http://transformunify.org/ontologies/TURBO_0010213)
- http://transformunify.org/ontologies/TURBO_0010213)


The following classes were also found to be unsatisfiable and were removed


[conclusionation frequency threshold value specification](http://transformunify.org/ontologies/TURBO_0001541)

[assay of specimen from organism](http://transformunify.org/ontologies/TURBO_0022089)


LOINC 94563-4, 'SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay' was not modelled. All assays with methods or challenges had been excluded.



| **PartNumber** | **PartName**              | **PartTypeName** | **Property**                           |
| -------------- | ------------------------- | ---------------- | -------------------------------------- |
| LP417915-8     | SARS coronavirus 2 Ab.IgG | COMPONENT        | http://loinc.org/property/COMPONENT    |
| LP417540-4     | SARS coronavirus 2        | COMPONENT        | http://loinc.org/property/analyte-core |
| LP217197-5     | IA                        | METHOD           | http://loinc.org/property/METHOD_TYP   |
| LP217195-9     | PrThr                     | PROPERTY         | http://loinc.org/property/PROPERTY     |
| LP7751-3       | Ord                       | SCALE            | http://loinc.org/property/SCALE_TYP    |
| LP7576-4       | Ser/Plas                  | SYSTEM           | http://loinc.org/property/SYSTEM       |
| LP6960-1       | Pt                        | TIME             | http://loinc.org/property/TIME_ASPCT   |

