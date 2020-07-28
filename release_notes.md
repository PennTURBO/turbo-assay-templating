## 2020-07-28, unversioned

- Only ~ 400 rows, due to removal of any TURBO terms
- Free free from mentions of non-OBO authorities 

# TODO
- add versioning to 
    - **upstream** TURBO ontology artifact (version IRI) and python script (repo tag)
    - assay template generating code including execute_robot.sh (repo tag)
    - turbo assay turtle file... version IRI from annotations file?
- improve readability of labels and definitions

# Additional Challenges
- 

# Causes of failures to build assay ontology from tamplate:
- Mentioned terms not in ontology
- Single quotes in term labels, especially in combination with parentheses
- Core analytes that have single quotes in thier labels, like 5' nucleotidase

# The TURBO ontology with OBO imports was found inconsistent by Hermit
Solved by removing several terms about birth

[DOB about day of birth](http://transformunify.org/ontologies/TURBO_0010198)
[DOB TMD 946](http://transformunify.org/ontologies/TURBO_0010238)
[birth instant TPO birth day](http://transformunify.org/ontologies/TURBO_0010202)
[birth instant starts NTR](http://transformunify.org/ontologies/TURBO_0010245)
[neonate stage OTR NTR](http://transformunify.org/ontologies/TURBO_0010242)
[SNS OTR birth instant](http://transformunify.org/ontologies/TURBO_0010203)
[Patient 946](http://transformunify.org/ontologies/TURBO_0010210)
[Day of birth 946](http://transformunify.org/ontologies/TURBO_0010239)
[Birth instant 946](http://transformunify.org/ontologies/TURBO_0010240)
[Neonate stage (process) 946](http://transformunify.org/ontologies/TURBO_0010212)
[Neonate stage temporal region 946](http://transformunify.org/ontologies/TURBO_0010243)
[birth instant TPO NTR](http://transformunify.org/ontologies/TURBO_0010244)
[birth instant](http://transformunify.org/ontologies/TURBO_0010201)
[day of birth](http://transformunify.org/ontologies/TURBO_0010199)
[neonate temporal region](http://transformunify.org/ontologies/TURBO_0010200)
[Neonate stage temporal region 946](http://transformunify.org/ontologies/TURBO_0010243)
[Start of neonate stage/birth 946](http://transformunify.org/ontologies/TURBO_0010213)

The following classes were also foudn to be unsatifiable and were removed
[conclusionation frequency threshold value specification](http://transformunify.org/ontologies/TURBO_0001541)
[assay of specimen from organism](http://transformunify.org/ontologies/TURBO_0022089)


LOINC 94563-4, 'SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay' was not modelled. All assys with methods or challenges had been excluded.

| **LoincNumber** | **LongCommonName**                                           | **PartNumber** | **PartName**              | **PartCodeSystem** | **PartTypeName** | **LinkTypeName**  | **Property**                           |
| --------------- | ------------------------------------------------------------ | -------------- | ------------------------- | ------------------ | ---------------- | ----------------- | -------------------------------------- |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP417915-8     | SARS coronavirus 2 Ab.IgG | http://loinc.org   | COMPONENT        | Primary           | http://loinc.org/property/COMPONENT    |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP417540-4     | SARS coronavirus 2        | http://loinc.org   | COMPONENT        | SyntaxEnhancement | http://loinc.org/property/analyte-core |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP217197-5     | IA                        | http://loinc.org   | METHOD           | Primary           | http://loinc.org/property/METHOD_TYP   |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP217195-9     | PrThr                     | http://loinc.org   | PROPERTY         | Primary           | http://loinc.org/property/PROPERTY     |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP7751-3       | Ord                       | http://loinc.org   | SCALE            | Primary           | http://loinc.org/property/SCALE_TYP    |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP7576-4       | Ser/Plas                  | http://loinc.org   | SYSTEM           | Primary           | http://loinc.org/property/SYSTEM       |
| 94563-4         | SARS-CoV-2 (COVID-19) IgG Ab [Presence] in Serum or Plasma by Immunoassay | LP6960-1       | Pt                        | http://loinc.org   | TIME             | Primary           | http://loinc.org/property/TIME_ASPCT   |
