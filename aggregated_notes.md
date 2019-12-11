# From https://github.com/PennTURBO/Turbo-Ontology/issues/4

Here are tabular and RDF representations of LOINC's _glucose level_. OBI has [measuring glucose concentration in blood serum assay](http://purl.obolibrary.org/obo/OBI_0000418) but no corresponding datum, so I have created [blood glucose measurement datum](http://transformunify.org/ontologies/TURBO_0010153) in the TURBO ontology. The TURBO ontology asserts the linkage between the LOINC term and the TURBO term via the `has OMOP concept ID` (TURBO_0010147) of 8840. Glucose level is roughly the 20th most common lab in our "PDS"  clinical data warehouse, out of roughly 30,000.

Up until now, I have only instantiated datums (and not assays) in the Synthea-OMOP-inspired Acorn files, but I think we **will** need to instantiate assays (and/or orders) in order to capture context like "who ordered the assay?" and "where was it performed?".

We can use this as a staring point for templating the laboratory assays (with requester context, etc.) and output datums./

According to https://loinc.org/get-started/loinc-term-basics/, the essential LOINC parts are 
- Component (Analyte)
    - The substance or entity being measured or observed.
- Property
    - The characteristic or attribute of the analyte.
- Time
    - The interval of time over which an observation was made.
- System (Specimen)
    - The specimen or thing upon which the observation was made.
- Scale
    - How the observation value is quantified or expressed: quantitative, ordinal, nominal.
- Method
    - _OPTIONAL_ A high-level classification of how the observation was made. Only needed when the technique affects the clinical interpretation of the results.


## From `Loinc_2.66_PartFile_3.0-Alpha.1/LoincPartLink.csv` at https://loinc.org/file-access/download-id/17948/

_I don't understand the **LinkTypeName**s yet_

LoincNumber | LongCommonName | PartNumber | PartName | PartCodeSystem | PartTypeName | LinkTypeName | Property
-- | -- | -- | -- | -- | -- | -- | --
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP14635-4 | Glucose | http://loinc.org | COMPONENT | Primary | http://loinc.org/property/COMPONENT
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP6827-2 | MCnc | http://loinc.org | PROPERTY | Primary | http://loinc.org/property/PROPERTY
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP6960-1 | Pt | http://loinc.org | TIME | Primary | http://loinc.org/property/TIME_ASPCT
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7576-4 | Ser/Plas | http://loinc.org | SYSTEM | Primary | http://loinc.org/property/SYSTEM
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7753-9 | Qn | http://loinc.org | SCALE | Primary | http://loinc.org/property/SCALE_TYP
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP14635-4 | Glucose | http://loinc.org | COMPONENT | DetailedModel | http://loinc.org/property/analyte
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP6827-2 | MCnc | http://loinc.org | PROPERTY | DetailedModel | http://loinc.org/property/PROPERTY
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP6960-1 | Pt | http://loinc.org | TIME | DetailedModel | http://loinc.org/property/time-core
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7576-4 | Ser/Plas | http://loinc.org | SYSTEM | DetailedModel | http://loinc.org/property/system-core
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7753-9 | Qn | http://loinc.org | SCALE | DetailedModel | http://loinc.org/property/SCALE_TYP
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP14635-4 | Glucose | http://loinc.org | COMPONENT | SyntaxEnhancement | http://loinc.org/property/analyte-core
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7479-1 | Plas | http://loinc.org | SYSTEM | Search | http://loinc.org/property/search
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7567-3 | Ser | http://loinc.org | SYSTEM | Search | http://loinc.org/property/search
2345-7 | Glucose [Mass/volume] in Serum or Plasma | LP7786-9 | CHEM | http://loinc.org | CLASS | Metadata | http://loinc.org/property/CLASS


## From http://data.bioontology.org/ontologies/LOINC/submissions/17/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb via https://bioportal.bioontology.org/ontologies/LOINC

```Turtle
<http://purl.bioontology.org/ontology/LNC/2345-7> a <http://www.w3.org/2002/07/owl#Class>;
  <http://bioportal.bioontology.org/ontologies/umls/cui> "C0484731";
  <http://bioportal.bioontology.org/ontologies/umls/hasSTY> <http://purl.bioontology.org/ontology/STY/T201>;
  <http://bioportal.bioontology.org/ontologies/umls/tui> "T201";
  <http://purl.bioontology.org/ontology/LNC/COMMON_ORDER_RANK> "120";
  <http://purl.bioontology.org/ontology/LNC/COMMON_SI_TEST_RANK> "0";
  <http://purl.bioontology.org/ontology/LNC/COMMON_TEST_RANK> "4";
  <http://purl.bioontology.org/ontology/LNC/CONSUMER_NAME> "Glucose lab";
  <http://purl.bioontology.org/ontology/LNC/EXAMPLE_UCUM_UNITS> "mg/dL";
  <http://purl.bioontology.org/ontology/LNC/EXAMPLE_UNITS> "mg/dL";
  <http://purl.bioontology.org/ontology/LNC/IEEE_CF_CODE10> "160197", "160369";
  <http://purl.bioontology.org/ontology/LNC/IEEE_DESCRIPTION> "Plasma glucose concentration taken from undetermined sample source",
    "Plasma glucose concentration taken from venous";
  <http://purl.bioontology.org/ontology/LNC/IEEE_DIM> "ML-3";
  <http://purl.bioontology.org/ontology/LNC/IEEE_REFID> "MDC_CONC_GLU_UNDETERMINED_PLASMA",
    "MDC_CONC_GLU_VENOUS_PLASMA";
  <http://purl.bioontology.org/ontology/LNC/IEEE_UOM_UCUM> "mg/dL g/L";
  <http://purl.bioontology.org/ontology/LNC/LCD> "Y";
  <http://purl.bioontology.org/ontology/LNC/LCL> "CHEM";
  <http://purl.bioontology.org/ontology/LNC/LCN> "1";
  <http://purl.bioontology.org/ontology/LNC/LCS> "ACTIVE";
  <http://purl.bioontology.org/ontology/LNC/LCT> "MIN";
  <http://purl.bioontology.org/ontology/LNC/LOINC_COMPONENT> "Glucose";
  <http://purl.bioontology.org/ontology/LNC/LOINC_PROPERTY> "MCnc";
  <http://purl.bioontology.org/ontology/LNC/LOINC_SCALE_TYP> "Qn";
  <http://purl.bioontology.org/ontology/LNC/LOINC_SYSTEM> "Ser/Plas";
  <http://purl.bioontology.org/ontology/LNC/LOINC_TIME_ASPECT> "Pt";
  <http://purl.bioontology.org/ontology/LNC/LOR> "Both";
  <http://purl.bioontology.org/ontology/LNC/LRN2> "Chemistry; Glu; Gluc; Glucoseur; Level; Mass concentration; Pl; Plasma; Plsm; Point in time; QNT; Quan; Quant; Quantitative; Random; SerP; SerPl; SerPlas; Serum; Serum or plasma; SR";
  <http://purl.bioontology.org/ontology/LNC/LUR> "Y";
  <http://purl.bioontology.org/ontology/LNC/TOP_2000_LAB_RESULTS_US> "4";
  <http://purl.bioontology.org/ontology/LNC/UNITS_AND_RANGE> "mg/dL:[74,106];";
  <http://purl.bioontology.org/ontology/LNC/UNIVERSAL_LAB_ORDERS_VALUE_SET> "TRUE";
  <http://purl.bioontology.org/ontology/LNC/VERSION_FIRST_RELEASED> "1.0";
  <http://purl.bioontology.org/ontology/LNC/VERSION_LAST_CHANGED> "2.34";
  <http://purl.bioontology.org/ontology/LNC/analyzes> <http://purl.bioontology.org/ontology/LNC/MTHU001012>;
  <http://purl.bioontology.org/ontology/LNC/has_class> <http://purl.bioontology.org/ontology/LNC/LP7786-9>;
  <http://purl.bioontology.org/ontology/LNC/has_component> <http://purl.bioontology.org/ontology/LNC/LP14635-4>;
  <http://purl.bioontology.org/ontology/LNC/has_property> <http://purl.bioontology.org/ontology/LNC/LP6827-2>;
  <http://purl.bioontology.org/ontology/LNC/has_scale> <http://purl.bioontology.org/ontology/LNC/LP7753-9>;
  <http://purl.bioontology.org/ontology/LNC/has_system> <http://purl.bioontology.org/ontology/LNC/LP7479-1>,
    <http://purl.bioontology.org/ontology/LNC/LP7567-3>, <http://purl.bioontology.org/ontology/LNC/LP7576-4>;
  <http://purl.bioontology.org/ontology/LNC/has_time_aspect> <http://purl.bioontology.org/ontology/LNC/LP6960-1>;
  <http://purl.bioontology.org/ontology/LNC/measures> <http://purl.bioontology.org/ontology/LNC/MTHU001675>;
  <http://purl.bioontology.org/ontology/LNC/member_of> <http://purl.bioontology.org/ontology/LNC/24320-4>,
    <http://purl.bioontology.org/ontology/LNC/24321-2>, <http://purl.bioontology.org/ontology/LNC/24322-0>,
    <http://purl.bioontology.org/ontology/LNC/24323-8>, <http://purl.bioontology.org/ontology/LNC/24353-5>,
    <http://purl.bioontology.org/ontology/LNC/24362-6>, <http://purl.bioontology.org/ontology/LNC/54232-4>,
    <http://purl.bioontology.org/ontology/LNC/54233-2>, <http://purl.bioontology.org/ontology/LNC/55399-0>,
    <http://purl.bioontology.org/ontology/LNC/70219-1>, <http://purl.bioontology.org/ontology/LNC/88843-8>,
    <http://purl.bioontology.org/ontology/LNC/89044-2>;
  <http://www.w3.org/2000/01/rdf-schema#subClassOf> <http://purl.bioontology.org/ontology/LNC/LP42107-0>,
    <http://purl.bioontology.org/ontology/LNC/MTHU000049>;
  <http://www.w3.org/2004/02/skos/core#altLabel> "Glucose SerPl-mCnc"@en, "Glucose [Mass/volume] in Serum or Plasma"@en,
    "Glucose:Mass Concentration:Point in time:Serum/Plasma:Quantitative"@en;
  <http://www.w3.org/2004/02/skos/core#notation> "2345-7";
  <http://www.w3.org/2004/02/skos/core#prefLabel> "Glucose:MCnc:Pt:Ser/Plas:Qn"@en .
```

## Same thing for LOINC's _patient body height_, which I have indirectly associated with OBI's _length measurement datum_ via OMOP concept 3036277.

LoincNumber | LongCommonName | PartNumber | PartName | PartCodeSystem | PartTypeName | LinkTypeName | Property
-- | -- | -- | -- | -- | -- | -- | --
8302-2 | Body height | LP64598-3 | Body height | http://loinc.org | COMPONENT | Primary | http://loinc.org/property/COMPONENT
8302-2 | Body height | LP6822-3 | Len | http://loinc.org | PROPERTY | Primary | http://loinc.org/property/PROPERTY
8302-2 | Body height | LP6960-1 | Pt | http://loinc.org | TIME | Primary | http://loinc.org/property/TIME_ASPCT
8302-2 | Body height | LP310005-6 | ^Patient | http://loinc.org | SYSTEM | Primary | http://loinc.org/property/SYSTEM
8302-2 | Body height | LP7753-9 | Qn | http://loinc.org | SCALE | Primary | http://loinc.org/property/SCALE_TYP
8302-2 | Body height | LP64598-3 | Body height | http://loinc.org | COMPONENT | DetailedModel | http://loinc.org/property/analyte
8302-2 | Body height | LP6822-3 | Len | http://loinc.org | PROPERTY | DetailedModel | http://loinc.org/property/PROPERTY
8302-2 | Body height | LP6960-1 | Pt | http://loinc.org | TIME | DetailedModel | http://loinc.org/property/time-core
8302-2 | Body height | LP7753-9 | Qn | http://loinc.org | SCALE | DetailedModel | http://loinc.org/property/SCALE_TYP
8302-2 | Body height | LP6985-8 | Patient | http://loinc.org | SUPER SYSTEM | DetailedModel | http://loinc.org/property/super-system
8302-2 | Body height | LP64598-3 | Body height | http://loinc.org | COMPONENT | SyntaxEnhancement | http://loinc.org/property/analyte-core
8302-2 | Body height | LP72681-7 | Height | http://loinc.org | COMPONENT | Search | http://loinc.org/property/search
8302-2 | Body height | LP281538-1 | - | http://loinc.org | SYSTEM | Search | http://loinc.org/property/search
8302-2 | Body height | LP7768-7 | BDYHGT.ATOM | http://loinc.org | CLASS | Metadata | http://loinc.org/property/CLASS

```Turtle
<http://purl.bioontology.org/ontology/LNC/8302-2> a <http://www.w3.org/2002/07/owl#Class>;
  <http://bioportal.bioontology.org/ontologies/umls/cui> "C0487985";
  <http://bioportal.bioontology.org/ontologies/umls/hasSTY> <http://purl.bioontology.org/ontology/STY/T201>;
  <http://bioportal.bioontology.org/ontologies/umls/tui> "T201";
  <http://purl.bioontology.org/ontology/LNC/ASSOC_OBSERVATIONS> "89263-8";
  <http://purl.bioontology.org/ontology/LNC/COMMON_ORDER_RANK> "0";
  <http://purl.bioontology.org/ontology/LNC/COMMON_SI_TEST_RANK> "0";
  <http://purl.bioontology.org/ontology/LNC/COMMON_TEST_RANK> "0";
  <http://purl.bioontology.org/ontology/LNC/CONSUMER_NAME> "Body height";
  <http://purl.bioontology.org/ontology/LNC/EXAMPLE_UCUM_UNITS> "[in_i]";
  <http://purl.bioontology.org/ontology/LNC/EXAMPLE_UNITS> "inches; cm";
  <http://purl.bioontology.org/ontology/LNC/IEEE_CF_CODE10> "68060";
  <http://purl.bioontology.org/ontology/LNC/IEEE_DIM> "L";
  <http://purl.bioontology.org/ontology/LNC/IEEE_REFID> "MDC_ATTR_PT_HEIGHT";
  <http://purl.bioontology.org/ontology/LNC/IEEE_UOM_UCUM> "cm [in_i]";
  <http://purl.bioontology.org/ontology/LNC/LCL> "BDYHGT.ATOM";
  <http://purl.bioontology.org/ontology/LNC/LCN> "2";
  <http://purl.bioontology.org/ontology/LNC/LCS> "ACTIVE";
  <http://purl.bioontology.org/ontology/LNC/LCT> "MIN";
  <http://purl.bioontology.org/ontology/LNC/LOINC_COMPONENT> "Body height";
  <http://purl.bioontology.org/ontology/LNC/LOINC_PROPERTY> "Len";
  <http://purl.bioontology.org/ontology/LNC/LOINC_SCALE_TYP> "Qn";
  <http://purl.bioontology.org/ontology/LNC/LOINC_SYSTEM> "^Patient";
  <http://purl.bioontology.org/ontology/LNC/LOINC_TIME_ASPECT> "Pt";
  <http://purl.bioontology.org/ontology/LNC/LRN2> "Axial length; bod; Bodies; BODY HEIGHT(LENGTH).ATOM; Body length; Length; Point in time; QNT; Quan; Quant; Quantitative; Random";
  <http://purl.bioontology.org/ontology/LNC/UNITS_AND_RANGE> "[in_us]:;[ft_us]:;cm:;m:;";
  <http://purl.bioontology.org/ontology/LNC/VERSION_FIRST_RELEASED> "1.0h(2)";
  <http://purl.bioontology.org/ontology/LNC/VERSION_LAST_CHANGED> "2.64";
  <http://purl.bioontology.org/ontology/LNC/has_class> <http://purl.bioontology.org/ontology/LNC/LP7768-7>;
  <http://purl.bioontology.org/ontology/LNC/has_component> <http://purl.bioontology.org/ontology/LNC/LP64598-3>;
  <http://purl.bioontology.org/ontology/LNC/has_fragments_for_synonyms> <http://purl.bioontology.org/ontology/LNC/LP28569-9>;
  <http://purl.bioontology.org/ontology/LNC/has_property> <http://purl.bioontology.org/ontology/LNC/LP6822-3>;
  <http://purl.bioontology.org/ontology/LNC/has_scale> <http://purl.bioontology.org/ontology/LNC/LP7753-9>;
  <http://purl.bioontology.org/ontology/LNC/has_supersystem> <http://purl.bioontology.org/ontology/LNC/LP6985-8>;
  <http://purl.bioontology.org/ontology/LNC/has_time_aspect> <http://purl.bioontology.org/ontology/LNC/LP6960-1>;
  <http://purl.bioontology.org/ontology/LNC/measures> <http://purl.bioontology.org/ontology/LNC/MTHU014711>;
  <http://purl.bioontology.org/ontology/LNC/member_of> <http://purl.bioontology.org/ontology/LNC/34566-0>,
    <http://purl.bioontology.org/ontology/LNC/54126-8>, <http://purl.bioontology.org/ontology/LNC/55177-0>,
    <http://purl.bioontology.org/ontology/LNC/55418-8>, <http://purl.bioontology.org/ontology/LNC/58446-6>,
    <http://purl.bioontology.org/ontology/LNC/62819-8>, <http://purl.bioontology.org/ontology/LNC/67795-5>,
    <http://purl.bioontology.org/ontology/LNC/67871-4>, <http://purl.bioontology.org/ontology/LNC/70300-9>,
    <http://purl.bioontology.org/ontology/LNC/72513-5>, <http://purl.bioontology.org/ontology/LNC/74028-2>,
    <http://purl.bioontology.org/ontology/LNC/74293-2>, <http://purl.bioontology.org/ontology/LNC/74728-7>,
    <http://purl.bioontology.org/ontology/LNC/75854-0>, <http://purl.bioontology.org/ontology/LNC/76046-2>,
    <http://purl.bioontology.org/ontology/LNC/85057-8>, <http://purl.bioontology.org/ontology/LNC/85353-1>,
    <http://purl.bioontology.org/ontology/LNC/86641-8>, <http://purl.bioontology.org/ontology/LNC/87825-6>,
    <http://purl.bioontology.org/ontology/LNC/89543-3>;
  <http://www.w3.org/2000/01/rdf-schema#subClassOf> <http://purl.bioontology.org/ontology/LNC/LP248771-0>,
    <http://purl.bioontology.org/ontology/LNC/MTHU000030>;
  <http://www.w3.org/2004/02/skos/core#altLabel> "Body height"@en, "Body height:Length:Point in time:^Patient:Quantitative"@en;
  <http://www.w3.org/2004/02/skos/core#notation> "8302-2";
  <http://www.w3.org/2004/02/skos/core#prefLabel> "Body height:Len:Pt:^Patient:Qn"@en .
```

----

```SPARQL
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX LNC: <http://purl.bioontology.org/ontology/LNC/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
select 
(str(?componentlab) as ?comp) 
(str(?propertylab) as ?prop) 
(str(?timelab) as ?time) 
(str(?systemlab) as ?sys) 
(str(?scalelab) as ?scale) 
(str(?methodlab) as ?meth) 
(str(?classlab) as ?class)
where {
    graph LNC: {
        ?LoincTerm skos:notation "2345-7" ;
                                          LNC:has_class ?LoincClass ;
                                          LNC:has_component ?LoincComponent ;
                                          LNC:has_property ?LoincProperty ;
                                          LNC:has_scale ?LoincScale ;
                                          LNC:has_system ?LoincSystem ;
                                          LNC:has_time_aspect ?LoincTime .
        optional {
            ?LoincTerm LNC:has_method ?LoincMethod .
            ?LoincMethod skos:prefLabel ?methodlab .
        }
        ?LoincComponent skos:prefLabel ?componentlab .
        ?LoincProperty skos:prefLabel ?propertylab .
        ?LoincTime skos:prefLabel ?timelab .
        ?LoincSystem skos:prefLabel ?systemlab .
        ?LoincScale skos:prefLabel ?scalelab .
        ?LoincClass skos:prefLabel ?classlab .
    }
}
```
>Showing results from 1 to 3 of 3. Query took 0.1s, moments ago.

row  | comp | prop | time | sys | scale | meth | class
-- | -- | -- | -- | -- | -- | -- | --
1 | Glucose | MCnc | Pt | Plas | Qn |   | CHEM
2 | Glucose | MCnc | Pt | Ser | Qn |   | CHEM
3 | Glucose | MCnc | Pt | Ser/Plas | Qn |   | CHEM

Compare to fields from https://search.loinc.org/searchLOINC/ (requires free LOINC account)

See also https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/LNC/sourcerepresentation.html

Search Res Field | Order | Included in SPARQL above | Notes
-- | -- | -- | --
search score | 0 | _NA_
LOINC | 1 | TRUE
LongName | 2 | **FALSE** | skos:altLabel ... but there may be multiple `skos:altLabel`s 
Component | 3 | TRUE
Property | 4 | TRUE
Timing | 5 | TRUE
System | 6 | TRUE
Scale | 7 | TRUE
Method | 8 | TRUE
exUCUMunits | 9 | **FALSE** | LNC:EXAMPLE_UCUM_UNITS
exUnits | 10 | **FALSE** | LNC:EXAMPLE_UNITS
Lforms | 11 | **FALSE** |  LNC:member_of ?
Rank | 12 | **FALSE** | LNC:COMMON_ORDER_RANK or LNC:COMMON_TEST_RANK ?
SIRank | 13 | **FALSE** | LNC:COMMON_SI_TEST_RANK
Class | 14 | TRUE | using LNC:has_class. see also LNC:class_of and rdfs:subClassOf.
ShortName | 15 | **FALSE** | skos:altLabel ... but there may be multiple `skos:altLabel`s 
DisplayName | 16 | **FALSE** | skos:prefLabel, skos:altLabel, LNC:LRN2,  LNC:LCL or LNC:CONSUMER_NAME?
Type | 17 | **FALSE** | rdf:type always owl:Class
OrderObs | 18 | **FALSE** | LNC:LOR
Copyright | 19 | **FALSE** | LNC:COPYRIGHT or LNC:EXTERNAL_COPYRIGHT_LINK ?

```SPARQL
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX LNC: <http://purl.bioontology.org/ontology/LNC/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
select 
?p (count( ?LoincTerm ) as ?count)
where {
    graph LNC: {
        ?LoincTerm ?p ?o 
    }
}
group by ?p 
order by desc (count( ?LoincTerm))
```

row  | p | count
-- | -- | --
1 | skos:altLabel | 256165
2 | umls:hasSTY | 249014
3 | umls:tui | 249014
4 | rdf:type | 217363
5 | skos:notation | 217206
6 | skos:prefLabel | 217206
7 | umls:cui | 217079
8 | LNC:component_of | 213936
9 | LNC:has_component | 213936
10 | rdfs:subClassOf | 200886
11 | LNC:has_system | 123784
12 | LNC:system_of | 123784
13 | LNC:class_of | 94933
14 | LNC:has_class | 94933
15 | LNC:fragments_for_synonyms_of | 91277
16 | LNC:has_fragments_for_synonyms | 91277
17 | LNC:COMMON_ORDER_RANK | 86217
18 | LNC:COMMON_SI_TEST_RANK | 86217
19 | LNC:COMMON_TEST_RANK | 86217
20 | LNC:LCL | 86217
21 | LNC:LCN | 86217
22 | LNC:LCS | 86217
23 | LNC:LCT | 86217
24 | LNC:LOINC_COMPONENT | 86217
25 | LNC:LOINC_PROPERTY | 86217
26 | LNC:LOINC_SCALE_TYP | 86217
27 | LNC:LOINC_SYSTEM | 86217
28 | LNC:LOINC_TIME_ASPECT | 86217
29 | LNC:LRN2 | 86217
30 | LNC:VERSION_FIRST_RELEASED | 86217
31 | LNC:VERSION_LAST_CHANGED | 86217
32 | LNC:has_time_aspect | 85880
33 | LNC:time_aspect_of | 85880
34 | LNC:has_scale | 83096
35 | LNC:scale_of | 83096
36 | LNC:has_property | 82442
37 | LNC:property_of | 82442
38 | LNC:answer_to | 78544
39 | LNC:has_answer | 78544
40 | LNC:measured_by | 75772
41 | LNC:measures | 75772
42 | LNC:LOR | 71297
43 | LNC:analyzed_by | 69485
44 | LNC:analyzes | 69485
45 | LNC:has_method | 61707
46 | LNC:method_of | 61707
47 | LNC:LUR | 55228
48 | LNC:LOINC_METHOD_TYP | 45275
49 | LNC:has_member | 40344
50 | LNC:member_of | 40344
51 | LNC:EXAMPLE_UNITS | 34650
52 | LNC:EXAMPLE_UCUM_UNITS | 34413
53 | LNC:has_suffix | 31446
54 | LNC:suffix_of | 31446
55 | LNC:CHANGE_REASON_PUBLIC | 27531
56 | LNC:OID | 23021
57 | LNC:has_multipart | 19998
58 | LNC:multipart_of | 19998
59 | LNC:challenge_of | 18620
60 | LNC:has_challenge | 18620
61 | LNC:has_supersystem | 17881
62 | LNC:supersystem_of | 17881
63 | LNC:ANSWER_LIST_ID | 17151
64 | LNC:ANSWER_LIST_NAME | 17151
65 | LNC:LSU | 14225
66 | LNC:SOS | 10444
67 | LNC:ANSWER_CODE | 9706
68 | LNC:evaluation_of | 9200
69 | LNC:has_evaluation | 9200
70 | LNC:LQS | 8661
71 | LNC:ASSOC_OBSERVATIONS | 7296
72 | LNC:LQT | 7175
73 | LNC:HL7_ATTACHMENT_STRUCTURE | 7019
74 | LNC:has_imaged_location | 6270
75 | LNC:is_imaged_location_for | 6270
76 | LNC:IMAGING_DOCUMENT_VALUE_SET | 5971
77 | LNC:has_modality_type | 5967
78 | LNC:is_modality_type_for | 5967
79 | LNC:COPYRIGHT | 5819
80 | LNC:EXTERNAL_COPYRIGHT_LINK | 5819
81 | LNC:has_imaging_focus | 5385
82 | LNC:is_imaging_focus_of | 5385
83 | LNC:divisor_of | 5299
84 | LNC:has_divisor | 5299
85 | LNC:LEA | 4759
86 | LNC:PANEL_TYPE | 3622
87 | LNC:has_lateral_anatomic_location | 3464
88 | LNC:is_lateral_anatomic_location_of | 3464
89 | LNC:has_lateral_location_presence | 3446
90 | LNC:is_presence_of_lateral_location | 3446
91 | LNC:has_aggregation_view | 3289
92 | LNC:is_aggregation_view_of | 3289
93 | LNC:has_timing_of | 2665
94 | LNC:is_timing_for | 2665
95 | LNC:associated_with | 2490
96 | LNC:has_given_pharmaceutical_substance | 2392
97 | LNC:is_given_pharmaceutical_substance_for | 2392
98 | LNC:TOP_2000_LAB_RESULTS_SI | 2192
99 | LNC:TOP_2000_LAB_RESULTS_US | 2178
100 | LNC:has_pharmaceutical_route | 1909
101 | LNC:is_pharmaceutical_route_for | 1909
102 | LNC:has_view_type | 1668
103 | LNC:is_view_type_for | 1668
104 | LNC:UNIVERSAL_LAB_ORDERS_VALUE_SET | 1522
105 | LNC:STATUS_TEXT | 1315
106 | LNC:LOINC_SCORE | 1115
107 | LNC:has_modality_subtype | 1069
108 | LNC:is_modality_subtype_for | 1069
109 | LNC:has_action_guidance | 1021
110 | LNC:is_action_guidance_for | 1021
111 | LNC:has_presence_guidance | 999
112 | LNC:is_presence_guidance_for | 999
113 | LNC:LMP | 887
114 | LNC:mapped_from | 887
115 | LNC:mapped_to | 887
116 | LNC:LCD | 824
117 | LNC:IEEE_CF_CODE10 | 631
118 | LNC:IEEE_REFID | 631
119 | LNC:LFO | 590
120 | LNC:IEEE_UOM_UCUM | 589
121 | LNC:IEEE_DIM | 558
122 | LNC:IEEE_DESCRIPTION | 488
123 | LNC:has_exam | 431
124 | LNC:is_exam_for | 431
125 | skos:definition | 389
126 | LNC:has_maneuver_type | 380
127 | LNC:is_maneuver_type_for | 380
128 | LNC:has_object_guidance | 354
129 | LNC:is_object_guidance_for | 354
130 | LNC:has_approach_guidance | 340
131 | LNC:is_approach_guidance_for | 340
132 | LNC:has_quotient | 311
133 | LNC:quotient_of | 311
134 | LNC:adjustment_of | 283
135 | LNC:has_adjustment | 283
136 | LNC:LSP | 259
137 | LNC:UNITS_AND_RANGE | 211
138 | LNC:CONSUMER_NAME | 193
139 | LNC:has_time_modifier | 191
140 | LNC:time_modifier_of | 191
141 | rdfs:comment | 157
142 | rdfs:label | 157
143 | LNC:count_of | 147
144 | LNC:has_count | 147
145 | LNC:ANSWER_CODE_SYSTEM | 119
146 | LNC:HL7_FIELD_SUBFIELD_ID | 92
147 | LNC:HL7_ATTACHMENT_REQUEST | 80
148 | LNC:SUBSEQUENT_TEXT_PROMPT | 60
149 | LNC:has_subject | 33
150 | LNC:is_subject_of | 33
151 | LNC:ASK_AT_ORDER_ENTRY | 3
152 | LNC:EXAMPLE_SI_UCUM_UNITS | 2
153 | LNC:FROMRSAB | 1
154 | LNC:FROMVSAB | 1
155 | LNC:MAPSETRSAB | 1
156 | LNC:MAPSETVERSION | 1
157 | LNC:MAPSETVSAB | 1
158 | LNC:MTH_MAPFROMCOMPLEXITY | 1
159 | LNC:MTH_MAPFROMEXHAUSTIVE | 1
160 | LNC:MTH_MAPSETCOMPLEXITY | 1
161 | LNC:MTH_MAPTOCOMPLEXITY | 1
162 | LNC:MTH_MAPTOEXHAUSTIVE | 1
163 | LNC:TORSAB | 1
164 | LNC:TOVSAB | 1
165 | LNC:mth_expanded_form_of | 1
166 | LNC:mth_has_expanded_form | 1
167 | owl:imports | 1
168 | owl:versionInfo | 1

# Creation of https://github.com/pennbiobank/turbo/blob/master/docs/private_stage/PDS_labs_LOINC_annotation_and_tabulation.csv

```SQL
SELECT
	lr.FK_Lab_Result_Item_ID,
	rlri.LAB_RESULT_ITEM_CODE,
	rlri.LAB_RESULT_ITEM_DESCRIPTION,
	rlri.LOINC,
	COUNT(1)
FROM
	mdm.LAB_RESULTS lr
JOIN mdm.R_LAB_RESULT_ITEM rlri ON
	lr.FK_LAB_RESULT_ITEM_ID = rlri.PK_LAB_RESULT_ITEM_ID
GROUP BY
	lr.FK_Lab_Result_Item_ID,
	rlri.LAB_RESULT_ITEM_CODE,
	rlri.LAB_RESULT_ITEM_DESCRIPTION,
	rlri.LOINC
ORDER BY
	COUNT(1) DESC
```

Roughly 9200 of 30.000 PDS labs have LOINC annotations, It looks like roughly 30 of those are invalid.

```R
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$prop, useNA = 'always'), decreasing = TRUE)

   MCnc   PrThr    ACnc    Prid    SCnc     NFr    Titr     Imp  ArVRat    NCnc     Arb    CCnc    MRto   Naric    Type    SRto 
   8562    7533    5065    3869    1833    1588    1431    1218     879     780     598     535     483     481     460     439 
    CFr   Ratio     MFr   LaCnc    MRat RelTime  EntNum   LnCnc     Len    NRto    Time  RelRto   Score    Find RelMCnc RelACnc 
    392     381     366     315     262     236     216     210     172     163     155     139     123     117      96      93 
     ID     Vol  EntLen     MoM   PPres    Susc     VFr    MCnt    CCnt  ArEnrg  EntVol     Num   TmStp   LsCnc   LenFr RelCCnc 
     86      84      78      71      70      69      56      52      45      44      43      42      42      41      39      38 
 Zscore    VRat EntMass    Pres     SFr   Angle    Mass    Aper   Osmol    SRat    TRto  MFr.DF     Vel    NRat    Rden      Hx 
     36      33      32      31      30      28      27      23      20      20      20      10      10       6       6       5 
   SCnt      Pn     Txt    Visc    CRat    <NA> 
      5       4       4       4       3       0 
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$time, useNA = 'always'), decreasing = TRUE)

   Pt   24H   XXX   12H    5H   48H   72H  <NA> 
39633   633   139    30     8     2     2     0 
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$sys, useNA = 'always'), decreasing = TRUE)

                      Ser                      Plas                  Ser/Plas                       Bld 
                    12088                      5420                      5361                      4655 
                    Urine                       XXX                      Tiss                  Bld/Tiss 
                     2589                      1897                      1803                      1717 
                    Stool                       PPP                       CSF                 Urine sed 
                      639                       635                       444                       376 
                 Body fld                       RBC                       Cvx                      Bone 
                      327                       263                       240                       200 
                 Bone mar                  Plas/Bld              Ser/Plas/Bld                   Isolate 
                      200                       139                       139                       133 
                 Specimen                      BldA                       PRP                   Bld.dot 
                      102                        83                        74                        70 
             Bld/Bone mar                      BldV        Reference lab test                 Amnio fld 
                       56                        53                        50                        46 
               Ventilator                  Synv fld                       Vag              Bld/Tiss/Sal 
                       44                        43                        41                        36 
                 Meconium                      BldC               Nuchal fold         Nuchal fold^fetus 
                       35                        27                        24                        24 
                  Ser+CSF                     Semen                Tiss fixed                       WBC 
                       21                        20                        20                        20 
                  Urethra            Urine+Ser/Plas                     cfDNA                Plas.cfDNA 
                       18                        18                        16                        16 
                 Calculus                       Gas                   Plr fld                     BldCo 
                       14                        13                        13                        12 
                {Setting}                      Dose                   Cvx/Vag                       Cvm 
                       11                        11                        10                         9 
                 Inhl gas              Ser/Plas+CSF            Current sample                       Eye 
                        9                         9                         8                         8 
              Respiratory        Respiratory system                     Liver                       Nph 
                        7                         7                         6                         6 
                   BldCoA                    BldCoV                Hereditary               Periton fld 
                        5                         5                         5                         5 
                     Cnjt                  Exhl gas                      Hair                  Provider 
                        4                         4                         4                         4 
                   Report                      Anal                  Facility                  NBS card 
                        4                         3                         3                         3 
                Ser^donor                    Sputum           Arterial system              Pericard fld 
                        3                         3                         2                         2 
                    Retic                    Saliva                  Tiss.FNA          Ventilator.alarm 
                        2                         2                         2                         2 
                Bld^donor                       CVS          Isolate/Specimen    Oxygen delivery system 
                        1                         1                         1                         1 
Respiratory system.airway                      <NA> 
                        1                         0 
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$scale, useNA = 'always'), decreasing = TRUE)

   Qn   Ord   Nom   Nar   Doc OrdQn  <NA> 
26461  8220  3985  1669    57    55     0 
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$meth, useNA = 'always'), decreasing = TRUE)

                                   <NA>                                  Molgen                                      IA 
                                  20806                                    4811                                    2388 
                                  Probe                               Probe.amp                           Probe.amp.tar 
                                   1358                                    1358                                    1272 
                                     IF         Creatinine-based formula (MDRD)                         Electrophoresis 
                                   1135                                     729                                     473 
                                 Manual                            Manual count                              Microscopy 
                                    408                                     408                                     395 
                       Microscopy.light                                    Coag                                      IB 
                                    395                                     345                                     337 
                    Thromboelastography                               Automated                         Automated count 
                                    301                                     298                                     273 
                                   FISH                              Calculated                                 Confirm 
                                    255                                     215                                     195 
                                 Screen                                    Neut                        Immune diffusion 
                                    127                                     114                                      95 
                      Non-probe.amp.tar                           Probe.amp.sig                    Microscopy.light.HPF 
                                     86                                      86                                      82 
                        High resolution                                Comp fix                          Flow cytometry 
                                     81                                      78                                      73 
                   Microscopy.light.LPF                              Test strip                        High sensitivity 
                                     72                                      63                                      60 
                                 Chromo                                 Culture           Detection limit <= 0.01 ng/mL 
                                     58                                      54                                      54 
   Dosage of chromosome specific cf DNA                       Computer assisted                               Line blot 
                                     48                                      40                                      40 
          Reaction: lactate to pyruvate                    Calculated.FibroSure                          Immunofixation 
                                     38                                      36                                      36 
                            With P-5'-P               Organism specific culture                              Genotyping 
                                     36                                      33                                      31 
     Creatinine-based formula (CKD-EPI)                                     MIC                                     XXX 
                                     30                                      30                                      28 
                              XXX stain                                Dialysis                                    Aggr 
                                     28                                      27                                      25 
                           Immune stain                           Platelet aggr                            Direct assay 
                                     25                                      25                                      24 
                              Multidisk                              Cyto stain                                 Elution 
                                     24                                      23                                      18 
                                    RPR                         Solubility test                                    VDRL 
                                     18                                      18                                      18 
                                   Aggl                              Immunobead                     Sudan black B stain 
                                     17                                      16                                      16 
                                   HPLC                          Agar diffusion                   Crithidia luciliae IF 
                                     14                                      13                                      13 
Calculated from oxygen partial pressure                    Isoelectric focusing                                      LA 
                                     12                                      12                                      12 
                               Measured                                      US                             US.measured 
                                     12                                      12                                      12 
                            US+Measured              Detection limit <= 20 mg/L                                      MG 
                                     12                                      11                                      11 
                     3D cold incubation                     Pathologist comment                    Test strip.automated 
                                     10                                      10                                      10 
                        Trichrome stain                    Cyto stain.thin prep             Esterase stain.non-specific 
                                     10                                       9                                       8 
                                    ISE                               Pathology              Periodic acid-Schiff stain 
                                      8                                       8                                       8 
                                    RIA                                     BCG                       C1q binding assay 
                                      8                                       6                                       6 
          Detection limit <= 0.05 mIU/L                                      HA                          Immediate spin 
                                      6                                       6                                       6 
                  Myeloperoxidase stain                         Raji cell assay                              Westergren 
                                      6                                       6                                       6 
                          Concentration                            Nephelometry                               Estimated 
                                      5                                       5                                       4 
                         Estimated from      Estimated from glycated hemoglobin                              Gram stain 
                                      4                                       4                                       4 
                        Kleihauer-Betke                             Phenotyping                              Rees-Ecker 
                                      4                                       4                                       4 
                                 {Role}                         Aerobic culture                     Bacterial subtyping 
                                      3                                       3                                       3 
                         Biopsy culture                              Glucometer                            Wright stain 
                                      3                                       3                                       3 
                    37 deg C incubation                  3D 37 deg C incubation                         Acid fast stain 
                                      2                                       2                                       2 
                         Hep2 substrate                           Malaria smear           Reaction: pyruvate to lactate 
                                      2                                       2                                       2 
                                  Smear                                     SSA                                Estimate 
                                      2                                       2                                       1 
                               Oximetry                            Rosette test                       Screen>1000 ng/mL 
                                      1                                       1                                       1 
                        Screen>50 ng/mL                         Wet preparation 
                                      1                                       1 
> sort(table(PDS_labs_promising_LOINC_annotated_only_and_tabulation_report$class, useNA = 'always'), decreasing = TRUE)

             CHEM             MICRO           MOLPATH          DRUG/TOX            HEM/BC              SERO           ALLERGY 
            11522              8832              2668              2476              2412              2003              1803 
      MOLPATH.MUT              COAG              CHAL                UA          CELLMARK      CHAL.ROUTINE              SPEC 
             1744              1411               977               906               521               438               432 
   MOLPATH.TRNLOC             BLDBK               HLA              MISC           ABXBACT MOLPATH.NUCREPEAT              PATH 
              396               381               367               203               167               138                94 
     HL7.GENETICS           TUMRRGT              PULM             OB.US       MOLPATH.DEL              FERT              CYTO 
               93                90                62                48                47                36                32 
      HL7.CYTOGEN      MOLPATH.MISC   MOLPATH.TRISOMY              CLIN     HEMODYN.MOLEC       MOLPATH.INV          DRUGDOSE 
               29                27                24                18                14                12                 8 
           H&P.HX      DOC.ONTOLOGY            ATTACH    ATTACH.CLINRPT           BP.ATOM        PANEL.CHEM              <NA> 
                5                 3                 2                 2                 2                 2                 0 
```
