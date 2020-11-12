# turbo-assay-templating

Here's a first step in modelling assays by using terms from [OBO foundry ontologies](http://www.obofoundry.org/). It generates a template that can be executed by the [ROBOT](http://robot.obolibrary.org/) application, leading to an RDF data model. Assays currently in scope have a specimen collected from an organism as the evaluant and can be used in EHR secondary research or in vivo basic research. 

Where we say "assays", healthcare providers might say "lab tests". However, we do not mean to imply that these tests are only relevant in the course of diagnosing disease or planning and monitoring therapy. They would be equally applicable to modelling in vivo research.

When this repo was created, we had been bootstrapping the assay models from LOINC. See details below.  

**After closer review of the LOINC terms and conditions, we have decided to put the previous approach on hold** and switch to an NLP + manual review approach for matching lab results from the EHR to terms from OBO foundry ontologies. See file XXX.



----



The LOINC-based approach required that

- The user has a LOINC account
- https://loinc.org/download/loinc-and-relma-complete-download-file/ had been downloaded. Our experience was mostly with version 2.68.
- The downloaded archive was unpacked (as a folder called `LOINC_2`) into `pipeline/data/`

Then the template could the be created and executed like this:

```Bash
cd pipeline
Rscript turbo_assay_template_generation.R
./execute_template.sh
```

There are several dependencies which are not discussed here yet.

