library(devtools)

# requires a properly formatted "turbo_R_setup.yaml" in the home directory of the user who started this script
# see https://gist.github.com/turbomam/a3915d00ee55d07510493a9944f96696 for template
source_gist(id = "https://gist.github.com/turbomam/f082295aafb95e71d109d15ca4535e46",
            sha1 = "dbc656aaf63b23dfdd35d875f6772e7c468170a4",
            filename = "turbo_R_setup.R")

library(tidyverse)

# clarify requirement for LOINC folder in location XXX

# # monitor
# ehr.bad.loincs
# ehr.na.loincs

# deal with multi-gene assays
# gene.check

# assays that are dropped because they include a numerator or a normalization

# ehr results or loincs with too few unique EMPIs
# we're not currently counting unique EMPIS per LOINC

# parts that are silently (?) dropped or less-annotated because they have multiple paths to root

####

# looks like some of the mapping work isn't going in a usefully sorted order
# could do cosine or jaccard between query and ols label, synonyms and truncatees all concatenated together

# expect 'reviewed part frame' from staged.loinc.mapping to be empty unless passed as an argument
# but start collecting mappings for all parts into one file
# bp mappings and ols search results should have
#  loinc code and label
#  external (obo?) id (IRI?) and label
#  hierarchy if available
#  mapping quality metadata

# staged.loinc.mapping variable shouldn't contain the word 'system'... they are part agnostic


# > library(tidyverse)
# ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
# ✓ tibble  3.0.1     ✓ purrr   0.3.4
# ✓ tidyr   1.1.0     ✓ forcats 0.5.0
# ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# x NLP::annotate()         masks ggplot2::annotate()
# x tibble::as_data_frame() masks igraph::as_data_frame(), dplyr::as_data_frame()
# x randomForest::combine() masks dplyr::combine()
# x purrr::compose()        masks igraph::compose()
# x NLP::content()          masks httr::content()
# x tidyr::crossing()       masks igraph::crossing()
# x dplyr::filter()         masks stats::filter()
# x purrr::flatten()        masks jsonlite::flatten()
# x igraph::groups()        masks dplyr::groups()
# x dplyr::lag()            masks stats::lag()
# x purrr::lift()           masks caret::lift()
# x randomForest::margin()  masks ggplot2::margin()
# x caret::progress()       masks httr::progress()
# x purrr::simplify()       masks igraph::simplify()

# setwd("~/loinc_in_obo/Loinc_2/")
# setwd("~/loinc_in_obo/")

# VPN and tunnel may be required
# set that up outside of this script

repeat.ehr.loinc.query <- FALSE
min.empis <- 2

PartRelatedCodeMapping.fn <-
  "Loinc_2/AccessoryFiles/PartFile/PartRelatedCodeMapping.csv"

hor.path <- "/Users/markampa/loinc_in_obo/"

bioontologies.base <- "http://data.bioontology.org/ontologies/"

obo.base <- "http://purl.obolibrary.org/obo/"

bp.loinc.prefix <- "http://purl.bioontology.org/ontology/LNC/"

SYSTEM_ols_search_res_reviewed.fn <-
  "SYSTEM_ols_search_res_reviewed.csv"

reviewed.part.filename = "SYSTEM_to_specimen_labels.csv"

Part.min.cols <-
  c("PartNumber", "PartName", "PartDisplayName", "Status")

excluded.Properties <- c("http://loinc.org/property/category")

####

MultiAxialHierarchy <-
  read_csv("Loinc_2/AccessoryFiles/MultiAxialHierarchy/MultiAxialHierarchy.csv")

LoincPartLink_Primary <-
  read_csv("Loinc_2/AccessoryFiles/PartFile/LoincPartLink_Primary.csv")
LoincPartLink_Supplementary <-
  read_csv("Loinc_2/AccessoryFiles/PartFile/LoincPartLink_Supplementary.csv")
LoincPartLink <-
  rbind.data.frame(LoincPartLink_Primary, LoincPartLink_Supplementary)

Part <- read_csv("Loinc_2/AccessoryFiles/PartFile/Part.csv")

Loinc <- read.csv("Loinc_2/LoincTable/Loinc.csv")

# use a single frame for systems, components and methods
reviewed.part.frame <-
  read_csv("~/loinc_in_obo/SYSTEM_to_specimen_labels.csv")
print(head(reviewed.part.frame))

PartRelatedCodeMapping <-
  read_csv(PartRelatedCodeMapping.fn)

SYSTEM_ols_search_res_reviewed <-
  read_csv(SYSTEM_ols_search_res_reviewed.fn)

Part.min <- Part[, Part.min.cols]

harmonized_ontology_rankings <-
  read.csv(paste0(hor.path, config$harmonized_ontology_rankings),
           stringsAsFactors = FALSE)

# say rank, not priority
ontology.rankings.min <-
  cbind.data.frame(
    ontology.iri = paste(bioontologies.base,
                         harmonized_ontology_rankings$rhs,
                         sep = ""),
    ontology.rank  = harmonized_ontology_rankings$priority
  )

onto.list.for.ols <-
  paste0(sort(tolower(harmonized_ontology_rankings$rhs)), collapse = ',')

ehr.driver <-
  JDBC(driverClass = "oracle.jdbc.OracleDriver",
       classPath = config$oracle.jdbc.path)

ehr.con.string <- paste0(
  "jdbc:oracle:thin:@//",
  config$pds.host,
  ":",
  config$pds.port,
  "/",
  config$pds.database
)

ehr.connection <- NULL

####

ehr.test.and.refresh <- function() {
  local.q <- "select 1 from dual"
  tryCatch({
    dbGetQuery(ehr.connection, local.q)
  }, warning = function(w) {
    
  }, error = function(e) {
    print(e)
    print("trying to reconnect")
    ehr.connection <<-
      dbConnect(ehr.driver,
                ehr.con.string,
                config$pds.user,
                config$pds.pw)
    dbGetQuery(ehr.connection, local.q)
  }, finally = {
    
  })
}

loinc.strip.url.lhs <- function(loinc.urls, lhs) {
  temp <-
    make.names(
      names = sub(pattern = lhs,
                  replacement = "",
                  loinc.urls),
      unique = TRUE,
      allow_ = TRUE
    )
  return(temp)
}

lpl.wide.builder <-
  function(names.from,
           values.from,
           modify.colnames = TRUE,
           to.strip = "http://loinc.org/property/") {
    temp <-
      unique(LoincPartLink[LoincPartLink$Property %in% accepted.properties &
                             LoincPartLink$LoincNumber %in% accepted.assays,
                           c("LoincNumber", names.from, values.from)])
    
    temp <-
      temp %>% pivot_wider(names_from = names.from  ,
                           values_from = values.from)
    
    if (modify.colnames) {
      colnames(temp) <-
        loinc.strip.url.lhs(colnames(temp), to.strip)
    }
    return(temp)
  }

get.jaccard.from.tibble.cols <-
  function(current_tibble,
           current_colname,
           current_comparator) {
    print(current_colname)
    print(current_comparator)
    current_col.contents <-
      pull(current_tibble, current_colname)
    comparator.col.contents <-
      pull(current_tibble, current_comparator)
    current_jaccard.num <-
      length(intersect(current_col.contents, comparator.col.contents))
    current_jaccard.denom <-
      length(union(current_col.contents, comparator.col.contents))
    current_jaccard <-
      current_jaccard.num / current_jaccard.denom
    print(current_jaccard)
    cat("\n")
    return(list(a = current_colname,
                b = current_comparator,
                jaccard = current_jaccard))
  }

do.jaccard.calcs <- function(outer.tibble, first.col) {
  relevant.colnames <-
    sort(colnames(outer.tibble[first.col:ncol(outer.tibble)]))
  print(relevant.colnames)
  get.comparable.cols <-
    lapply(relevant.colnames, function(current_colname) {
      for.comparing <-
        relevant.colnames[relevant.colnames > current_colname]
      do.jaccard.prep <-
        lapply(for.comparing, function(current_comparator) {
          print(current_colname)
          print(current_comparator)
          jaccard.from.tibble.cols <-
            get.jaccard.from.tibble.cols(
              current_tibble = outer.tibble,
              current_colname = current_colname,
              current_comparator = current_comparator
            )
          return(jaccard.from.tibble.cols)
        })
      do.jaccard.prep <- do.call(rbind.data.frame, do.jaccard.prep)
      return(do.jaccard.prep)
    })
  comparable.cols <- do.call(rbind.data.frame, get.comparable.cols)
  comparable.cols <- comparable.cols[comparable.cols$jaccard > 0 , ]
  return(comparable.cols)
}

get.bp.labels <-
  function(mappings.frame,
           uri.col.name = "target.term",
           onto.abbrev.col.name = "target.prefix",
           api.family = "classes") {
    # mappings.frame <- SYSTEM.bp.mapres
    outer <-
      apply(
        X = mappings.frame,
        MARGIN = 1,
        FUN = function(current_row) {
          # current_row <- SYSTEM.bp.mapres[1,]
          current_full.uri <- current_row[[uri.col.name]]
          print(current_full.uri)
          encoded.term <-
            URLencode(current_full.uri, reserved = TRUE)
          # print(encoded.term)
          api.ontology.name <- current_row[[onto.abbrev.col.name]]
          # print(api.ontology.name)
          
          prepared.get <-
            paste(api.base.uri,
                  api.ontology.name,
                  api.family,
                  encoded.term,
                  sep = "/")
          
          print(prepared.get)
          
          mapping.res.list <-
            httr::GET(url = prepared.get,
                      add_headers(
                        Authorization = paste0("apikey token=", config$public.bioportal.api.key)
                      ))
          
          # print(mapping.res.list$status_code)
          
          if (mapping.res.list$status_code == 200) {
            # test for presence of content, validity of Char-ified JSON?
            retrieved.label <-
              fromJSON(rawToChar(mapping.res.list$content))$prefLabel
            print(retrieved.label)
            if (!(is.null(label = retrieved.label))) {
              return(list(target.term = current_full.uri, label = retrieved.label))
            }
            
          }
        }
      )
  }

vals.by.relation <- function(current_row, relation) {
  values <-
    strsplit(current_row[[relation]], split = "|", fixed = TRUE)
  values <- unlist(values)
  lp.vec <-
    rep(current_row[['loinc.part']], length(values))
  query.vec <-
    rep(current_row[['query']], length(values))
  obo.id.vec <-
    rep(current_row[['obo_id']], length(values))
  relation.vec <- rep(relation, length(values))
  values <-
    cbind.data.frame(
      loinc.part = lp.vec,
      query = query.vec,
      obo_id = obo.id.vec,
      relation = relation.vec,
      value = values
    )
  return(values)
}

# what about nodes that go though more than one path?
leaves.to.paths <- function(current.leaves) {
  # current.leaves <- needs.mapping
  
  current.leaves <-
    MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% current.leaves, ]
  
  current.leaves.ancestors <-
    lapply(sort(unique(current.leaves$CODE)), function(current.code) {
      # current.code <- "LP14304-7"
      print(current.code)
      current.leafs.paths <-
        current.leaves[current.leaves$CODE == current.code , ]
      
      if (nrow(current.leafs.paths) > 0) {
        current.leafs.paths <-
          current.leafs.paths[order(current.leafs.paths$SEQUENCE), ]
        current.leafs.ancestors <- apply(
          X = current.leafs.paths,
          MARGIN = 1,
          FUN = function(current_row) {
            # current_row <-
            current.path <- current_row[["PATH_TO_ROOT"]]
            current.sequence <- current_row[["SEQUENCE"]]
            current.nodes <-
              unlist(strsplit(current.path, ".", fixed = TRUE))
            current.frame <-
              cbind.data.frame(
                CODE = rep(current.code, length(current.nodes)),
                SEQUENCE = rep(current.sequence, length(current.nodes)),
                ancestor = current.nodes,
                position = 1:length(current.nodes)
              )
            return(current.frame)
          }
        )
        current.leafs.ancestors <-
          do.call(rbind.data.frame, current.leafs.ancestors)
        return(current.leafs.ancestors)
      }
    })
  
  if (length(current.leaves.ancestors) > 0) {
    current.leaves.ancestors <-
      do.call(rbind.data.frame, current.leaves.ancestors)
    
    current.leaves.ancestors <-
      left_join(x = current.leaves.ancestors,
                y = Part,
                by = c("ancestor" = "PartNumber"))
    
    unnamed.ancestors <-
      sort(unique(current.leaves.ancestors$ancestor[is.na(current.leaves.ancestors$PartDisplayName)]))
    
    named.ancestors <-
      current.leaves.ancestors[(!(current.leaves.ancestors$ancestor %in% unnamed.ancestors)) ,]
    
    unnamed.ancestors <-
      current.leaves.ancestors[current.leaves.ancestors$ancestor %in% unnamed.ancestors , c("CODE", "SEQUENCE", "ancestor", "position")]
    
    MultiAxialHierarchy.salvage <-
      unique(MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% unnamed.ancestors$ancestor ,
                                 c("CODE", "CODE_TEXT")])
    
    unnamed.ancestors  <-
      left_join(x = unnamed.ancestors,
                y = MultiAxialHierarchy.salvage,
                by = c("ancestor" = "CODE"))
    
    names(unnamed.ancestors) <-
      c("CODE",
        "SEQUENCE",
        "ancestor",
        "position",
        "PartDisplayName")
    current.leaves.ancestors <-
      plyr::rbind.fill(named.ancestors , unnamed.ancestors)
    
    muti.seq.check <-
      unique(current.leaves.ancestors[, c("CODE", "SEQUENCE")])
    muti.seq.check <-
      muti.seq.check[(!(is.na(muti.seq.check$SEQUENCE))), ]
    muti.seq.check <- make.table.frame(muti.seq.check$CODE)
    
    single.seq.leaves <-
      muti.seq.check$value[muti.seq.check$count == 1]
    single.seq.leaves <-
      current.leaves.ancestors[current.leaves.ancestors$CODE %in% single.seq.leaves , ]
    
    # what about nodes that go though more than one path?
    
    labeled.paths <-
      lapply(sort(unique(single.seq.leaves$CODE)), function(current_code) {
        current_path <-
          single.seq.leaves[single.seq.leaves$CODE == current_code ,]
        current_path <-
          current_path[order(current_path$position) ,]
        current_path <-
          paste(current_path$PartDisplayName, collapse = ">>>")
        # print(current_path)
        return(list(CODE = current_code, path = current_path))
      })
    labeled.paths <- do.call(rbind.data.frame, labeled.paths)
    return(labeled.paths)
  }
}

# document where these defaults come from
# move them into a config file
get.common.constrained.assays <-
  function(all.assays = max.one.gene.wide,
           component.code.list = c(),
           method.code.list = c(),
           property.name.list = c(
             "ACnc",
             "Arb",
             "ArVRat",
             "CCnc",
             "CFr",
             "LaCnc",
             "LnCnc",
             "LsCnc",
             "MCnc",
             "MFr",
             "MFr.DF",
             "MRto",
             "NCnc",
             "NFr",
             "NRto",
             "Num",
             "PPres",
             "PrThr",
             "Ratio",
             "RelACnc",
             "RelCCnc",
             "RelMCnc",
             "RelRto",
             "SatFr",
             "SCnc",
             "SCnt",
             "SFr",
             "SRto",
             "VFr"
           ),
           scale.name.list = c("Qn", "Ord", "Nom"),
           system.name.list = c(
             "Amnio fld",
             "Amnio fld/CVS",
             "Bld",
             "Bld/Bone mar",
             "BldA",
             "BldC",
             "BldCo",
             "BldV",
             "Body fld",
             "Bone mar",
             "CSF",
             "Cvm",
             "CVS",
             "Cvx",
             "Cvx/Vag",
             "Eye",
             "Liver",
             "Meconium",
             "Nph",
             "Pericard fld",
             "Periton fld",
             "Plas",
             "Plr fld",
             "RBC",
             "Respiratory system",
             "Respiratory.upper",
             "Retic",
             "Saliva",
             "Semen",
             "Ser",
             "Ser/Plas",
             "Ser/Plas/Bld",
             "Ser+Bld",
             "Ser+CSF",
             "Specimen",
             "Sputum",
             "Stool",
             "Synv fld",
             "Thrt",
             "Urethra",
             "Urine",
             "Vag",
             "WBC"
           ),
           time.name.list = c("Pt", "24H"),
           min.min.patients = 2,
           drop.na.cols = TRUE,
           useless.columns = c(
             "EXTERNAL_COPYRIGHT_LINK",
             "EXTERNAL_COPYRIGHT_NOTICE",
             "HL7_ATTACHMENT_STRUCTURE",
             "HL7_FIELD_SUBFIELD_ID",
             "PanelType",
             "ValidHL7AttachmentRequest",
             "analyte.divisor.suffix",
             "PartDisplayName.analyte.divisor.suffix",
             "PartName.analyte.divisor.suffix",
             "Status.analyte.divisor.suffix",
             "analyte.numerator",
             "PartDisplayName.analyte.numerator",
             "PartName.analyte.numerator",
             "Status.analyte.numerator",
             "adjustment",
             "PartDisplayName.adjustment",
             "PartName.adjustment",
             "Status.adjustment"
           ),
           drops.assay = c("analyte.numerator", "adjustment")) {
    common.constrained.assays <- all.assays
    if (length(component.code.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$analyte.core %in% component.code.list  |
                                    max.one.gene.wide$COMPONENT %in%  component.code.list , ]
    }
    if (length(method.code.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$METHOD_TYP %in% method.code.list ,]
    }
    if (length(property.name.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$PartName.PROPERTY %in% property.name.list ,]
    }
    if (length(scale.name.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$PartName.SCALE_TYP %in% scale.name.list ,]
    }
    if (length(system.name.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$PartName.SYSTEM %in% system.name.list ,]
    }
    if (length(time.name.list) > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$PartName.TIME_ASPCT %in% time.name.list ,]
    }
    if (is.numeric(min.min.patients) & min.min.patients > 0) {
      common.constrained.assays <-
        common.constrained.assays[common.constrained.assays$UNIQUE_EMPI_COUNT >= min.min.patients ,]
    }
    if (drop.na.cols) {
      all.na.flag <-
        sapply(common.constrained.assays, function(x)
          all(is.na(x)))
      print(table(all.na.flag))
      common.constrained.assays <-
        common.constrained.assays[, (!(all.na.flag))]
    }
    # drop assays first!
    if (length(drops.assay) > 0) {
      dropper.flag <- as.matrix(common.constrained.assays[, drops.assay])
      dropper.flag[dropper.flag == ""] <- NA
      dropper.flag <- (!(is.na(dropper.flag)))
      dropper.flag <- rowSums(dropper.flag)
      dropper.flag <- dropper.flag > 0
      print(table(dropper.flag))
      common.constrained.assays <-
        common.constrained.assays[(!(dropper.flag)), ]
    }
    if (length(useless.columns) > 0) {
      common.constrained.assays <-
        common.constrained.assays[, setdiff(colnames(common.constrained.assays), useless.columns)]
    }
    return(common.constrained.assays)
  }

get.column.usefulness <- function(data) {
  column.usefulness <-
    lapply(sort(colnames(data)), function(current_colname) {
      # print(current_colname)
      row_count <- nrow(data)
      unique.vals <-
        nrow(unique(data[, current_colname]))
      unique.fraction <- unique.vals / row_count
      na.count <- is.na(data[, current_colname])
      na.count <- sum(as.numeric(na.count))
      na.fraction <- na.count / row_count
      blank.count <-
        data[, current_colname] == ""
      blank.count[is.na(blank.count)] <- 0
      blank.count <- sum(as.numeric(blank.count))
      blank.fraction <- blank.count / row_count
      useless.fraction <- na.fraction + blank.fraction
      return(
        list(
          current_colname = current_colname,
          unique.vals = unique.vals,
          unique.fraction = unique.fraction,
          na.count = na.count,
          na.fraction = na.fraction,
          blank.count = blank.count,
          blank.fraction = blank.fraction,
          useless.fraction = useless.fraction
        )
      )
    })
  column.usefulness <-
    do.call(rbind.data.frame, column.usefulness)
  return(column.usefulness)
}

chebi.sys.iri <- "https://www.ebi.ac.uk/chebi"
ncbitaxon.sys.iri <- "https://www.ncbi.nlm.nih.gov/taxonomy"

staged.loinc.mapping <-
  function(needs.mapping,
           reviewed.part.filename = "",
           primary.ExtCodeSystems = c(chebi.sys.iri,
                                      ncbitaxon.sys.iri),
           secondary.ExtCodeSystems = c(
             "http://www.nlm.nih.gov/research/umls/rxnorm",
             "https://www.ncbi.nlm.nih.gov/gene",
             "http://pubchem.ncbi.nlm.nih.gov",
             "https://www.ncbi.nlm.nih.gov/clinvar"
           ),
           tertiary.ExtCodeSystems = c("http://snomed.info/sct"),
           skip.ExtCodeSystems = c("http://www.radlex.org"),
           already.reviewed = c(),
           ontology.rankings,
           current.code.col.name ,
           current.name.col.name) {
    # remember, tissues and cells can't bear an analyte role
    # this may be most directly applicable to COMPONENT subparts
    # and maybe METHOD?
    
    paths.from.leaves <- leaves.to.paths(needs.mapping)
    
    # test for file existence?
    # what if the filename is NA, not just blank
    # test if frame has at least one column and the expected columns
    
    # if (nchar(reviewed.part.filename) > 0) {
    #   loinc_to_obo_mapping_reviewed <-
    #     read_csv("loinc_to_obo_mapping_reviewed.csv")
    #   needs.mapping <-
    #     sort(setdiff(needs.mapping, loinc_to_obo_mapping_reviewed$PartNumber))
    #
    # }
    
    if (length(already.reviewed) > 0) {
      needs.mapping <-
        sort(setdiff(needs.mapping, already.reviewed))
    }
    
    # now check for LOINC provided mappings into acceptable
    loinc.provided.primary <-
      PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% needs.mapping &
                               PartRelatedCodeMapping$ExtCodeSystem %in% primary.ExtCodeSystems , ]
    
    # what if there are multiple mappings from different (primary) target ontologies
    needs.mapping <-
      sort(setdiff(needs.mapping, loinc.provided.primary$PartNumber))
    
    if (nrow(loinc.provided.primary) > 0) {
      loinc.provided.primary <-
        left_join(
          x = loinc.provided.primary,
          y = Part,
          by = c("PartNumber" = "PartNumber")
        )
      loinc.provided.primary <-
        left_join(x = loinc.provided.primary,
                  y = paths.from.leaves,
                  by = c("PartNumber" = "CODE"))
      
      # generalize this
      loinc.provided.primary$ontology.abbreviation <- NA
      loinc.provided.primary$ontology.abbreviation[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri] <-
        "CHEBI"
      loinc.provided.primary$ontology.abbreviation[loinc.provided.primary$ExtCodeSystem == ncbitaxon.sys.iri] <-
        "NCBITAXON"
      
      loinc.provided.primary <-
        left_join(
          x = loinc.provided.primary,
          y = harmonized_ontology_rankings,
          by = c("ontology.abbreviation" = "rhs")
        )
      
      loinc.provided.primary$target.term <- NA
      loinc.provided.primary$target.term[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri] <-
        loinc.provided.primary$ExtCodeId[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri]
      loinc.provided.primary$target.term[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri] <-
        sub(
          pattern = ":",
          replacement = "_",
          x = loinc.provided.primary$ExtCodeId[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri],
          fixed = TRUE
        )
      loinc.provided.primary$target.term[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri] <-
        paste("http://purl.obolibrary.org/obo/",
              loinc.provided.primary$target.term[loinc.provided.primary$ExtCodeSystem == chebi.sys.iri],
              sep = "")
      loinc.provided.primary$target.term[loinc.provided.primary$ExtCodeSystem == ncbitaxon.sys.iri] <-
        paste(
          "http://purl.obolibrary.org/obo/NCBITaxon_",
          loinc.provided.primary$ExtCodeId[loinc.provided.primary$ExtCodeSystem == ncbitaxon.sys.iri],
          sep = ""
        )
      
      loinc.provided.primary <-
        loinc.provided.primary[, c(
          "PartNumber",
          "PartDisplayName",
          "target.term",
          "ExtCodeDisplayName",
          "domain",
          "path"
        )]
      
      colnames(loinc.provided.primary) <- c(
        "LOINC_PartNumber",
        "PartDisplayName",
        "target.term",
        "target.label",
        "target.domain",
        "path.to.LOINC.root"
      )
      
      loinc.provided.primary$mapping.method <- "LOINC provided"
    }
    
    loinc.provided.secondary <-
      PartRelatedCodeMapping[PartRelatedCodeMapping$PartNumber %in% needs.mapping &
                               PartRelatedCodeMapping$ExtCodeSystem %in% secondary.ExtCodeSystems , ]
    
    # print(head(loinc.provided.secondary))
    
    # don't pursue tertiary or "skip" now
    
    ####
    
    # now get BioPortal mappings and constrain to top priority targets
    # takes codes, not IRIs as input... that should be generalized
    api.ontology.name <- "LOINC"
    
    current.bp.mapres <- bp.map.retreive.and.parse(needs.mapping)
    current.bp.mapres <-
      do.call(rbind.data.frame, current.bp.mapres)
    
    # end of unfiltered
    
    if (nrow(current.bp.mapres) > 0) {
      current.bp.mapres <-
        unique(current.bp.mapres[current.bp.mapres$mapping.methods == "LOOM" &
                                   current.bp.mapres$target.ontology %in% ontology.rankings.min$ontology.iri ,
                                 c("source.term", "target.term")])
      
      current.bp.mapres$source.bare <-
        sub(
          pattern = bp.loinc.prefix,
          replacement = "",
          x = current.bp.mapres$source.term
        )
      
      ####
      
      current.bp.mapres$target.pattern <- NA
      current.bp.mapres$rhs <- NA
      current.bp.mapres$reassembled <- NA
      
      ####
      
      obo.base <- "http://purl.obolibrary.org/obo/"
      obo.flag <-
        grepl(pattern = obo.base, x = current.bp.mapres$target.term)
      current.bp.mapres$target.pattern[obo.flag] <- "OBO"
      current.bp.mapres$rhs[obo.flag] <-
        current.bp.mapres$target.term[obo.flag]
      current.bp.mapres$rhs[obo.flag] <-
        sub(pattern = obo.base,
            replacement = "",
            x = current.bp.mapres$rhs[obo.flag])
      current.bp.mapres$reassembled[obo.flag] <-
        current.bp.mapres$target.term[obo.flag]
      
      # add more patterns, and parse out IDs, too
      bp.ncbitaxon.base <-
        "http://purl.bioontology.org/ontology/NCBITAXON/"
      bp.ncbitaxon.flag <-
        grepl(pattern = bp.ncbitaxon.base, x = current.bp.mapres$target.term)
      current.bp.mapres$target.pattern[bp.ncbitaxon.flag] <-
        "BP NCBITAXON"
      current.bp.mapres$rhs[bp.ncbitaxon.flag] <-
        current.bp.mapres$target.term[bp.ncbitaxon.flag]
      current.bp.mapres$rhs[bp.ncbitaxon.flag] <-
        sub(pattern = bp.ncbitaxon.base,
            replacement = "",
            x = current.bp.mapres$rhs[bp.ncbitaxon.flag])
      current.bp.mapres$reassembled[bp.ncbitaxon.flag] <-
        paste(obo.base, "NCBITaxon_", current.bp.mapres$rhs[bp.ncbitaxon.flag], sep = "")
      
      bp.fma.base <-
        "http://purl.org/sig/ont/fma/fma"
      bp.fma.flag <-
        grepl(pattern = bp.fma.base, x = current.bp.mapres$target.term)
      current.bp.mapres$target.pattern[bp.fma.flag] <-
        "BP FMA"
      current.bp.mapres$rhs[bp.fma.flag] <-
        current.bp.mapres$target.term[bp.fma.flag]
      current.bp.mapres$rhs[bp.fma.flag] <-
        sub(pattern = bp.fma.base,
            replacement = "",
            x = current.bp.mapres$rhs[bp.fma.flag])
      current.bp.mapres$reassembled[bp.fma.flag] <-
        paste(obo.base, "FMA_", current.bp.mapres$rhs[bp.fma.flag], sep = "")
      
      ####
      
      # OBO prefixes
      current.bp.mapres$target.prefix <-
        sub(
          pattern = obo.base,
          replacement = "",
          x = current.bp.mapres$target.term
        )
      current.bp.mapres$target.prefix <-
        sub(
          pattern = "_.*$",
          replacement = "",
          x = current.bp.mapres$target.prefix
        )
      
      # filter out target IRIs that didn't match any of those patterns
      current.bp.mapres <-
        current.bp.mapres[(!(
          grepl(pattern = "^http", x = current.bp.mapres$target.prefix)
        )),]
      
      current.bp.mapres$target.prefix.uc <-
        toupper(current.bp.mapres$target.prefix)
      
      current.bp.mapres <-
        left_join(
          x = current.bp.mapres,
          y = harmonized_ontology_rankings,
          by = c("target.prefix.uc" = "rhs")
        )
      
      current.bp.mapres <-
        current.bp.mapres[order(current.bp.mapres$target.term),]
      
      current.bp.mapres$bp.uri.col <- current.bp.mapres$target.term
      current.bp.mapres$bp.uri.col <-
        sub(
          pattern = "http://purl.obolibrary.org/obo/NCBITaxon_",
          replacement = "http://purl.bioontology.org/ontology/NCBITAXON/",
          x = current.bp.mapres$bp.uri.col
        )
      
      # not currently working with some GO terms
      current.bp.mapres.labels <-
        get.bp.labels(
          mappings.frame = current.bp.mapres,
          uri.col.name = "bp.uri.col",
          onto.abbrev.col.name = "target.prefix.uc",
          api.family = "classes"
        )
      
      named.flag <- sapply(current.bp.mapres.labels, length) == 2
      
      current.bp.mapres.labels <-
        current.bp.mapres.labels[named.flag]
      
      current.bp.mapres.labels <-
        do.call(rbind.data.frame, current.bp.mapres.labels)
      
      # get more labels with ols?
      unlabeled.bp.mappees <-
        unique(current.bp.mapres[(!(
          current.bp.mapres$bp.uri.col %in% current.bp.mapres.labels$target.term
        )), c("target.prefix", "target.term")])
      
      
      ols.labels <-
        apply(
          X = unlabeled.bp.mappees,
          MARGIN = 1,
          FUN = function(current_row) {
            current_iri <- current_row[['target.term']]
            current_prefix <-
              tolower(current_row[['target.prefix']])
            ols.base <-
              "https://www.ebi.ac.uk/ols/api/ontologies/"
            ols.incremental <-
              paste0(ols.base, current_prefix, "/terms/")
            first.encoding <-
              URLencode(current_iri, reserved = TRUE, repeated = TRUE)
            second.encoding <-
              URLencode(first.encoding,
                        reserved = TRUE,
                        repeated = TRUE)
            prepared.url <- paste0(ols.incremental, second.encoding)
            ols.result <- httr::GET(prepared.url)
            print(ols.result$status_code)
            if (ols.result$status_code == 200) {
              ols.list <- fromJSON(rawToChar(ols.result$content))
              return(list("target.term" = current_iri, "label" = ols.list$label))
            }
          }
        )
      
      ols.labels <- do.call(rbind.data.frame, ols.labels)
      
      any.label <-
        unique(rbind.data.frame(current.bp.mapres.labels, ols.labels))
      
      current.bp.mapres <-
        left_join(
          x = current.bp.mapres,
          y = any.label,
          by = c("bp.uri.col" = "target.term")
        )
      
      ####
      
      current.bp.mapres <-
        left_join(
          x = current.bp.mapres,
          y = Part,
          by = c("source.bare" = "PartNumber")
        )
      
      if (is.data.frame(paths.from.leaves)) {
        current.bp.mapres <-
          left_join(x = current.bp.mapres,
                    y = paths.from.leaves,
                    by = c("source.bare" = "CODE"))
      } else {
        current.bp.mapres$path = ""
      }
      
      needs.mapping <-
        sort(setdiff(needs.mapping, current.bp.mapres$source.bare))
      
      current.bp.mapres <-
        current.bp.mapres[, c("source.bare",
                              "PartDisplayName",
                              "reassembled",
                              "label",
                              "domain",
                              "path")]
      
      colnames(current.bp.mapres) <- c(
        "LOINC_PartNumber",
        "PartDisplayName",
        "target.term",
        "target.label",
        "target.domain",
        "path.to.LOINC.root"
      )
      
      current.bp.mapres$mapping.method <-
        "Pre-prioritized BioPortal mapping"
      
    }
    
    ####
    
    current.dispnames <-
      unique(max.one.gene.wide[, c(current.code.col.name , current.name.col.name)])
    
    part.vector <-
      unlist(current.dispnames[, 1])
    
    current.dispnames <-
      current.dispnames[part.vector %in% needs.mapping , ]
    
    name.vector <-
      unlist(current.dispnames[, 2])
    current.dispnames <-
      current.dispnames[order(name.vector),]
    
    # refactor
    # need better sorting
    current.ols.search.res.plural <-
      apply(
        X = current.dispnames,
        MARGIN = 1,
        FUN = function(current_row) {
          inner <-  ols.serch.term.labels.universal(
            current.string = current_row[[current.name.col.name]],
            current.id =  current_row[[current.code.col.name]],
            strip.final.s =  FALSE,
            ontology.filter = paste0("&ontology=", onto.list.for.ols),
            kept.row.count =  30,
            req.exact = 'false'
          )
          if (is.data.frame(inner)) {
            return(inner)
          }
        }
      )
    current.ols.search.res.plural <-
      do.call(plyr::rbind.fill, current.ols.search.res.plural)
    
    # need to report total failures
    
    current.dispnames.singular <-
      current.dispnames[grepl(pattern = "s$",
                              x = unlist(current.dispnames[, 2])) ,]
    
    current.dispnames.singular.original <-
      unlist(current.dispnames.singular[, 2])
    current.dispnames.singular.replacement <- sub(pattern = "s$",
                                                  replacement = "",
                                                  x = current.dispnames.singular.original)
    
    current.dispnames.singular[, 2] <-
      current.dispnames.singular.replacement
    
    current.dispnames.singular <-
      current.dispnames.singular[order(current.dispnames.singular[, 2]), ]
    
    current.ols.search.res.singular <-
      apply(
        X = current.dispnames.singular,
        MARGIN = 1,
        FUN = function(current_row) {
          inner <-  ols.serch.term.labels.universal(
            current.string = current_row[[current.name.col.name]],
            current.id =  current_row[[current.code.col.name]],
            strip.final.s =  FALSE,
            ontology.filter = paste0("&ontology=", onto.list.for.ols),
            kept.row.count =  30,
            req.exact = 'false'
          )
          if (is.data.frame(inner)) {
            return(inner)
          }
        }
      )
    
    current.ols.search.res.singular <-
      do.call(plyr::rbind.fill, current.ols.search.res.singular)
    
    
    current.ols.search.res <-
      rbind.data.frame(current.ols.search.res.plural,
                       current.ols.search.res.singular)
    
    current.ols.search.res$synonym <-
      sapply(current.ols.search.res$synonym, function(current_syn) {
        paste(current_syn, collapse = "|")
      })
    
    current.ols.search.res$annotations_trimmed <-
      sapply(current.ols.search.res$annotations_trimmed, function(current_syn) {
        paste(current_syn, collapse = "|")
      })
    
    current.ols.search.res <-
      left_join(x = current.ols.search.res,
                y = harmonized_ontology_rankings,
                by = c("ontology_prefix" = "rhs"))
    
    q.vs.label.cosine <-
      stringdist(a = current.ols.search.res$query,
                 b = current.ols.search.res$label,
                 method = 'cosine')
    
    current.ols.search.res$q.vs.label.cosine <- q.vs.label.cosine
    
    current.ols.search.res.synonyms <-
      apply(
        X = current.ols.search.res,
        MARGIN = 1,
        FUN = function(current_row) {
          synonyms <- vals.by.relation(current_row, 'synonym')
          return(synonyms)
        }
      )
    
    current.ols.search.res.synonyms <-
      do.call(rbind.data.frame, current.ols.search.res.synonyms)
    
    current.ols.search.res.trimmed <-
      apply(
        X = current.ols.search.res,
        MARGIN = 1,
        FUN = function(current_row) {
          synonyms <- vals.by.relation(current_row, 'annotations_trimmed')
          return(synonyms)
        }
      )
    current.ols.search.res.trimmed <-
      do.call(rbind.data.frame, current.ols.search.res.trimmed)
    current.ols.search.res.non.label <-
      rbind.data.frame(current.ols.search.res.synonyms,
                       current.ols.search.res.trimmed)
    
    current.ols.search.res.non.label$other.cosine <-
      stringdist(a = current.ols.search.res.non.label$query,
                 b = current.ols.search.res.non.label$value,
                 method = 'cosine')
    
    current.ols.search.res.non.label.best.cosine <-
      aggregate(
        current.ols.search.res.non.label$other.cosine,
        by = list(
          current.ols.search.res.non.label$loinc.part,
          current.ols.search.res.non.label$obo_id
        ),
        FUN = min
      )
    
    colnames(current.ols.search.res.non.label.best.cosine) <-
      c("loinc.part", "obo_id", "other.cosine")
    
    current.ols.search.res.non.label <-
      inner_join(x = current.ols.search.res.non.label, y = current.ols.search.res.non.label.best.cosine)
    
    #
    
    current.ols.search.res.non.label <-
      current.ols.search.res.non.label[, c("loinc.part",
                                           "obo_id",
                                           "relation",
                                           "value",
                                           "other.cosine")]
    
    current.ols.search.res <-
      inner_join(x = current.ols.search.res, y = current.ols.search.res.non.label)
    
    # tm <- as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S")
    # tsfn <- strftime(tm , "%Y-%m-%d-%H-%M-%S-UTC.Rdata")
    
    # current.ols.search <-
    #   analyte.core.any.map.meth$SYSTEM.ols.search.res
    
    if (is.data.frame(paths.from.leaves)) {
      current.ols.search.res <-
        left_join(x = current.ols.search.res,
                  y = paths.from.leaves,
                  by = c("loinc.part" = "CODE"))
    } else {
      current.ols.search.res$path <- ""
    }
    
    # This isn't triaged yet (but the BP mapping is)
    
    
    current.ols.search.res.metadata <- current.ols.search.res[, c(
      "synonym",
      "annotations_trimmed",
      "rank",
      "q.vs.label.cosine",
      "relation",
      "value",
      "other.cosine",
      "priority"
    )]
    
    
    current.ols.search.res <-
      current.ols.search.res[, c("loinc.part", "query", "iri", "label", "domain", "path")]
    
    
    colnames(current.ols.search.res) <- c(
      "LOINC_PartNumber",
      "PartDisplayName",
      "target.term",
      "target.label",
      "target.domain",
      "path.to.LOINC.root"
    )
    
    current.ols.search.res$mapping.method <- "OLS search"
    
    current.ols.search.res <-
      cbind.data.frame(current.ols.search.res, current.ols.search.res.metadata)
    
    best.cosine <-
      current.ols.search.res[, c("q.vs.label.cosine", "other.cosine")]
    best.cosine <- apply(best.cosine, 1, FUN = min)
    
    current.ols.search.res$best.cosine <- best.cosine
    
    current.ols.search.res <-
      lapply(sort(unique(
        current.ols.search.res$LOINC_PartNumber
      )), function(current_part) {
        print(current_part)
        one.part.frame <-
          current.ols.search.res[current.ols.search.res$LOINC_PartNumber == current_part,]
        max.cosine.per.part <-
          max(one.part.frame$best.cosine, na.rm = TRUE)
        min.cosine.per.part <-
          min(one.part.frame$best.cosine, na.rm = TRUE)
        sclaed.cosines <-
          (max.cosine.per.part - one.part.frame$best.cosine) / (max.cosine.per.part - min.cosine.per.part)
        one.part.frame$best.cosine.scaled <- sclaed.cosines
        return(one.part.frame)
      })
    
    current.ols.search.res <-
      do.call(rbind.data.frame, current.ols.search.res)
    
    # some NA best cosine scaled values
    
    any.method <-
      plyr::rbind.fill(loinc.provided.primary,
                       current.bp.mapres,
                       current.ols.search.res)
    
    any.method$onto.abbrev <-
      sub(
        pattern = "http://purl.obolibrary.org/obo/",
        replacement = "",
        x = any.method$target.term,
        fixed = TRUE
      )
    any.method$onto.abbrev <-
      sub(pattern = "_.*$",
          replacement = "",
          x = any.method$onto.abbrev)
    
    any.method$q.vs.label.cosine <-
      stringdist(
        a = tolower(any.method$PartDisplayName),
        b = tolower(any.method$target.label),
        method = "cosine"
      )
    
    
    ####
    
    return(
      list(
        reviewed.part.frame = reviewed.part.frame ,
        # loinc.provided.primary = loinc.provided.primary,
        loinc.provided.secondary = loinc.provided.secondary,
        # SYSTEM.bp.mapres = SYSTEM.bp.mapres,
        # SYSTEM.ols.search.res = SYSTEM.ols.search.res,
        # paths.from.leaves = paths.from.leaves
        any.method = any.method
      )
    )
  }

####

if (repeat.ehr.loinc.query) {
  ehr.test.and.refresh()
  
  ehr.reference.labs.query <- "
SELECT
	rlri.PK_LAB_RESULT_ITEM_ID ,
	RLRI.LAB_RESULT_ITEM_CODE ,
	RLRI.LAB_RESULT_ITEM_DESCRIPTION ,
	RLRI.LOINC ,
	count(DISTINCT pe.EMPI ) AS UNIQUE_EMPI_COUNT
FROM
	MDM.R_LAB_RESULT_ITEM rlri
JOIN MDM.LAB_RESULTS lr ON
	rlri.PK_LAB_RESULT_ITEM_ID = lr.FK_LAB_RESULT_ITEM_ID
JOIN mdm.PATIENT_ENCOUNTER pe ON
	lr.FK_PATIENT_ENCOUNTER_ID = pe.PK_PATIENT_ENCOUNTER_ID
GROUP BY
	rlri.PK_LAB_RESULT_ITEM_ID ,
	RLRI.LAB_RESULT_ITEM_CODE ,
	RLRI.LAB_RESULT_ITEM_DESCRIPTION ,
	RLRI.LOINC"
  
  print(Sys.time())
  timed.system <- system.time(
    ehr.reference.labs.result <-
      dbGetQuery(conn = ehr.connection, statement = ehr.reference.labs.query)
  )
  print(Sys.time())
  print(timed.system)
  
  # "
  # SELECT
  # 	lr1.PK_LAB_RESULT_ID
  # FROM
  # 	mdm.LAB_RESULTS lr1 MINUS
  # SELECT
  # 	lr2.PK_LAB_RESULT_ID
  # FROM
  # 	mdm.LAB_RESULTS lr2,
  # 	mdm.R_LAB_RESULT_ITEM rlri2
  # WHERE
  # 	lr2.FK_LAB_RESULT_ITEM_ID = rlri2.PK_LAB_RESULT_ITEM_ID"
  # # takes 15 minutes
  # # so far, have not any observed lab results that don't map to a R lab result item
  
  # # Close connection
  # dbDisconnect(ehr.connection)
  save(ehr.reference.labs.result,
       "ehr_reference_labs_result.Rdata")
} else {
  load("ehr_reference_labs_result.Rdata")
}

####

ehr.bad.loincs <-
  ehr.reference.labs.result[(!(is.na(ehr.reference.labs.result$LOINC))) &
                              (!(
                                ehr.reference.labs.result$LOINC %in% LoincPartLink$LoincNumber
                              )) , ]

ehr.na.loincs <-
  ehr.reference.labs.result[is.na(ehr.reference.labs.result$LOINC) , ]

# ####

lpl.min <-
  LoincPartLink[, c("PartTypeName", "LinkTypeName", "Property")]

lpl.wide <-
  lpl.min %>% pivot_wider(names_from = PartTypeName,
                          values_from = LinkTypeName,
                          values_fn = length)

lpl.long <-
  lpl.wide %>% pivot_longer(-Property, names_to = "PartTypeName", values_to = "count")
lpl.long <- lpl.long[complete.cases(lpl.long),]

lpl.long$search.flag <-
  lpl.long$Property == "http://loinc.org/property/search"
lpl.long$doc.flag <-
  grepl(pattern = "http://loinc.org/property/document-",
        x = lpl.long$Property,
        fixed = TRUE)
lpl.long$rad.flag <-
  grepl(pattern = "http://loinc.org/property/rad-",
        x = lpl.long$Property,
        fixed = TRUE)

flag.flag <- grepl(pattern = "flag", x = colnames(lpl.long))
flag.flag <- lpl.long[, flag.flag]
flag.flag[] <- lapply(flag.flag[], as.numeric)
flag.flag$rowSums <- rowSums(flag.flag)
flag.flag <- flag.flag$rowSums == 0
lpl.long <- lpl.long[flag.flag, ]

flag.flag <- grepl(pattern = "flag", x = colnames(lpl.long))
flag.flag <- colnames(lpl.long)[flag.flag]

lpl.long <- lpl.long[, setdiff(colnames(lpl.long), flag.flag)]

lpl.wide <-
  lpl.long %>% pivot_wider(names_from = PartTypeName ,
                           values_from = count,
                           values_fn = length)

lpl.wide$rowSums <-
  rowSums(lpl.wide[, 2:ncol(lpl.wide)], na.rm = TRUE)
print(sort(lpl.wide$rowSums))
# no Property is associated with more than one PartTypeName

lpl.wide.colSums <-
  colSums(lpl.wide[, 2:ncol(lpl.wide)], na.rm = TRUE)
print(sort(lpl.wide.colSums))
# some PartTypeNames are associated with more than one Property

# skip category... almost guaranteed many to many
# some columns may be redundant... check jaccards
# do names and codes???

gene.check <-
  unique(LoincPartLink[LoincPartLink$Property == "http://loinc.org/property/analyte-gene" , c("LoincNumber", "PartName")])
gene.check <- make.table.frame(gene.check$LoincNumber)
gene.check <- gene.check$value[gene.check$count > 1]

# (!(LoincPartLink %in% gene.check)) &

accepted.assays <- setdiff(LoincPartLink$LoincNumber, gene.check)
accepted.properties <-
  setdiff(lpl.wide$Property, excluded.Properties)
lpl <-
  unique(LoincPartLink[LoincPartLink$Property %in% accepted.properties &
                         LoincPartLink$LoincNumber %in% accepted.assays,
                       c("LoincNumber", "PartName", "Property")])

max.one.gene.wide <-
  lpl.wide.builder(names.from = 'Property',
                   values.from = 'PartNumber',
                   to.strip = "http://loinc.org/property/")

codes.jaccard.calcs <- do.jaccard.calcs(max.one.gene.wide, 2)

ehr.reference.labs.result.min.empi.count <-
  ehr.reference.labs.result[(!(is.na(ehr.reference.labs.result$LOINC))), c("LOINC", "UNIQUE_EMPI_COUNT")]
ehr.reference.labs.result.min.empi.count <-
  aggregate(
    ehr.reference.labs.result.min.empi.count$UNIQUE_EMPI_COUNT,
    by = list(ehr.reference.labs.result.min.empi.count$LOINC),
    FUN = max
  )

colnames(ehr.reference.labs.result.min.empi.count) <-
  c("LOINC", "UNIQUE_EMPI_COUNT")

ehr.reference.labs.result.min.empi.count <-
  ehr.reference.labs.result.min.empi.count[ehr.reference.labs.result.min.empi.count$UNIQUE_EMPI_COUNT >= min.empis ,]

Loinc.min <-
  Loinc[Loinc$LOINC_NUM %in% max.one.gene.wide$LoincNumber, c(
    "LOINC_NUM",
    "VersionLastChanged",
    "CHNG_TYPE",
    "DefinitionDescription",
    "STATUS",
    "CONSUMER_NAME",
    "CLASSTYPE",
    "FORMULA",
    "UNITSREQUIRED",
    "SUBMITTED_UNITS",
    "ORDER_OBS",
    "CDISC_COMMON_TESTS",
    "HL7_FIELD_SUBFIELD_ID",
    "EXTERNAL_COPYRIGHT_NOTICE",
    "EXAMPLE_UNITS",
    "LONG_COMMON_NAME",
    "UnitsAndRange",
    "EXAMPLE_UCUM_UNITS",
    "EXAMPLE_SI_UCUM_UNITS",
    "STATUS_REASON",
    "STATUS_TEXT",
    "CHANGE_REASON_PUBLIC",
    "COMMON_TEST_RANK",
    "COMMON_ORDER_RANK",
    "COMMON_SI_TEST_RANK",
    "HL7_ATTACHMENT_STRUCTURE",
    "EXTERNAL_COPYRIGHT_LINK",
    "PanelType",
    "AskAtOrderEntry",
    "AssociatedObservations",
    "VersionFirstReleased",
    "ValidHL7AttachmentRequest"
  )]

####

part.names.sorted <-
  sort(colnames(max.one.gene.wide[2:ncol(max.one.gene.wide)]))
Part.min.cols.merge.prep <- setdiff(Part.min.cols, "PartNumber")
placeholder <-
  lapply(Part.min.cols.merge.prep, function(current_col) {
    print(current_col)
    max.one.gene.wide <<-
      max.one.gene.wide %>% add_column(!!(current_col) := NA)
  })

for (i in 1:length(part.names.sorted)) {
  joinee = "base"
  if (i > 1) {
    joinee <- part.names.sorted[i - 1]
  }
  joiner <- part.names.sorted[i]
  by.trick <- "PartNumber"
  names(by.trick) <- joiner
  print(by.trick)
  max.one.gene.wide <<-
    left_join(
      x = max.one.gene.wide,
      y = Part.min,
      by = by.trick,
      all.x = TRUE,
      suffix = c(
        paste0(".", joinee, collapse = ""),
        paste0(".", joiner, collapse = "")
      )
    )
}

max.one.gene.wide <-
  left_join(x = max.one.gene.wide,
            y = Loinc.min,
            by = c("LoincNumber" = "LOINC_NUM"))

max.one.gene.wide <-
  inner_join(x = max.one.gene.wide,
             y = ehr.reference.labs.result.min.empi.count,
             by = c("LoincNumber" = "LOINC"))

dot.base.cols <-
  grepl(pattern = ".base$", x = colnames(max.one.gene.wide))
not.dot.base.cols <- colnames(max.one.gene.wide)[(!(dot.base.cols))]
max.one.gene.wide <- max.one.gene.wide[, not.dot.base.cols]

common.constrained.assays <-
  get.common.constrained.assays()

####

turbo.url <-
  "https://raw.githubusercontent.com/PennTURBO/Turbo-Ontology/master/ontologies/robot_implementation/turbo_merged.ttl"

ont.tf <- tempfile()
download.file(turbo.url, ont.tf)

turbo.model <-
  rrdf::load.rdf(filename = ont.tf, format = "TURTLE")

turbo.label.query <-
  "select
distinct
?s ?o
where
{
?s <http://www.w3.org/2000/01/rdf-schema#label> ?prel .
# bind(lcase(str(?prel)) as ?o)
bind(str(?prel) as ?o)
filter(isiri(?s))
}"
  
  Sys.time()
  system.time(turbo.labels <-
                rrdf::sparql.rdf(model = turbo.model, sparql = turbo.label.query))
  
  turbo.labels <- as.data.frame(turbo.labels)
  
  check.turbo.labels <- make.table.frame(turbo.labels$s)
  check.turbo.labels <-
    check.turbo.labels$value[check.turbo.labels$count == 1]
  
  turbo.labels <-
    turbo.labels[turbo.labels$s %in% check.turbo.labels ,]
  

####    ####    ####    ####

loinc_to_obo_mapping_reviewed <-
  read_csv("loinc_to_obo_mapping_reviewed.csv")

# loinc_to_obo_mapping_reviewed <- left_join(loinc_to_obo_mapping_reviewed, Part, by = "PartNumber")
# loinc_to_obo_mapping_reviewed <- left_join(loinc_to_obo_mapping_reviewed, turbo.labels, by = c("iris" = "s"))
#
# write.csv(loinc_to_obo_mapping_reviewed, "loinc_to_obo_mapping_reviewed_new.csv", row.names = FALSE)

common.constrained.assays <-
  get.common.constrained.assays(component.code.list = loinc_to_obo_mapping_reviewed$PartNumber)

temp <- Part[Part$PartNumber %in% analytes.from.reviewed.mapping ,]

analytes.from.reviewed.mapping <-
  intersect(loinc_to_obo_mapping_reviewed$PartNumber, Part$PartNumber[Part$PartTypeName == "COMPONENT"])

# temp <- make.table.frame(loinc_to_obo_mapping_reviewed$PartNumber)


####

component.analyte.core.any.map.meth <-
  staged.loinc.mapping(
    needs.mapping = sort(unique(max.one.gene.wide$analyte.core)),
    # reviewed.part.filename = "loinc_to_obo_mapping_reviewed.csv",
    already.reviewed = analytes.from.reviewed.mapping,
    ontology.rankings = ontology.rankings.min,
    current.code.col.name = "analyte.core",
    current.name.col.name = "PartName.analyte.core"
  )

temp <- component.analyte.core.any.map.meth$any.method
temp$plus <-
  grepl(pattern = "+",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$period <-
  grepl(pattern = ".",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$gene <-
  grepl(pattern = "gene",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$cell <-
  grepl(pattern = "cell",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$free <-
  grepl(pattern = "free",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$total <-
  grepl(pattern = "total",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp$bound <-
  grepl(pattern = "bound",
        x = temp$PartDisplayName,
        fixed = TRUE)
temp <-
  unique(temp[, c(
    "LOINC_PartNumber",
    "PartDisplayName",
    "cell",
    "gene",
    "period",
    "plus",
    "free",
    "total",
    "bound",
    "target.term",
    "target.label",
    "onto.abbrev",
    "path.to.LOINC.root",
    "mapping.method",
    "rank",
    "q.vs.label.cosine",
    "other.cosine",
    "best.cosine",
    "best.cosine.scaled",
    "priority",
    "target.domain"
  )])

write_excel_csv(temp, "core_analyte_mapping_needs_review.csv")

####

  
  ####
  
  cd_markers <- read_csv("cd_markers.csv")
  
  temp <-
    apply(
      X = cd_markers,
      MARGIN = 1,
      FUN = function(current_marker) {
        print(current_marker[["PartNumber"]])
        print(current_marker[["PartDisplayName"]])
        
        inner_temp <- ols.serch.term.labels.universal(
          current.string = current_marker[["PartDisplayName"]],
          current.id = current_marker[["PartNumber"]],
          strip.final.s = FALSE,
          ontology.filter = paste0("&ontology=", onto.list.for.ols),
          kept.row.count = 30,
          req.exact = 'false'
        )
        return(inner_temp)
        
      }
    )
  
  temp <- do.call(rbind.data.frame, temp)
  
  
  temp <-
    temp[, c(
      "loinc.part",
      "query",
      "iri",
      "label",
      "rank",
      "short_form",
      "obo_id",
      "ontology_name",
      "ontology_prefix"
    )]
  
  write.csv(temp, file = "cd_markers_needs_review.csv", row.names = FALSE)
  
  ####
  
  loinc_combo_analytes <- read_csv("loinc_combo_analytes.csv")
  loinc_combo_analytes <-
    left_join(loinc_combo_analytes, turbo.labels, by = c("obo1" = "s"))
  loinc_combo_analytes <-
    left_join(loinc_combo_analytes, turbo.labels, by = c("obo2" = "s"))
  
  loinc_combo_analytes$axiom.labels <-
    paste0("'(",
           loinc_combo_analytes$o.x,
           "' and '",
           loinc_combo_analytes$o.y,
           "')")
  
  write.csv(loinc_combo_analytes, "loinc_combo_analytes_new.csv", row.names = FALSE)
  
  # sort(table(PartRelatedCodeMapping$ExtCodeSystem))
  
  ####
  
  # component (core analyte + ?)
  # method (no hierarchy, few good mappings)
  # property
  
  # scale
  # do Qn, (Ord & Nom both "categorical"?... need to create instances)
  # > sort(table(max.one.gene.wide$PartDisplayName.SCALE_TYP))
  # - OrdQn   Doc   Nar   Nom   Ord    Qn 
  # 17    21    25    60   251   721  2313 
  
  # system: GOOD!
  
  # time: make axioms but don't bother searching?
  # do PT, 24 (one EHR-relevant assay each for 72, 5, 48, 12)
  # > sort(table(max.one.gene.wide$PartDisplayName.TIME_ASPCT))
  # 12 hours             48 hours              5 hours             72 hours                    *          Unspecified 
  # 1                    1                    1                    1                    2                   12 
  # 24 hours Point in time (spot) 
  # 120                 3270                 12 
  
  # map
  # constrain
  # add more component subparts
  # no hierarchy found for systems or methods
  # very little useful method mapping
  
  # property
  # scale... using Qn, Ord and Nom only now
  # map property and scale to (output) datum types?
  
  
  # time... using Pt and 24 only right now
  
  # units.brainstorm <-
  #   as.data.frame.matrix(
  #     table(
  #       common.constrained.assays$PartName.PROPERTY,
  #       common.constrained.assays$PartName.SCALE_TYP
  #     )
  #   )
  # units.brainstorm$PROPERTY <- rownames(units.brainstorm)
  # units.brainstorm <-
  #   pivot_longer(data = units.brainstorm,
  #                -PROPERTY,
  #                names_to = "SCALE" ,
  #                values_to = "count")
  # units.brainstorm <- units.brainstorm[units.brainstorm$count > 0 , ]
  # 
  # prop.sum <-
  #   aggregate(units.brainstorm$count,
  #             by = list(units.brainstorm$PROPERTY),
  #             FUN = sum)
  # 
  # names(prop.sum) <- c("PROPERTY", "count.sum")
  # 
  # prop.sum$prop.lc <- tolower(prop.sum$PROPERTY)
  # 
  # absent_needed_props <- read_excel("absent_needed_props.xlsx",
  #                                   col_names = FALSE)
  # colnames(absent_needed_props) <- c("absent_needed_props")
  # absent_needed_props <-
  #   pull(absent_needed_props, absent_needed_props)
  # 
  # prop.sum <- prop.sum[prop.sum$prop.lc %in% absent_needed_props ,]
  # 
  # turbo_loinc_properties_reviewed <- read_csv("turbo_loinc_properties_reviewed.csv")
  # 
  # understood.properties <-
  #   turbo_loinc_properties_reviewed$PartName[(!(is.na(turbo_loinc_properties_reviewed$use))) &
  #                                              turbo_loinc_properties_reviewed$use == "y"]
  
  ####
  
  # # # saves columns, not tied together in any way
  # prop.sum.json <- jsonlite::toJSON(prop.sum)
  # attr(prop.sum.json, "tag") <- "tagattr"
  # yaml::write_yaml(prop.sum.json, "properties_followup.yaml")
  
  # goal here: get the reviewed mappings, the LOINC provided mappings, 
  #   the (highly scrutinized) bp mappings and the minimially scrutinized ols search results 
  #   into one table that's easy to review
  
  # need better column ordering
  # query (loinc) ID, query (loinc) name, OBO:ID, OBO name, quality metrics...
  
  # these don't show ow important the parts are in terms of number of loinc assays,  
  #   number of EHR assays, or (minimum?) number of patients
  
  ####
  
  # column.usefulness <- get.column.usefulness(common.constrained.assays)
  # write.csv(column.usefulness, "column_usefulness.csv", row.names = FALSE)
  
  
  