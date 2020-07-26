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
  # current.leaves <- analytes.from.reviewed.mapping
  
  current.leaves <-
    MultiAxialHierarchy[MultiAxialHierarchy$CODE %in% current.leaves, ]
  
  current.leaves.ancestors <-
    lapply(sort(unique(current.leaves$CODE)), function(current.code) {
      # current.code <- "LP14304-7"
      print(paste0("Gettering heirarchy for ", current.code))
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
    
    multi.seq.check <-
      unique(current.leaves.ancestors[, c("CODE", "SEQUENCE")])
    multi.seq.check <-
      multi.seq.check[(!(is.na(multi.seq.check$SEQUENCE))), ]
    multi.seq.check <- make.table.frame(multi.seq.check$CODE)
    
    single.seq.leaves <-
      multi.seq.check$value[multi.seq.check$count == 1]
    multi.seq.check <-
      multi.seq.check$value[multi.seq.check$count > 1]
    multi.seq.check <-
      current.leaves.ancestors[current.leaves.ancestors$CODE %in% multi.seq.check , ]
    single.seq.leaves <-
      current.leaves.ancestors[current.leaves.ancestors$CODE %in% single.seq.leaves , ]
    
    multi.seq.check$SEQUENCE <- as.numeric(multi.seq.check$SEQUENCE)
    
    # could also get longest path or shortest
    first.among.equals <-
      aggregate(
        multi.seq.check$SEQUENCE,
        by = list(multi.seq.check$CODE),
        FUN = min
      )
    names(first.among.equals) <- c("CODE", "SEQUENCE")
    
    multi.seq.check <-
      base::merge(x = multi.seq.check, y = first.among.equals)
    
    single.seq.leaves <-
      rbind.data.frame(single.seq.leaves, multi.seq.check)
    
    
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

staged.loinc.mapping <-
  function(needs.mapping,
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
    
    return(
      list(
        loinc.provided.secondary = loinc.provided.secondary,
        paths.from.leaves = paths.from.leaves,
        any.method = any.method
        # reviewed.part.frame = reviewed.part.frame ,
        # loinc.provided.primary = loinc.provided.primary,
        # SYSTEM.bp.mapres = SYSTEM.bp.mapres,
        # SYSTEM.ols.search.res = SYSTEM.ols.search.res,
      )
    )
  }