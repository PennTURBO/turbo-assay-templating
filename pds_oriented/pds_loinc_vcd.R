library(assertr)
library(vcd)
library(readr)

input <-
  read_csv(
    "C:/Users/Mark Miller/Desktop/PDS_labs_promising_LOINC_annotated_only_and_tabulation_report.csv"
  )

# the analyses below don't weight anything by thePDS order counts yet
# assays may appear on multiple rows
#  for example if they can take blood ,serum, or plasma as their system
#  and the component could be "glucose" and "post"!

common.loinc.parts <-
  unique(input[,
               c("LOINC",
                 "ordercount",
                 "prop",
                 "time",
                 "sys",
                 "scale")])

# common.loinc.parts <-
#   common.loinc.parts[# common.loinc.parts$time %in% c("Pt", "24H", "XXX") &
#     # common.loinc.parts$scale %in% c("Pt", "24H", "XXX") &
#     common.loinc.parts$prop %in% c("MCnc",
#                                    "PrThr",
#                                    "ACnc",
#                                    "Prid",
#                                    "SCnc",
#                                    "NFr",
#                                    "Titr",
#                                    "Imp",
#                                    "ArVRat") ,]

common.loinc.parts <-
  common.loinc.parts[common.loinc.parts$prop %in% c("MCnc",
                                                    "PrThr",
                                                    "ACnc",
                                                    "Prid",
                                                    "SCnc",
                                                    "NFr",
                                                    "Titr",
                                                    "Imp",
                                                    "ArVRat") , ]


# viscols <-
#   setdiff(names(common.loinc.parts), c("LOINC", "ordercount"))
#
# viscols <- c("prop", "time", "sys", "scale")

viscols <- c("time", "prop", "scale")

common.loinc.parts <- common.loinc.parts[, viscols]

dim(common.loinc.parts)
head(common.loinc.parts)

tabulated <- xtabs(data = common.loinc.parts)

length(tabulated)
object.size(tabulated)

pdf("loinc_mosaic.pdf")
mosaic(tabulated, zero_size = 0)
dev.off()


# concatted <- col_concat(temp[, tocat], sep = ";")
#
# result <- sort(table(concatted), decreasing = TRUE)
# result.names <- names(result)
# result <- as.numeric(result)
# names(result) <- result.names
#
# result <- result[result > 100]
#
# result <- result[result$`as.numeric(result)` > 100 ,]
#
# barplot(result)

temp <-
  input[input$time == "Pt" &
          input$scale == "Qn" &
          grepl(pattern = "^[A-Z]Cnc$", x = input$prop) &
          is.na(input$meth) &
          input$sys %in% c("Ser", "Plas", "Ser/Plas", "Bld", "Urine"),]

head(sort(table(temp$comp, useNA = "always"), decreasing = TRUE))
