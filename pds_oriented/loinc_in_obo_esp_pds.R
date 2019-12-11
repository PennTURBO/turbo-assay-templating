library(SPARQL)
library(config)
library(httr)
library(jsonlite)
library(rJava)
library(RJDBC)

config.yaml.file <- "~/labs_config.yaml"

jdbcDriver <-
  JDBC("oracle.jdbc.OracleDriver", classPath = "C:/ojdbc8.jar")

config.bootstrap <-
  config::get(file = config.yaml.file)

selected.relational.ehr.configuration <-
  config::get(config = "pds",
              file = config.yaml.file)

selected.gdb.configuration <-
  config::get(config = config.bootstrap$selected.gdb.configuration,
              file = config.yaml.file)

# does this ever get used?
bp.api.key <- config.bootstrap$bp.api.key

graphdb.address.port <-
  selected.gdb.configuration$graphdb.address.port
selected.repo <- selected.gdb.configuration$selected.repo
api.user <- selected.gdb.configuration$api.user
api.pass <- selected.gdb.configuration$api.pass

debug.flag <- config.bootstrap$debug.flag

host.name <- selected.relational.ehr.configuration$host
database.name <- selected.relational.ehr.configuration$dbname
db.user <- selected.relational.ehr.configuration$dbuser
# postgres.passwd <- rstudioapi::askForPassword("Database password")
db.passwd <- selected.relational.ehr.configuration$dbpass
db.port <- selected.relational.ehr.configuration$port

jdbcConnection  <- dbConnect(
  jdbcDriver,
  paste0(
    "jdbc:oracle:thin:@",
    host.name,
    ":",
    db.port,
    "/",
    database.name
  ),
  db.user,
  db.passwd
)

table2 <- dbGetQuery(jdbcConnection, 'SELECT
	*
FROM
	mdm.R_LOCATION rl')
