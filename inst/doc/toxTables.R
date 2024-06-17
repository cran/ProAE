## ----include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

knitr::opts_chunk$set(eval = TRUE, message = FALSE, results = 'asis', comment='')
options(width = 200)


## ----load-data----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ProAE)
require(knitr)

require(kableExtra)

data(tox_acute)


## ----results='markup'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

acute <- tox_acute

str(acute)
  

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_1 <- toxTables(dsn = acute,
                    id_var="id",
                    cycle_var="Cycle",
                    baseline_val = 1,
                    arm="arm")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

knitr::kable(table_1$individual)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

knitr::kable(table_1$individual,
             col.names = c('Item/Attribute', 'Drug', 'Placebo', 'Drug', 'Placebo', "p", "Drug", "Placebo", "p"),
             caption = "Table 1. Frequency distributions of patients with nausea score >= 1 and >= 3, by treatment arm") %>%               
             kableExtra::add_header_above(c(" ", "n" = 2,  "Score >= 1" = 2, " ", "Score >= 3" = 2, " ")) %>% 
             kableExtra::add_footnote("p: result from a standard chi square test comparing score threshold frequency distributions between arms",
                          notation = "symbol")


