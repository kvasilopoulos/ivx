
library(readxl)
library(dplyr)


# KMS ---------------------------------------------------------------------


kms <- readxl::read_excel("data-raw/kms-monthly.xlsx",
                      col_types = c("text", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric")) %>%
  mutate(Date = Date %>% lubridate::ymd(truncated = 1)) %>%
  rename(DE = "D/E", DY = "D/Y", DP = "D/P", Ret = "LOG_EXCESS_VW",
         EP = "E/P", BM = "B/M")

kms_quarterly <- readxl::read_excel("data-raw/kms-quarterly.xlsx", col_names = FALSE,
                            col_types = c("text",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric")) %>%
  setNames(c(colnames(kms), "dont")) %>%
  select(-dont) %>%
  mutate(Date = Date  %>%
              zoo::as.yearqtr(format = "%Y%q") %>%
              zoo::as.Date())

usethis::use_data(kms, overwrite = T)
usethis::use_data(kms_quarterly, overwrite = T)


# YLPC --------------------------------------------------------------------

rdiffx <- function(x) c(0, x[-1]/x[-length(x)] - 1)

ylpc <- readr::read_csv("data-raw/ylpc-data.csv") %>%
  mutate_at(vars(hpi), rdiffx) %>%
  dplyr::select(date, cpi, def, gdp, inc, ind, int, inv, mog, res, une, hpi)

# TODO update ordering of the data
usethis::use_data(ylpc, overwrite = T)
