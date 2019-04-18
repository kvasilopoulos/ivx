
library(readxl)
library(dplyr)

monthly <- readxl::read_excel("data-raw/monthly.xlsx",
                      col_types = c("text", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric")) %>%
  mutate(Date = Date %>% lubridate::ymd(truncated = 1)) %>%
  rename(DE = "D/E", DY = "D/Y", DP = "D/P", Ret = "LOG_EXCESS_VW",
         EP = "E/P", BM = "B/M")

quarterly <- readxl::read_excel("data-raw/quarterly.xlsx", col_names = FALSE,
                            col_types = c("text",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric")) %>%
  setNames(c(colnames(datam), "dont")) %>%
  select(-dont) %>%
  mutate(Date = Date  %>%
              zoo::as.yearqtr(format = "%Y%q") %>%
              zoo::as.Date())

usethis::use_data(monthly, overwrite = T)
usethis::use_data(quarterly, overwrite = T)
