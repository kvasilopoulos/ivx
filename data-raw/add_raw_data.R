
library(readxl)
library(dplyr)

datam <- readxl::read_excel("data-raw/monthly.xlsx",
                      col_types = c("text", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric")) %>%
  mutate(Date = Date %>% lubridate::ymd(truncated = 1)) %>%
  rename(D_E = "D/E", D_Y = "D/Y", D_P = "D/P", Ret = "LOG_EXCESS_VW",
         E_P = "E/P", B_M = "B/M")

dataq <- readxl::read_excel("data-raw/quarterly.xlsx", col_names = FALSE,
                            col_types = c("text",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric")) %>%
  setNames(c(colnames(datam), "dont")) %>%
  mutate(Date = Date  %>%
              zoo::as.yearqtr(format = "%Y%q") %>%
              zoo::as.Date())

# usethis::use_data(datam, overwrite = T)
# usethis::use_data(dataq, overwrite = T)
