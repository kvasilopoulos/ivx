library(readxl)
monthly <- read_excel("data-raw/monthly.xlsx",
                      col_types = c("text", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric")) %>% as.data.frame() %>%
  mutate(Date %>% ymd(truncated = 1))


# usethis::use_data(monthly, overwrite = T)
