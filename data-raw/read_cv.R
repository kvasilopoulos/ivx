library(pdftools)
library(stringr)
library(dplyr)

options(stringsAsFactors = FALSE)
text <- pdf_text("./data-raw/crit.pdf")
text2 <- strsplit(text, "\n")
dfex <- list()
delta <- matrix()
j <- 1

for (i in 9:28) {
  temp <- text2[[i]] # when you will sort out the df-gls remove 1col
  if (i %% 2 == 1) {
    temp2 <- temp[-c(1:3, 34:35)]
    delta[j] <- text2[[i]][2]
  }else{
    temp2 <- temp[-c(1:2, 34)]
    delta[j] <- text2[[i]][1]
  }
  extr <- str_extract_all(temp2, "[+-0-9.0-9]+")
  dfex[[j]] <- do.call(rbind.data.frame, extr) %>%
    sapply(as.numeric)
  colnames(dfex[[j]]) <- NULL
  j <- j + 1
}
c1 <- do.call(cbind, dfex[seq(1,19,2)])
c2 <- do.call(cbind, dfex[seq(2,20,2)])
c12 <- rbind(c1,c2)#[, -rm_gls]
rm_gls <- seq(1, NCOL(c12),9)
c_crit <- c12[, -rm_gls]


c_odd <- seq(1, NCOL(c_crit) - 1, 2)
c_even <- seq(2, NCOL(c_crit), 2)

delta_sim <- c(-0.999, seq(-0.975, -0.025, 0.025))
dfgls_sim <- seq(1, -5, -0.1)

c_low_crit <- matrix(c_crit[, c_odd], length(dfgls_sim), length(delta_sim),
                     dimnames = list(dfgls_sim, delta_sim))
c_up_crit <- matrix(c_crit[, c_even], length(dfgls_sim), length(delta_sim),
                    dimnames = list(dfgls_sim, delta_sim))
options(stringsAsFactors = TRUE)


# usethis::use_data(c_low_crit, overwrite = T)
# usethis::use_data(c_up_crit, overwrite = T)
