# ivx 1.1.0

* Added `ivx_ar` that implements Yang, B., Long, W., Peng, L., & Cai, Z. (2020) 
new instrumental variable based  Wald statistic which accounts for serial 
correlation and heteroscedasticity in the error terms of the linear predictive regression model.
* Added the Yang et al. (2020) dataset named `ylpc`.
* Renamed the `monthly` and `quarterly` dataset into `kms` and `kms_quarterly`
* Remove dependecy on `tibble` and `magrittr`.
* Added `texreg` functionality that converts regression output to LaTeX or HTML tables.
Specifically added `extract` methods for `ivx` and `ivx_ar`, which coefficients and GOF measures 
from a statistical object. 

# ivx 1.0.0

* Added a `NEWS.md` file to track changes to the package.
