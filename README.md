# genDFM

genDFM is a general package to estimate dynamic factor models in R. It provides a set of functions to estimate the model, plot the components, and compute forecasts. The package includes a wide range of information criteria for the choice of factors and allows for the estimation of Factor Augmented Vector Autoregressions and computation of IRFs. The core of the package is written in C++ and uses the Armadillo library.

## Installation

You can install the development version of genDFM from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("radziwanowskii/genDFM")
```

## Load Package 
When the package is intalled, you can load the package
```r
library(genDFM)
```

## Load the data
Load the MacroPL dataset and remove the first column containing the dates
```r
data(MacroPL)
df <- MacroPL[,-1]
```

