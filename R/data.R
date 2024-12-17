#' MacroPL Dataset
#'
#' This dataset contains macroeconomic variables from Poland, collected from the \url{https://stat.gov.pl} and \url{https://nbp.pl}
#' Each variable has an abbreviation and specific transformation applied for time series analysis.
#'
#' @format Time series object containing 36 monthly macroeconomic indicators from Poland from January 2002 to July 2024.
#' \describe{
#'   \item{Date}{Date of the observation.}
#'   \item{SBR}{State budget revenues. Transform (difference). }
#'   \item{SBE}{State budget expenditures. Transform (difference).}
#'   \item{BSB}{Balance of the state budget. Transform (difference).}
#'   \item{ExC}{Exports of goods (current prices). Transform (difference).}
#'   \item{ExP}{Exports of goods (constant prices). Transform (difference).}
#'   \item{ImC}{Imports of goods (current prices). Transform (difference).}
#'   \item{ImP}{Imports of goods (constant prices). Transform (difference).}
#'   \item{Man}{Manufacturing. Transform (difference).}
#'   \item{Bld}{Building construction. Transform (difference).}
#'   \item{UnE}{Unemployment. Transform (difference).}
#'   \item{UnC}{Unemployment (comparison with the last year). No transformation.}
#'   \item{CPI}{Consumer price indices. Transform (difference).}
#'   \item{PFB}{Prices (food and non-alcoholic beverages). Transform (difference).}
#'   \item{PAT}{Prices (alcoholic beverages and tobacco products). Transform (difference).}
#'   \item{PCF}{Prices (clothing and footwear). Transform (difference).}
#'   \item{PHE}{Prices (housing or home use and energy carriers). Transform (difference).}
#'   \item{PHH}{Prices (home furnishings and housekeeping). Transform (difference).}
#'   \item{PH}{Prices (health). Transform (difference).}
#'   \item{PT}{Prices (transportation). Transform (difference).}
#'   \item{PC}{Prices (communications). Transform (difference).}
#'   \item{PRC}{Prices (recreation and culture). Transform (difference).}
#'   \item{PE}{Prices (education). Transform (difference).}
#'   \item{DPU}{Dwellings put into use. No transformation.}
#'   \item{TEx}{Transaction price indices of exports. No transformation.}
#'   \item{TIm}{Transaction price indices of imports. No transformation.}
#'   \item{ToT}{Terms of trade. No transformation.}
#'   \item{M0}{M0 money aggregate. Transform (difference).}
#'   \item{M1}{M1 money aggregate. Transform (difference).}
#'   \item{M2}{M2 money aggregate. Transform (difference).}
#'   \item{M3}{M3 money aggregate. Transform (difference).}
#'   \item{ERE}{Exchange rate Euro (end of the month). No transformation.}
#'   \item{ERA}{Exchange rate Euro (monthly average). No transformation.}
#'   \item{USE}{Exchange rate USD (end of the month). No transformation.}
#'   \item{USA}{Exchange rate USD (monthly average). No transformation.}
#'   \item{Ref}{Reference rate. No transformation.}
#'   \item{Lom}{Lombard rate. No transformation.}
#' }
#'
#' @source \url{https://stat.gov.pl}, \url{https://nbp.pl}
#' @note First 22 variables can be treated as slow-moving.
#' @name MacroPL
#' @docType data
#' @keywords datasets
"MacroPL"