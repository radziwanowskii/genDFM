#' MacroPL Dataset
#'
#' This dataset contains macroeconomic variables from Poland, collected from various official sources.
#' Each variable has an abbreviation and specific transformation applied for time series analysis.
#'
#' @format A data frame with several columns:
#' \describe{
#'   \item{Date}{Date of the observation.}
#'   \item{SBR}{State budget revenues. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{SBE}{State budget expenditures. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{BSB}{Balance of the state budget. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{ExC}{Exports of goods (current prices). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{ExP}{Exports of goods (constant prices). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{ImC}{Imports of goods (current prices). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{ImP}{Imports of goods (constant prices). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{Man}{Manufacturing. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{Bld}{Building construction. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{UnE}{Unemployment. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{UnC}{Unemployment (comparison with the last year). No transformation. Source: \url{https://stat.gov.pl}}
#'   \item{CPI}{Consumer price indices. Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PFB}{Prices (food and non-alcoholic beverages). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PAT}{Prices (alcoholic beverages and tobacco products). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PCF}{Prices (clothing and footwear). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PHE}{Prices (housing or home use and energy carriers). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PHH}{Prices (home furnishings and housekeeping). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PH}{Prices (health). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PT}{Prices (transportation). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PC}{Prices (communications). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PRC}{Prices (recreation and culture). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{PE}{Prices (education). Transform (difference). Source: \url{https://stat.gov.pl}}
#'   \item{DPU}{Dwellings put into use. No transformation. Source: \url{https://stat.gov.pl}}
#'   \item{TEx}{Transaction price indices of exports. No transformation. Source: \url{https://stat.gov.pl}}
#'   \item{TIm}{Transaction price indices of imports. No transformation. Source: \url{https://stat.gov.pl}}
#'   \item{ToT}{Terms of trade. No transformation. Source: \url{https://stat.gov.pl}}
#'   \item{M0}{M0 money aggregate. Transform (difference). Source: \url{https://nbp.pl}}
#'   \item{M1}{M1 money aggregate. Transform (difference). Source: \url{https://nbp.pl}}
#'   \item{M2}{M2 money aggregate. Transform (difference). Source: \url{https://nbp.pl}}
#'   \item{M3}{M3 money aggregate. Transform (difference). Source: \url{https://nbp.pl}}
#'   \item{ERE}{Exchange rate Euro (end of the month). No transformation. Source: \url{https://nbp.pl}}
#'   \item{ERA}{Exchange rate Euro (monthly average). No transformation. Source: \url{https://nbp.pl}}
#'   \item{USE}{Exchange rate USD (end of the month). No transformation. Source: \url{https://nbp.pl}}
#'   \item{USA}{Exchange rate USD (monthly average). No transformation. Source: \url{https://nbp.pl}}
#'   \item{Ref}{Reference rate. No transformation. Source: \url{https://nbp.pl}}
#'   \item{Lom}{Lombard rate. No transformation. Source: \url{https://nbp.pl}}
#' }
#'
#' @source \url{https://stat.gov.pl}, \url{https://nbp.pl}
#' @note First 22 variables can be treated as slow-moving.
#' @name MacroPL
#' @docType data
#' @keywords datasets
"MacroPL"