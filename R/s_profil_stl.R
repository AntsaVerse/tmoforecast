#' Extract the Monthly Seasonal Profile from an STL Decomposition
#'
#' Performs an STL (Seasonal-Trend decomposition using Loess) decomposition
#' of a monthly time series and extracts its seasonal component. The
#' seasonal component is then aggregated by calendar month to obtain an
#' average monthly seasonal profile.
#'
#' @param y A univariate time-series object of class \code{ts}, preferably
#'   with monthly frequency (\code{frequency = 12}).
#'
#' @return A list containing:
#' \describe{
#'   \item{stl}{The complete STL decomposition returned by
#'   [stats::stl()].}
#'   \item{profil}{A data frame containing the average seasonal component
#'   for each calendar month. It contains two columns:
#'   \code{mois}, the month number from 1 to 12, and \code{seasonal},
#'   the mean seasonal component for that month.}
#' }
#'
#' @details
#' The function uses a periodic seasonal component by setting
#' \code{s.window = "periodic"} in [stats::stl()]. This assumes that the
#' seasonal pattern is stable over time.
#'
#' The seasonal component extracted from the STL decomposition is matched
#' with the corresponding calendar month. The monthly seasonal profile is
#' then calculated as the arithmetic mean of the seasonal component over
#' all available years.
#'
#' For an additive decomposition, the seasonal component represents the
#' average seasonal contribution to the observed series. The resulting
#' profile can therefore be used to identify months with positive or
#' negative seasonal effects and to construct seasonal adjustments or
#' seasonal forecasts.
#'
#' @examples
#' \dontrun{
#' # Monthly time series
#' y <- ts(
#'   AirPassengers,
#'   start = c(1949, 1),
#'   frequency = 12
#' )
#'
#' # Extract STL decomposition and seasonal profile
#' result <- s_profil_stl(y)
#'
#' # STL decomposition
#' result$stl
#'
#' # Monthly seasonal profile
#' result$profil
#' }
#'
#' @importFrom stats stl time
#' @importFrom zoo as.yearmon
#' @importFrom dplyr mutate group_by summarise
#' @importFrom lubridate month
#'
#' @export
s_profil_stl <- function(y) {

  #==========================================================
  # 1. Décomposition STL
  #==========================================================

  stl_fit <- stats::stl(
    y,
    s.window = "periodic"
  )

  #==========================================================
  # 2. Extraction de la composante saisonnière
  #==========================================================

  seasonal <- stl_fit$time.series[, "seasonal"]

  #==========================================================
  # 3. Création de la base historique
  #==========================================================

  season_df <- data.frame(
    date = as.Date(
    format(
      zoo::as.yearmon(stats::time(y)),
      "%Y-%m-01"
    )
  ),
    seasonal = as.numeric(seasonal)
  ) |>
    dplyr::mutate(
      mois = lubridate::month(date)
    )

  #==========================================================
  # 4. Profil saisonnier mensuel
  #==========================================================

  profil_saisonnier <- season_df |>
    dplyr::group_by(mois) |>
    dplyr::summarise(
      seasonal = mean(
        seasonal,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  return(
    list(

      stl = stl_fit,

      profil = profil_saisonnier

    )
  )
}