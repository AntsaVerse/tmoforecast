#' Forecast a Monthly Time Series Using Rolling-Window Model Selection
#'
#' Automatically compares several univariate time-series forecasting
#' models using rolling-window validation, selects the best-performing
#' model according to a composite performance score, and produces
#' forecasts for a user-defined horizon.
#'
#' The function first cleans the input series using [forecast::tsclean()].
#' It then evaluates several forecasting models using a rolling-window
#' validation procedure. The models compared are:
#' \itemize{
#'   \item STL-ARIMA
#'   \item ARIMA
#'   \item STL-ETS
#'   \item ETS
#'   \item BATS
#'   \item TBATS
#'   \item NNETAR
#'   \item Theta
#'   \item Naive
#'   \item Seasonal Naive
#' }
#'
#' For each model, RMSE, MAE and MAPE are calculated over the validation
#' windows. A composite score is then computed as the sum of the
#' standardized RMSE, MAE and MAPE. The model with the lowest score is
#' selected and re-estimated using the complete time series.
#'
#' Finally, the selected model is used to generate forecasts and several
#' residual diagnostics are computed, including the Ljung-Box test,
#' Jarque-Bera test and the autocorrelation function.
#'
#' @param serie_data A univariate time-series object of class \code{ts}.
#'   The series should have a monthly frequency (frequency = 12).
#' @param h_validation Integer. Number of periods used as the forecasting
#'   horizon during rolling-window validation. Default is 12.
#' @param window_size Integer. Number of observations included in each
#'   training window. Default is 60.
#' @param h_forecast Integer. Number of periods to forecast after selecting
#'   the best model. Default is 24.
#'
#' @return A list containing:
#' \describe{
#'   \item{serie}{The cleaned input time series after applying
#'   [forecast::tsclean()].}
#'   \item{performance}{A data frame containing the performance indicators
#'   RMSE, MAE, MAPE and the composite Score for each candidate model,
#'   ordered from best to worst model.}
#'   \item{best_model}{Character string giving the name of the selected
#'   forecasting model.}
#'   \item{fit}{The final fitted model corresponding to the selected model.}
#'   \item{forecast}{A forecast object containing the forecasts and
#'   prediction intervals.}
#'   \item{diagnostics}{A list containing the Ljung-Box test,
#'   Jarque-Bera test and autocorrelation function of the residuals.}
#'   \item{plot}{A \code{ggplot} object displaying the historical series
#'   and the forecasts.}
#'   \item{database}{A data frame containing the historical and forecast
#'   values, with the columns \code{Date}, \code{Type} and \code{Valeur}.}
#' }
#'
#' @details
#' The rolling-window validation starts after \code{window_size}
#' observations and continues until the end of the series while preserving
#' a validation horizon of \code{h_validation} observations.
#'
#' The composite performance score is calculated as:
#'
#' \deqn{
#' Score = Z(RMSE) + Z(MAE) + Z(MAPE)
#' }
#'
#' where \eqn{Z(.)} denotes the standardized value of each performance
#' indicator. Lower values of the score indicate better forecasting
#' performance.
#'
#' The final model is re-estimated using the complete cleaned series before
#' producing the final forecasts.
#'
#' @section Models:
#' The following forecasting models are evaluated:
#' \describe{
#'   \item{STL_ARIMA}{STL decomposition followed by ARIMA modelling.}
#'   \item{ARIMA}{Automatic seasonal ARIMA model selection.}
#'   \item{STL_ETS}{STL decomposition followed by ETS modelling.}
#'   \item{ETS}{Exponential smoothing state-space model.}
#'   \item{BATS}{BATS model for complex seasonal patterns.}
#'   \item{TBATS}{TBATS model for complex and multiple seasonal patterns.}
#'   \item{NNETAR}{Neural-network autoregressive model.}
#'   \item{THETA}{Theta forecasting method.}
#'   \item{NAIVE}{Naive forecasting method.}
#'   \item{SNAIVE}{Seasonal naive forecasting method.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example monthly time series
#' data <- ts(
#'   AirPassengers,
#'   start = c(1949, 1),
#'   frequency = 12
#' )
#'
#' # Run the forecasting procedure
#' result <- forecast_series(
#'   serie_data = data,
#'   h_validation = 12,
#'   window_size = 60,
#'   h_forecast = 24
#' )
#'
#' # Selected model
#' result$best_model
#'
#' # Model performance
#' result$performance
#'
#' # Forecast object
#' result$forecast
#'
#' # Historical and forecast database
#' head(result$database)
#'
#' # Forecast plot
#' result$plot
#'
#' # Ljung-Box test
#' result$diagnostics$Ljung_Box
#' }
#'
#' @importFrom forecast tsclean stlm auto.arima ets bats tbats nnetar
#'   thetaf naive snaive forecast residuals autoplot
#' @importFrom Metrics rmse mae
#' @importFrom stats Box.test acf as.numeric na.omit
#'   time window residuals
#' @importFrom zoo as.yearmon
#' @importFrom ggplot2 ggtitle
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
forecast_series <- function(
    serie_data,
    h_validation = 12,
    window_size = 60,
    h_forecast = 24
) {

  t0 <- Sys.time()

  y <- forecast::tsclean(serie_data)

  #==============================================================
  # Paramètres de la validation
  #==============================================================

  h <- h_validation

  # Liste des modèles comparés

  models <- c(
    "STL_ARIMA",
    "ARIMA",
    "STL_ETS",
    "ETS",
    "BATS",
    "TBATS",
    "NNETAR",
    "THETA",
    "NAIVE",
    "SNAIVE"
  )

  #==============================================================
  # Barre de progression
  #==============================================================

  n_windows <- length(window_size:(length(y) - h))
  total_steps <- n_windows * length(models)

  pb <- utils::txtProgressBar(
    min = 0,
    max = total_steps,
    style = 3
  )

  counter <- 0

  cat("\n")
  cat(
    "FIKAROHANA NY MODELY MAMARITRA NY FIVOARAHAN'ILAY ANTOTAN'ISA\n"
  )
  cat("Isan'ny modely kajiana :", length(models), "\n")
  cat("Isan'ny windows hanaovana test :", n_windows, "\n\n")

  #==============================================================
  # Fonction permettant d'estimer un modèle et produire
  # une prévision
  #==============================================================

  fit_forecast <- function(train, model, h) {

    switch(

      model,

      "STL_ARIMA" = forecast::forecast(
        forecast::stlm(
          train,
          s.window = "periodic",
          modelfunction = forecast::auto.arima
        ),
        h = h
      ),

      "ARIMA" = forecast::forecast(
        forecast::auto.arima(
          train,
          seasonal = TRUE,
          stepwise = FALSE,
          approximation = FALSE
        ),
        h = h
      ),

      "STL_ETS" = forecast::forecast(
        forecast::stlm(
          train,
          s.window = "periodic",
          modelfunction = forecast::ets
        ),
        h = h
      ),

      "ETS" = forecast::forecast(
        forecast::ets(train),
        h = h
      ),

      "BATS" = forecast::forecast(
        forecast::bats(train),
        h = h
      ),

      "TBATS" = forecast::forecast(
        forecast::tbats(train),
        h = h
      ),

      "NNETAR" = forecast::forecast(
        forecast::nnetar(train),
        h = h
      ),

      "THETA" = forecast::thetaf(
        train,
        h = h
      ),

      "NAIVE" = forecast::naive(
        train,
        h = h
      ),

      "SNAIVE" = forecast::snaive(
        train,
        h = h
      )
    )
  }

  #==============================================================
  # Validation Rolling Window
  #==============================================================

  errors <- vector("list", length(models))
  actuals <- vector("list", length(models))
  residuals_models <- vector("list", length(models))

  names(errors) <- models
  names(actuals) <- models
  names(residuals_models) <- models

  cat("Miketrika ny modely za zao namana.\n")

  for (i in window_size:(length(y) - h)) {

    train <- window(
      y,
      start = time(y)[i - window_size + 1],
      end = time(y)[i]
    )

    test <- y[(i + 1):(i + h)]

    for (m in models) {

      counter <- counter + 1

      utils::setTxtProgressBar(
        pb,
        counter
      )

      fc <- fit_forecast(
        train,
        m,
        h
      )

      errors[[m]] <- c(
        errors[[m]],
        as.numeric(test - fc$mean)
      )

      actuals[[m]] <- c(
        actuals[[m]],
        as.numeric(test)
      )

      if (!is.null(stats::residuals(fc))) {

        residuals_models[[m]] <- c(
          residuals_models[[m]],
          stats::residuals(fc)
        )
      }
    }
  }

  close(pb)

  cat("\n")
  cat("\n")
  cat("Fandrefesana ireo indicateurs de performance...\n\n")

  #==============================================================
  # Evaluation des performances
  #==============================================================

  perf <- data.frame(

    Modele = models,

    RMSE = sapply(
      errors,
      function(x)
        sqrt(mean(x^2, na.rm = TRUE))
    ),

    MAE = sapply(
      errors,
      function(x)
        mean(abs(x), na.rm = TRUE)
    ),

    MAPE = sapply(
      seq_along(models),
      function(i) {

        mean(
          abs(errors[[i]] / actuals[[i]]),
          na.rm = TRUE
        ) * 100
      }
    )
  )

  perf$Score <-
    scale(perf$RMSE)[, 1] +
    scale(perf$MAE)[, 1] +
    scale(perf$MAPE)[, 1]

  results <- perf[
    order(perf$Score),
  ]

  print(results)

  #==============================================================
  # Sélection du meilleur modèle
  #==============================================================

  best_model <- results$Modele[1]

  #==============================================================
  # Réestimation finale
  #==============================================================

  cat(
    "Sempotra ihany fa ny modely tsara indrindra zany dia :",
    best_model,
    "...\n"
  )

  cat(
    "Hamerina hikajy ny :",
    best_model,
    " amin'ny antontan'isa manontolo...\n"
  )

  best_fit <- switch(

    best_model,

    "ARIMA" =
      forecast::auto.arima(
        y,
        seasonal = TRUE,
        stepwise = FALSE,
        approximation = FALSE
      ),

    "ETS" =
      forecast::ets(y),

    "TBATS" =
      forecast::tbats(y),

    "BATS" =
      forecast::bats(y),

    "NNETAR" =
      forecast::nnetar(y),

    "STL_ARIMA" =
      forecast::stlm(
        y,
        s.window = "periodic",
        modelfunction = forecast::auto.arima
      ),

    "STL_ETS" =
      forecast::stlm(
        y,
        s.window = "periodic",
        modelfunction = forecast::ets
      ),

    NULL
  )

  #==============================================================
  # Prévision finale
  #==============================================================

  forecast12 <- switch(

    best_model,

    "THETA" =
      forecast::thetaf(
        y,
        h = h_forecast,
        level = 95
      ),

    "NAIVE" =
      forecast::naive(
        y,
        h = h_forecast,
        level = 95
      ),

    "SNAIVE" =
      forecast::snaive(
        y,
        h = h_forecast,
        level = 95
      ),

    forecast::forecast(
      best_fit,
      h = h_forecast,
      level = 95
    )
  )

  #==============================================================
  # Diagnostics des résidus
  #==============================================================

  if (is.null(best_fit)) {

    ljung_box <- NULL
    jarque_bera <- NULL
    acf_res <- NULL

  } else {

    res <- stats::residuals(best_fit)

    res <- na.omit(res)

    ljung_box <- stats::Box.test(
      res,
      lag = 12,
      type = "Ljung"
    )

    jarque_bera <- tseries::jarque.bera.test(
      as.numeric(res)
    )

    acf_res <- stats::acf(
      res,
      plot = FALSE
    )
  }

  #==============================================================
  # Graphique
  #==============================================================

  cat("Faminavinana ny ho avy...\n")

  forecast_plot <- forecast::autoplot(
    forecast12
  ) +
    ggplot2::ggtitle(
      paste(
        "Prévision - Modèle retenu :",
        best_model
      )
    )

  #==============================================================
  # Export historique + prévision
  #==============================================================

  historical_dates <- time(y)

  historical_data <- data.frame(

    Date = as.Date(
      zoo::as.yearmon(historical_dates)
    ),

    Type = "Historique",

    Valeur = as.numeric(y)
  )

  forecast_dates <- time(
    forecast12$mean
  )

  forecast_data <- data.frame(

    Date = as.Date(
      zoo::as.yearmon(forecast_dates)
    ),

    Type = "Prévision",

    Valeur = as.numeric(
      forecast12$mean
    )
  )

  forecast_database <- rbind(
    historical_data,
    forecast_data
  )

  forecast_database <- forecast_database[
    order(forecast_database$Date),
  ]

  #==============================================================
  # Sortie de la fonction
  #==============================================================

  cat("\n")
  cat("=============================================\n")
  cat("Vita ara-dalana ny vinavina.\n")
  cat("Ny modely voasivana dia :", best_model, "\n")
  cat("=============================================\n\n")

  print(forecast_plot)

  elapsed <- Sys.time() - t0

  cat(
    "\n"
  )

  cat(
    "Fotoana lany :",
    round(
      as.numeric(elapsed, units = "mins"),
      1
    ),
    "minitra\n"
  )

  return(

    list(

      serie = y,

      performance = results,

      best_model = best_model,

      fit = best_fit,

      forecast = forecast12,

      diagnostics = list(
        Ljung_Box = ljung_box,
        Jarque_Bera = jarque_bera,
        ACF = acf_res
      ),

      plot = forecast_plot,

      database = forecast_database
    )
  )
}