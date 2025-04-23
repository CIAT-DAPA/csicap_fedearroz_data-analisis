## ------------------------------------------ ##
## CSICAP - Fill monthly gaps with satellite
## By: Harold Achicanoy
## Alliance Bioversity CIAT
## March 2025
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
library(pacman)
pacman::p_load(terra, tidyverse, lubridate, MASS, RobustLinearReg, future, furrr)

robust_r2_score <- function(obs, pred) {
  # Formulation of a robust r2. in this case, we have a metric that tells us how
  # much better our model is than a median of our observations using median squared
  # deviation rather than variance.
  #
  # Args:
  #   obs: a sequence of observations
  #   pred: a sequence of expected values from a model
  
  # median squared deviation total
  msdtot <- function(x) {
    median((x - median(x))^2)
  }
  
  # median squared deviation err (or residual)
  msderr <- function(x, x_hat) {
    median((x - x_hat)^2)
  }
  
  obs <- as.numeric(obs)
  prd <- as.numeric(pred)
  
  msd_res <- msderr(obs, prd)
  msd_tot <- msdtot(obs)
  
  return(1 - (msd_res/msd_tot))
}
rMSE <- function(obs, pred){sqrt(mean((obs - pred)^2))}

root <- '//CATALOGUE/WFP_ClimateRiskPr1'
outd <- paste0(root,'/7.Results/CSICAP')
datd <- paste0(root,'/1.Data')

## Load observed data ----
### IDEAM data ----
ideam_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideam_prec_monthly.parquet')) |> base::as.data.frame()
ideam_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideam_feat_stations.parquet')) |> base::as.data.frame()
ideam_unq <- ideam_unq[,c(c('longitude','latitude'),base::setdiff(names(ideam_unq),c('longitude','latitude')))]
names(ideam_unq)[3] <- 'station'
names(ideam_mnt)[2] <- 'station'

### Fedearroz data ----
fdrrz_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly.parquet')) |> base::as.data.frame()
fdrrz_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/fedearroz/fedearroz_feat_stations.parquet')) |> base::as.data.frame()

### Cenicafe data ----
cencf_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly.parquet')) |> base::as.data.frame()
cencf_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/cenicafe/cenicafe_feat_stations.parquet')) |> base::as.data.frame()

### IDEAM new data ----
ideamnew_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideamnew_prec_monthly.parquet')) |> base::as.data.frame()
ideamnew_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideamnew_feat_stations.parquet')) |> base::as.data.frame()

## Load quality index by station ----
ideam_fdrrz_evaluation_mth <- utils::read.csv(file.path(outd,'results/fedearroz/fedearroz_calidad_estaciones.csv'), fileEncoding = 'latin1')
ideam_cencf_evaluation_mth <- utils::read.csv(file.path(outd,'results/cenicafe/cenicafe_calidad_estaciones.csv'), fileEncoding = 'latin1')
ideam_ideamnew_evaluation_mth <- utils::read.csv(file.path(outd,'results/ideam/ideam_full_calidad_estaciones.csv'), fileEncoding = 'latin1')

ideam_quality <- ideam_fdrrz_evaluation_mth[ideam_fdrrz_evaluation_mth$station %in% ideam_unq$station,]
fdrrz_quality <- ideam_fdrrz_evaluation_mth[ideam_fdrrz_evaluation_mth$station %in% fdrrz_unq$station,]; rm(ideam_fdrrz_evaluation_mth)
cencf_quality <- ideam_cencf_evaluation_mth[ideam_cencf_evaluation_mth$station %in% cencf_unq$station,]; rm(ideam_cencf_evaluation_mth)
ideamnew_quality <- ideam_ideamnew_evaluation_mth[ideam_ideamnew_evaluation_mth$station %in% ideamnew_unq$station,]; rm(ideam_ideamnew_evaluation_mth)

ideam_mnt$fecha <- as.Date(ideam_mnt$fecha)
fdrrz_mnt$fecha <- as.Date(fdrrz_mnt$fecha)
cencf_mnt$fecha <- as.Date(cencf_mnt$fecha)
ideamnew_mnt$fecha <- as.Date(ideamnew_mnt$fecha)

## Monthly data imputation ----
# Function to fill monthly precipitation gaps using their corresponding satellite source
impute_monthly_missing <- function(station, qlt_station, crd_station, unq_station, mnt_station) {
  
  mnt_station <- mnt_station
  dates <- range(mnt_station$fecha)
  dates_seq <- seq(from = as.Date(dates[1]), to = as.Date(dates[2]), by = 'month')
  if (nrow(mnt_station) < length(dates_seq) | any(is.na(mnt_station$valor_observado))) {
    # Identify temporal missing data
    mss_dts <- base::setdiff(as.character(dates_seq), as.character(mnt_station$fecha))
    if (length(mss_dts) > 0) {
      aux <- data.frame(
        fecha = as.Date(mss_dts),
        station = station,
        valor_observado = NA
      ); rm(mss_dts)
      cat(paste0(station,' reports: ',nrow(aux),' missing data ... imputing\n'))
      mnt_station <- rbind(mnt_station, aux); rm(aux)
      mnt_station <- mnt_station |> dplyr::arrange(fecha)
    } else {
      cat(paste0(station,' reports: ',sum(is.na(mnt_station$valor_observado)),' missing data ... imputing\n'))
      mnt_station <- mnt_station |> dplyr::arrange(fecha)
    }
    # Identify best satellite data source according to the reported data
    best_source <- qlt_station[qlt_station$station == station,'best_source']
    # Source directory
    srcd <- dplyr::case_when(best_source == 'AgERA5' ~ paste0(datd,'/monthly_AgERA5/precipitation_flux'),
                             best_source == 'CHIRPS' ~ paste0(datd,'/monthly_CHIRPS'),
                             best_source == 'MSWEP' ~ paste0(datd,'/monthly_MSWEP'),
                             best_source == 'IMERG' ~ paste0(datd,'/monthly_IMERG'))
    
    fls <- list.files(path = srcd, pattern = ifelse(best_source %in% c('AgERA5','CHIRPS','IMERG'),'.tif$','.nc$'), full.names = T)
    if (best_source == 'MSWEP') { fls <- fls[-grep(pattern = 'new', x = fls)] }
    r <- terra::rast(fls)
    r_vls <- terra::extract(x = r, y = crd_station)
    src_dts <- dplyr::case_when(best_source == 'AgERA5' ~ gsub('.tif','',gsub('prec_','',basename(fls))),
                                best_source == 'CHIRPS' ~ gsub('.','-',paste0(gsub('.tif','',gsub('chirps-v2.0.','',basename(fls))),'.01'), fixed = T),
                                best_source == 'MSWEP' ~ as.character(as.Date(paste0(gsub('.nc','',basename(fls)),'01'),, format = "%Y%m%d")),
                                best_source == 'IMERG' ~ as.character(as.Date(str_sub(string = basename(fls), start = 1, end = 8), format = "%Y%m%d")))
    src_dfm <- data.frame(fecha = src_dts, valor_estimado = as.numeric(r_vls)[-1])
    src_dfm$fecha <- as.Date(src_dfm$fecha)
    
    mnt_station_mrg <- dplyr::left_join(x = mnt_station, y = src_dfm, by = 'fecha')
    if (best_source == 'IMERG') {
      mnt_station_mrg$valor_estimado <- mnt_station_mrg$valor_estimado * 730.5
    }
    mnt_station_mrg$year <- lubridate::year(mnt_station_mrg$fecha)
    mnt_station_mrg$month <- lubridate::month(mnt_station_mrg$fecha) |> factor(ordered = T)
    
    # mnt_station_mrg |>
    #   # dplyr::filter(month == 5) |>
    #   ggplot2::ggplot(aes(x = valor_estimado, y = valor_observado)) +
    #   ggplot2::geom_point() +
    #   ggplot2::stat_smooth(method = lm, se = T) +
    #   ggplot2::theme_minimal()
    
    # OLS
    # lm_fit <- lm(valor_observado ~ valor_estimado + month, data = mnt_station_mrg[complete.cases(mnt_station_mrg),])
    # mnt_station_mrg$linear <- mnt_station_mrg$valor_observado
    # mnt_station_mrg$linear[is.na(mnt_station_mrg$valor_observado)] <- predict(lm_fit, mnt_station_mrg)[is.na(mnt_station_mrg$valor_observado)] |> as.numeric() |> round(1)
    
    # # Theil-Sen linear regression
    # ts_fit <- mblm::mblm(valor_observado ~ valor_estimado, dataframe = mnt_station_mrg[complete.cases(mnt_station_mrg),])
    # mnt_station_mrg$theil_sen <- mnt_station_mrg$valor_observado
    # mnt_station_mrg$theil_sen[is.na(mnt_station_mrg$valor_observado)] <- predict(ts_fit, mnt_station_mrg)[is.na(mnt_station_mrg$valor_observado)] |> as.numeric() |> round(1)
    
    # Robust linear regression
    mnt_complete <- mnt_station_mrg[complete.cases(mnt_station_mrg),]
    tmp_cvr <- nrow(mnt_complete)/12 # Temporal match in years
    if (tmp_cvr > 2){ # Impute only if temporal match is higher than 2 years
      # Imputing missing data
      rlm_fit <- MASS::rlm(valor_observado ~ valor_estimado + month, data = mnt_complete); mnt_complete
      mnt_station_mrg$robust <- mnt_station_mrg$valor_observado
      mnt_station_mrg$robust[is.na(mnt_station_mrg$valor_observado)] <- predict(rlm_fit, mnt_station_mrg)[is.na(mnt_station_mrg$valor_observado)] |> as.numeric() |> round(1)
      if (any(mnt_station_mrg$robust < 0, na.rm = T)) {
        mnt_station_mrg$robust[mnt_station_mrg$robust < 0] <- 0
      }
      # Adding attributes
      mnt_station_mrg$best_source <- best_source
      mnt_station_mrg$status <- 'actual'
      mnt_station_mrg$status[is.na(mnt_station_mrg$valor_observado) & !is.na(mnt_station_mrg$robust)] <- 'imputed'
      mnt_station_mrg$status[is.na(mnt_station_mrg$valor_observado) & is.na(mnt_station_mrg$robust)] <- 'not_imputed'
      # Validation
      aux <- mnt_station_mrg[complete.cases(mnt_station_mrg),]; rownames(aux) <- 1:nrow(aux)
      n <- nrow(mnt_station_mrg) # Number of observations
      n_m <- sum(is.na(mnt_station_mrg$valor_observado)) # Number of missing
      p <- round(n_m/n, 2) # Percentage of missing
      set.seed(1235)
      seeds <- sample(x = 1:100000, size = 2000, replace = F)
      cross_validation_results <- 1:length(seeds) |>
        purrr::map(.f = function (i) {
          if (n_m > 10) {
            # If number of missing values is higher than 10 observations then
            # use percentage of missing to determine artificial missing.
            # Define training and validation datasets
            set.seed(seeds[i])
            train_ids <- sample(x = 1:nrow(aux), size = floor((1-p)*nrow(aux)))
            val_ids <- base::setdiff(1:nrow(aux), train_ids)
          } else {
            # Otherwise define the minimum number of missing as 10 observations.
            set.seed(seeds[i])
            train_ids <- sample(x = 1:nrow(aux), size = nrow(aux) - 10)
            val_ids <- base::setdiff(1:nrow(aux), train_ids)
          }
          # Train a robust regression over the training dataset
          res_dfm <- tryCatch(expr = {
            rlm_fit_train <- MASS::rlm(valor_observado ~ valor_estimado + month, data = aux[train_ids,])
            val_preds <- predict(rlm_fit_train, aux[val_ids,])
            val_preds[val_preds < 0] <- 0
            res <- data.frame(seed = seeds[i],
                              rmse = rMSE(aux$valor_observado[val_ids], val_preds),
                              pearson = cor(aux$valor_observado[val_ids], val_preds))
          },
          error = function(e) {
            res <- data.frame(seed = seeds[i],
                              rmse = NA,
                              pearson = NA)
          })
          return(res_dfm)
        }) |>
        dplyr::bind_rows()
      station_val_metrics <- data.frame(station = station, length_ts = n, missing = n_m, p_missing = p, model_r2 = robust_r2_score(obs = mnt_station_mrg$valor_observado[complete.cases(mnt_station_mrg)], pred = rlm_fit$fitted.values))
      station_val_metrics <- dplyr::bind_cols(station_val_metrics,
                                              data.frame(t(psych::describe(x = cross_validation_results[,-1])[,'median'])))
      names(station_val_metrics)[(ncol(station_val_metrics)-1):ncol(station_val_metrics)] <- c('impt_rmse', 'impt_pearson')
      station_val_metrics$temporal_match <- tmp_cvr
      station_val_metrics$median_year <- median(mnt_complete$year)
    } else {
      cat(paste0(station,' cannot be imputed, due to poor temporal coverage just ',round(tmp_cvr,2),' years\n'))
      mnt_station_mrg$robust <- NA
      mnt_station_mrg$best_source <- best_source
      mnt_station_mrg$status <- 'not imputed'
      station_val_metrics <- data.frame(station = station, length_ts = nrow(mnt_station), missing = NA, p_missing = NA, model_r2 = NA, impt_rmse = NA, impt_pearson = NA, temporal_match = tmp_cvr, median_year = median(mnt_complete$year))
    }
    
  } else {
    cat(paste0(station,' does not report missing data.\n'))
    mnt_station_mrg <- mnt_station
    mnt_station_mrg$valor_estimado <- NA
    mnt_station_mrg$year <- lubridate::year(mnt_station_mrg$fecha)
    mnt_station_mrg$month <- lubridate::month(mnt_station_mrg$fecha)
    mnt_station_mrg$robust <- NA
    mnt_station_mrg$best_source <- 'Not required'
    mnt_station_mrg$status <- 'actual'
    
    station_val_metrics <- data.frame(station = station, length_ts = nrow(mnt_station), missing = 0, p_missing = 0, model_r2 = NA, impt_rmse = NA, impt_pearson = NA, temporal_match = NA, median_year = NA)
  }
  
  mnt_station_mrg$month <- factor(x = mnt_station_mrg$month, levels = 1:12, ordered = T)
  return(list(imputed_values = mnt_station_mrg, val_metrics = station_val_metrics))
  
}

# IDEAM
if (.Platform$OS.type == "windows") {
  future::plan(multisession, workers = 10, gc = TRUE)
} else {
  future::plan(multicore, workers = 10)
}
stations <- ideam_quality$station
ideam_mnt_imputation <- 1:length(stations) |>
  furrr::future_map(.f = function(i) {
    
    station <- stations[i]
    mnt_station <- ideam_mnt
    unq_station <- ideam_unq
    qlt_station <- ideam_quality
    mnt_station <- mnt_station[mnt_station$station == station,]
    crd_station <- unq_station[unq_station$station == station, 1:2]
    res <- impute_monthly_missing(station, qlt_station, crd_station, unq_station, mnt_station)
    return(res)
    
  }, .progress = T)
future:::ClusterRegistry('stop')
future::plan(sequential)
ideam_imputed_values <- ideam_mnt_imputation |> purrr::map('imputed_values') |> dplyr::bind_rows()
ideam_val_metrics <- ideam_mnt_imputation |> purrr::map('val_metrics') |> dplyr::bind_rows()
arrow::write_parquet(x = ideam_imputed_values, sink = file.path(outd,'results/ideam/ideam_imputados_estaciones.parquet'), version = 'latest')
arrow::write_parquet(x = ideam_val_metrics, sink = file.path(outd,'results/ideam/ideam_metricas_imputacion.parquet'), version = 'latest')

# Fedearroz
if (.Platform$OS.type == "windows") {
  future::plan(multisession, workers = 10, gc = TRUE)
} else {
  future::plan(multicore, workers = 10)
}
stations <- fdrrz_quality$station
fdrrz_mnt_imputation <- 1:length(stations) |>
  furrr::future_map(.f = function(i) {
    
    station <- stations[i]
    mnt_station <- fdrrz_mnt
    unq_station <- fdrrz_unq
    qlt_station <- fdrrz_quality
    mnt_station <- mnt_station[mnt_station$station == station,]
    crd_station <- unq_station[unq_station$station == station, 1:2]
    res <- impute_monthly_missing(station, qlt_station, crd_station, unq_station, mnt_station)
    return(res)
    
  }, .progress = T)
future:::ClusterRegistry('stop')
future::plan(sequential)
fdrrz_imputed_values <- fdrrz_mnt_imputation |> purrr::map('imputed_values') |> dplyr::bind_rows()
fdrrz_val_metrics <- fdrrz_mnt_imputation |> purrr::map('val_metrics') |> dplyr::bind_rows()
arrow::write_parquet(x = fdrrz_imputed_values, sink = file.path(outd,'results/fedearroz/fedearroz_imputados_estaciones.parquet'), version = 'latest')
arrow::write_parquet(x = fdrrz_val_metrics, sink = file.path(outd,'results/fedearroz/fedearroz_metricas_imputacion.parquet'), version = 'latest')

# Cenicafe
if (.Platform$OS.type == "windows") {
  future::plan(multisession, workers = 10, gc = TRUE)
} else {
  future::plan(multicore, workers = 10)
}
stations <- cencf_quality$station
cencf_mnt_imputation <- 1:length(stations) |>
  furrr::future_map(.f = function(i) {
    
    station <- stations[i]
    mnt_station <- cencf_mnt
    unq_station <- cencf_unq
    qlt_station <- cencf_quality
    mnt_station <- mnt_station[mnt_station$station == station,]
    crd_station <- unq_station[unq_station$station == station, 1:2]
    res <- impute_monthly_missing(station, qlt_station, crd_station, unq_station, mnt_station)
    return(res)
    
  }, .progress = T)
future:::ClusterRegistry('stop')
future::plan(sequential)
cencf_imputed_values <- cencf_mnt_imputation |> purrr::map('imputed_values') |> dplyr::bind_rows()
cencf_val_metrics <- cencf_mnt_imputation |> purrr::map('val_metrics') |> dplyr::bind_rows()
arrow::write_parquet(x = cencf_imputed_values, sink = file.path(outd,'results/cenicafe/cenicafe_imputados_estaciones.parquet'), version = 'latest')
arrow::write_parquet(x = cencf_val_metrics, sink = file.path(outd,'results/cenicafe/cenicafe_metricas_imputacion.parquet'), version = 'latest')

# IDEAM new
if (.Platform$OS.type == "windows") {
  future::plan(multisession, workers = 10, gc = TRUE)
} else {
  future::plan(multicore, workers = 10)
}
stations <- ideamnew_quality$station
ideamnew_mnt_imputation <- 1:length(stations) |>
  furrr::future_map(.f = function(i) {
    
    station <- stations[i]
    mnt_station <- ideamnew_mnt
    unq_station <- ideamnew_unq
    qlt_station <- ideamnew_quality
    mnt_station <- mnt_station[mnt_station$station == station,]
    crd_station <- unq_station[unq_station$station == station, 1:2]
    res <- impute_monthly_missing(station, qlt_station, crd_station, unq_station, mnt_station)
    return(res)
    
  }, .progress = T)
future:::ClusterRegistry('stop')
future::plan(sequential)
ideamnew_imputed_values <- ideamnew_mnt_imputation |> purrr::map('imputed_values') |> dplyr::bind_rows()
ideamnew_val_metrics <- ideamnew_mnt_imputation |> purrr::map('val_metrics') |> dplyr::bind_rows()
arrow::write_parquet(x = ideamnew_imputed_values, sink = file.path(outd,'results/ideam/ideamnew_imputados_estaciones.parquet'), version = 'latest')
arrow::write_parquet(x = ideamnew_val_metrics, sink = file.path(outd,'results/ideam/ideamnew_metricas_imputacion.parquet'), version = 'latest')
