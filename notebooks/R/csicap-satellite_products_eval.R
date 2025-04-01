## ------------------------------------------ ##
## CSICAP - Satellite products evaluation
## By: Harold Achicanoy
## Alliance Bioversity CIAT
## March 2025
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)}
pacman::p_load(terra, arrow, googledrive, cloudml, data.table, tidyverse,
               lubridate, Metrics, FactoMineR, factoextra, scales, minerva,
               readxl, geodata, sf)
source('https://raw.githubusercontent.com/haachicanoy/r_scripts/refs/heads/master/calculate_index_by_pca.R')

## Define directories ----
root <- '//CATALOGUE/WFP_ClimateRiskPr1'
outd <- paste0(root,'/7.Results/CSICAP')
# root <- 'C:/Users/haachicanoy/Downloads'
# outd <- root

## Load observed data ----
### IDEAM data ----
if (file.exists(file.path(outd,'data/processed/ideam/ideam_prec_monthly.parquet'))|file.exists(file.path(outd,'data/processed/ideam/ideam_feat_stations.parquet'))) {
  ideam_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideam_prec_monthly.parquet')) |> base::as.data.frame()
  ideam_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideam_feat_stations.parquet')) |> base::as.data.frame()
} else {
  # Loading IDEAM data
  ideam <- arrow::read_parquet(file = file.path(outd,'data/raw/ideam/df_final_georeferencias.parquet'))
  ideam <- ideam |> dplyr::arrange(codigo_estacion, fecha)
  # IDEAM data into data.table format
  ideamDT <- data.table::setDT(x = ideam)
  ideamDT[,month:=lubridate::month(fecha)]
  ideamDT[,year:=lubridate::year(fecha)]
  ideamDT[,key:=paste0(year,'-',sprintf('%02d',month),'-01','__',codigo_estacion)]
  # Compute monthly values
  system.time(expr = { ideamDT_mnt <- ideamDT[,sum(valor_observado), by = key] })
  ideam_mnt <- base::as.data.frame(ideamDT_mnt); rm(ideamDT_mnt)
  ideam_mnt <- ideam_mnt |>
    tidyr::separate(col = key, into = c('fecha','codigo_estacion'), sep = '__', remove = F)
  ideam_mnt$codigo_estacion <- as.character(ideam_mnt$codigo_estacion)
  ideam_mnt$fecha <- as.Date(ideam_mnt$fecha)
  names(ideam_mnt)[ncol(ideam_mnt)] <- 'valor_observado'
  ideam_mnt$key <- NULL
  # Unique station's features
  ideam_unq <- ideam |>
    dplyr::select(codigo_estacion,etiqueta_variale,descripcion_variable,
                  frecuencia_datos,nombre,departamento,categoria,
                  latitude,longitude,altitud) |>
    unique() |>
    base::as.data.frame()
  ideam_unq$longitude <- as.numeric(ideam_unq$longitude)
  ideam_unq$latitude <- as.numeric(ideam_unq$latitude)
  ideam_unq$altitud <- as.numeric(ideam_unq$altitud)
  # Save monthly data
  arrow::write_parquet(x = ideam_mnt, sink = file.path(outd,'data/processed/ideam/ideam_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = ideam_unq, sink = file.path(outd,'data/processed/ideam/ideam_feat_stations.parquet'), version = 'latest')
  rm(ideam, ideamDT); gc(T)
}
ideam_unq <- ideam_unq[,c(c('longitude','latitude'),base::setdiff(names(ideam_unq),c('longitude','latitude')))]
names(ideam_unq)[3] <- 'station'
names(ideam_mnt)[2] <- 'station'

### Fedearroz data ----
if (file.exists(file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly.parquet'))|file.exists(file.path(outd,'data/processed/fedearroz/fedearroz_feat_stations.parquet'))) {
  fdrrz_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly.parquet')) |> base::as.data.frame()
  fdrrz_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/fedearroz/fedearroz_feat_stations.parquet')) |> base::as.data.frame()
} else {
  fdrrz_mnt <- readxl::read_excel(path = file.path(outd,'data/raw/fedearroz/prec_monthly_obs_sat.xlsx'), sheet = 1) |> base::as.data.frame()
  fdrrz_mnt <- fdrrz_mnt[,-1]
  fdrrz_mnt <- fdrrz_mnt[fdrrz_mnt$fuente == 'fedearroz',]; rownames(fdrrz_mnt) <- 1:nrow(fdrrz_mnt)
  fdrrz_mnt$fuente <- NULL
  fdrrz_unq <- fdrrz_mnt[,c('longitud','latitud','station','dpto','mun')] |> unique()
  fdrrz_mnt <- fdrrz_mnt[,c('month_year','station','prec_month')]
  names(fdrrz_mnt) <- c('fecha','station','valor_observado')
  fdrrz_mnt$fecha <- as.Date(paste0(fdrrz_mnt$fecha,'-01'))
  arrow::write_parquet(x = fdrrz_mnt, sink = file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = fdrrz_unq, sink = file.path(outd,'data/processed/fedearroz/fedearroz_feat_stations.parquet'), version = 'latest')
}

### Cenicafe data ----
if (file.exists(file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly.parquet'))|file.exists(file.path(outd,'data/processed/cenicafe/cenicafe_feat_stations.parquet'))) {
  cencf_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly.parquet')) |> base::as.data.frame()
  cencf_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/cenicafe/cenicafe_feat_stations.parquet')) |> base::as.data.frame()
} else {
  cencf <- data.table::fread(file.path(outd,'data/raw/cenicafe/Cenicafe_lluvia-diaria_20240806.csv'), encoding = 'Latin-1') |> base::as.data.frame()
  cencf <- cencf |> dplyr::arrange(COD_ESTACION, FECHA)
  cencf$VALOR[which(cencf$VALOR == '')] <- NA
  cencf$VALOR <- gsub(pattern = ',',replacement = '.',x = cencf$VALOR) |> as.numeric()
  # IDEAM data into data.table format
  cencfDT <- data.table::setDT(x = cencf)
  cencfDT[,month:=lubridate::month(FECHA)]
  cencfDT[,year:=lubridate::year(FECHA)]
  cencfDT[,key:=paste0(year,'-',sprintf('%02d',month),'-01','__',COD_ESTACION)]
  # Compute monthly values
  system.time(expr = {
    cencfDT_mnt <- cencfDT[,valor_observado:=sum(VALOR, na.rm = T), by = key]
    cencfDT_NAs <- cencfDT[,faltantes:=sum(is.na(VALOR)), by = key]
  }); rm(cencfDT_NAs)
  cencf_mnt <- base::as.data.frame(cencfDT_mnt); rm(cencfDT_mnt)
  cencf_mnt$valor_observado[which(cencf_mnt$faltantes >= 4)] <- NA
  cencf_mnt <- cencf_mnt |>
    dplyr::select(key, valor_observado, faltantes) |>
    base::unique() |>
    tidyr::separate(col = key, into = c('fecha','station'), sep = '__', remove = F)
  cencf_mnt$station <- as.character(cencf_mnt$station)
  cencf_mnt$fecha <- as.Date(cencf_mnt$fecha)
  cencf_mnt$key <- NULL
  cencf_mnt <- cencf_mnt |>
    dplyr::arrange(station, fecha) |>
    base::as.data.frame()
  cencf_mnt$faltantes <- NULL
  # Unique station's features
  cencf_unq <- cencf |>
    dplyr::select(COD_ESTACION,NOMBRE,DEPARTAMENTO,MUNICIPIO,
                  TES_CODIGO,NOM_VARIABLE) |>
    base::unique() |>
    base::as.data.frame()
  names(cencf_unq) <- c('station','nombre','departamento','municipio','tes_codigo','nom_variable')
  cencf_unq$station <- as.character(cencf_unq$station)
  cencf_crd <- utils::read.csv(file.path(outd,'data/raw/cenicafe/Catalogo_precipitacion.csv'))
  names(cencf_crd) <- tolower(names(cencf_crd))
  names(cencf_crd)[1] <- 'station'
  cencf_crd$station <- as.character(cencf_crd$station)
  cencf_crd <- cencf_crd[cencf_crd$nom_variable == 'Precipitación acumulada día',]
  
  cencf_unq <- dplyr::left_join(x = cencf_unq,
                                y = cencf_crd |> dplyr::select(station, nombre, longitud, latitud, altura, corriente, estado, cat.inicial, cat.final, fech_inicio, fech_ult, años, n_registros_diarios),
                                by = c('station','nombre'))
  rm(cencf_crd)
  names(cencf_unq)[9] <- 'altitud'
  
  cencf_unq <- cencf_unq[,c('longitud','latitud','station',base::setdiff(names(cencf_unq),c('longitud','latitud','station')))] |> unique()
  arrow::write_parquet(x = cencf_mnt, sink = file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = cencf_unq, sink = file.path(outd,'data/processed/cenicafe/cenicafe_feat_stations.parquet'), version = 'latest')
}

### IDEAM new data ----
if (file.exists(file.path(outd,'data/processed/ideam/ideamnew_prec_monthly.parquet'))|file.exists(file.path(outd,'data/processed/ideam/ideamnew_feat_stations.parquet'))) {
  ideamnew_mnt <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideamnew_prec_monthly.parquet')) |> base::as.data.frame()
  ideamnew_unq <- arrow::read_parquet(file = file.path(outd,'data/processed/ideam/ideamnew_feat_stations.parquet')) |> base::as.data.frame()
} else {
  # Loading IDEAM data
  ideamnew <- utils::read.csv(file.path(outd,'data/raw/ideam/missing_station_filter.csv'))
  ideamnew <- ideamnew[,3:ncol(ideamnew)]
  ideamnew$fecha <- as.Date(ideamnew$fecha)
  ideamnew <- ideamnew |> dplyr::arrange(codigo, fecha)
  ideamnew$month <- lubridate::month(ideamnew$fecha)
  ideamnew$year <- lubridate::year(ideamnew$fecha)
  ideamnew$key <- paste0(ideamnew$year,'-',sprintf('%02d',ideamnew$month),'-01','__',ideamnew$codigo)
  # Compute monthly values
  ideamnew_mnt <- ideamnew |>
    dplyr::select(key, precipitacion) |>
    dplyr::group_by(key) |>
    dplyr::summarise(valor_observado = sum(precipitacion)) |>
    dplyr::ungroup() |>
    base::as.data.frame()
  ideamnew_mnt <- ideamnew_mnt |>
    tidyr::separate(col = key, into = c('fecha','codigo_estacion'), sep = '__', remove = F)
  ideamnew_mnt$codigo_estacion <- as.character(ideamnew_mnt$codigo_estacion)
  ideamnew_mnt$fecha <- as.Date(ideamnew_mnt$fecha)
  ideamnew_mnt$key <- NULL
  names(ideamnew_mnt)[2] <- 'station'
  # Unique station's features
  ideamnew_unq <- ideamnew |>
    dplyr::select(codigo,nombre,departamento,municipio,
                  latitud,longitud,altitud) |>
    unique() |>
    base::as.data.frame()
  ideamnew_unq$longitude <- as.numeric(ideamnew_unq$longitud); ideamnew_unq$longitud <- NULL
  ideamnew_unq$latitude <- as.numeric(ideamnew_unq$latitud); ideamnew_unq$latitud <- NULL
  ideamnew_unq$altitude <- as.numeric(ideamnew_unq$altitud); ideamnew_unq$altitud <- NULL
  names(ideamnew_unq)[1] <- 'station'
  ideamnew_unq$station <- as.character(ideamnew_unq$station)
  ideamnew_unq <- ideamnew_unq[,c('longitude','latitude','station','nombre','departamento','municipio','altitude')]
  # Save monthly data
  arrow::write_parquet(x = ideamnew_mnt, sink = file.path(outd,'data/processed/ideam/ideamnew_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = ideamnew_unq, sink = file.path(outd,'data/processed/ideam/ideamnew_feat_stations.parquet'), version = 'latest')
  rm(ideamnew); gc(F,T,T)
}

# Mapping altitude
ideam_unq |>
  ggplot2::ggplot(aes(x = longitude, y = latitude, color = altitud)) +
  ggplot2::geom_point() +
  ggplot2::coord_equal() +
  ggplot2::scale_color_gradientn(colours = rainbow(5)) +
  ggplot2::theme_minimal()

## Load satellite data ----
extract_satellite_data <- function(crds = ideam_unq, Data = ideam_mnt){
  
  stlls <- c('CHIRPS','AgERA5','MSWEP','IMERG')
  stlls_data <- 1:length(stlls) |>
    purrr::map(.f = function(i){
      stll <- stlls[i]
      dtad <- dplyr::case_when( # Data directory
        stll == 'CHIRPS' ~ paste0(root,'/1.Data/monthly_CHIRPS'),
        stll == 'AgERA5' ~ paste0(root,'/1.Data/monthly_AgERA5/precipitation_flux'),
        stll == 'MSWEP' ~ paste0(root,'/1.Data/monthly_MSWEPcorrected'),
        stll == 'IMERG' ~ paste0(root,'/1.Data/monthly_IMERG')
      )
      frmt <- ifelse(stll %in% c('CHIRPS','AgERA5','IMERG'), yes = '.tif$', no = '.nc$')
      fls <- list.files(path = dtad, pattern = frmt, full.names = T)
      src <- terra::rast(fls)
      names(src) <- dplyr::case_when( # Layer names
        stll == 'CHIRPS' ~ paste0('D',paste0(gsub('.','-',gsub('.tif','', gsub('chirps-v2.0.','',basename(fls))), fixed = T),'-01')),
        stll == 'AgERA5' ~ paste0('D',gsub('.tif','',gsub('prec_','',basename(fls)))),
        stll == 'MSWEP' ~ paste0('D',as.character(lubridate::ymd(gsub('.nc','01',basename(fls))))),
        stll == 'IMERG' ~ paste0('D',as.character(lubridate::ymd(paste0(gsub('000000_11.132km.tif','',basename(fls))))))
      )
      src_mnt <- terra::extract(x = src, y = crds[,1:2])
      src_mnt$station <- crds$station
      src_mnt$ID <- NULL
      src_mnt <- src_mnt |>
        tidyr::pivot_longer(cols = 1:(ncol(src_mnt)-1), names_to = 'fecha', values_to = stll) |>
        base::as.data.frame()
      src_mnt$fecha <- as.Date(gsub('D','',src_mnt$fecha))
      
      return(src_mnt)
      
    })
  
  tbls <- list()
  tbls[[1]] <- Data
  tbls[2:5] <- stlls_data
  
  Data_final <- tbls |> purrr::reduce(dplyr::left_join, by = c('fecha','station'))
  return(Data_final)
}
# Get satellite data for IDEAM
if(!file.exists(file.path(outd,'data/processed/ideam/ideam_prec_monthly_merged.parquet'))){
  ideam_mrg <- extract_satellite_data(crds = ideam_unq, Data = ideam_mnt)
  ideam_mrg$IMERG <- ideam_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  arrow::write_parquet(x = ideam_mrg, sink = file.path(outd,'ideam_prec_monthly_merged.parquet'), version = 'latest')
} else {
  ideam_mrg <- arrow::read_parquet(file.path(outd,'data/processed/ideam/ideam_prec_monthly_merged.parquet'))
}
# Get satellite data for Fedearroz
if(!file.exists(file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly_merged.parquet'))){
  fdrrz_mrg <- extract_satellite_data(crds = fdrrz_unq, Data = fdrrz_mnt)
  fdrrz_mrg$IMERG <- fdrrz_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  arrow::write_parquet(x = fdrrz_mrg, sink = file.path(outd,'fedearroz_prec_monthly_merged.parquet'), version = 'latest')
} else {
  fdrrz_mrg <- arrow::read_parquet(file.path(outd,'data/processed/fedearroz/fedearroz_prec_monthly_merged.parquet'))
}
# Get satellite data for Cenicafe
if(!file.exists(file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly_merged.parquet'))){
  cencf_mrg <- extract_satellite_data(crds = cencf_unq, Data = cencf_mnt)
  cencf_mrg$IMERG <- cencf_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  arrow::write_parquet(x = cencf_mrg, sink = file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly_merged.parquet'), version = 'latest')
} else {
  cencf_mrg <- arrow::read_parquet(file.path(outd,'data/processed/cenicafe/cenicafe_prec_monthly_merged.parquet'))
}
# Get satellite data for IDEAM new
if(!file.exists(file.path(outd,'data/processed/ideam/ideamnew_prec_monthly_merged.parquet'))){
  ideamnew_mrg <- extract_satellite_data(crds = ideamnew_unq, Data = ideamnew_mnt)
  ideamnew_mrg$IMERG <- ideamnew_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  arrow::write_parquet(x = ideamnew_mrg, sink = file.path(outd,'data/processed/ideam/ideamnew_prec_monthly_merged.parquet'), version = 'latest')
} else {
  ideamnew_mrg <- arrow::read_parquet(file.path(outd,'data/processed/ideam/ideamnew_prec_monthly_merged.parquet'))
}

## Evaluation metrics ----
# Kling–Gupta efficiency
kling_gupta <- function(obs, pred){
  # Pearson correlation
  cc <- cor(obs, pred, method = "pearson")
  std_obs <- sd(obs)
  std_pred <- sd(pred)
  # Avoid zero division
  if (std_obs == 0 || std_pred == 0) {
    alpha <- NaN  # Ratio of standard deviation indefinite
  } else {
    alpha <- std_pred / std_obs  # Ratio of standard deviation
  }
  mean_obs <- mean(obs)
  mean_pred <- mean(pred)
  # Avoid zero division
  if (mean_obs == 0) {
    beta <- NaN  # Mean's ratio indefinite
  } else {
    beta <- mean_pred / mean_obs  # Mean's ratio
  }
  # Compute KGE with NaN validation
  if (is.na(cc) || is.na(alpha) || is.na(beta)) {
    return(NaN)
  } else {
    kge <- 1 - sqrt((cc - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
    return(kge)
  }
}
# R-squared
r2 <- function(obs, pred){cor(obs, pred)^2}
# Adjusted R-squared
r2_adj <- function(r2, n, p = 1){
  if (n > p + 1) { return(1 - ((1 - r2) * (n - 1) / (n - p - 1))) } else { return(NaN) }
}
# Root mean squared error
rMSE <- function(obs, pred){sqrt(mean((obs - pred)^2))}
# Mean absolute percentage error
mape <- function(obs, pred){
  obs_nonzero <- obs[obs != 0]
  pred_nonzero <- pred[obs != 0]
  return(mean(abs((obs_nonzero - pred_nonzero) / obs_nonzero)) * 100)
}
# Mean Arc-tangent Absolute Percentage Error
maape <- function(obs, pred){
  obs_nonzero <- obs[obs != 0]
  pred_nonzero <- pred[obs != 0]
  return(mean(atan(abs((obs_nonzero - pred_nonzero) / obs_nonzero))))
}
# Standard error
std_error <- function(obs, pred){
  diff <- obs - pred
  if(all(is.na(diff))){
    return(NA)
  } else {
    stde <- sd(diff, na.rm = T)
    return(stde)
  }
}
# Bias
bias <- function(obs, pred){ mean(pred) -  mean(obs) }

## Evaluate satellite data ----
stllts_evaluation <- function(Data = ideam_mrg, analysis = 'general_mth'){
  
  Data$year <- paste0('Y',lubridate::year(Data$fecha))
  Data$month <- lubridate::month(Data$fecha, label = T, abbr = T)
  Data$quarter <- lubridate::quarter(Data$fecha)
  Data$quarter <- dplyr::case_when(Data$quarter == 1 ~ 'EFM',
                                   Data$quarter == 2 ~ 'AMJ',
                                   Data$quarter == 3 ~ 'JAS',
                                   Data$quarter == 4 ~ 'OND')
  
  # Satellite sources
  stlls <- c('CHIRPS','AgERA5','MSWEP','IMERG')
  
  # Data preparation for quarterly data
  if(analysis == 'general_qrt' | analysis == 'per_qrt'){
    Data <- Data |>
      dplyr::select(station,valor_observado,
                    CHIRPS,AgERA5,MSWEP,IMERG,
                    year, quarter) |>
      dplyr::group_by(station,year,quarter) |>
      dplyr::summarise_all(sum) |>
      dplyr::ungroup() |> base::as.data.frame()
  }
  
  # Stations IDs
  stations <- unique(Data$station)
  
  if(analysis == 'general_mth' | analysis == 'general_qrt'){
    
    all_metrics <- 1:length(stlls) |>
      purrr::map(.f = function(s){
        stll <- stlls[s]
        cat('Evaluating',stll,'data source. Getting',analysis,'metrics.\n')
        # Filter dataset by observed and satellite data of interest
        dfm <- Data |> dplyr::select('station','valor_observado',stll) |> tidyr::drop_na()
        # Get metrics for all stations
        metrics <- 1:length(stations) |>
          purrr::map(.f = function(i){
            cat('Comparing against',stations[i],'\n')
            aux <- dfm |>
              dplyr::select('station','valor_observado',stll) |>
              dplyr::filter(station == stations[i])
            if(nrow(aux) > 2){
              names(aux)[3] <- 'valor_registrado'
              obs <- aux$valor_observado
              prd <- aux$valor_registrado
              mtrcs <- data.frame(station = stations[i],
                                  n           = nrow(aux),
                                  rsq         = r2(obs = obs, pred = prd),
                                  rsq_adj     = r2_adj(r2 = r2(obs = obs, pred = prd), n = length(obs)),
                                  rmse        = rMSE(obs = obs, pred = prd),
                                  mae         = Metrics::mae(obs, prd),
                                  std_error   = std_error(obs = obs, pred = prd),
                                  kge         = kling_gupta(obs = obs, pred = prd),
                                  Spearman    = cor(x = obs, y = prd, method = 'spearman'),
                                  Kendall     = cor(x = obs, y = prd, method = 'kendall'),
                                  mic         = minerva::mine(x = obs, y = prd)$MIC,
                                  bias        = bias(obs = obs, pred = prd),
                                  mape        = mape(obs = obs, pred = prd),
                                  maape       = maape(obs = obs, pred = prd))
            } else {
              mtrcs <- data.frame(station = stations[i],
                                  n         = nrow(aux),
                                  rsq       = NA,
                                  rsq_adj   = NA,
                                  rmse      = NA,
                                  mae       = NA,
                                  std_error = NA,
                                  kge       = NA,
                                  Spearman  = NA,
                                  Kendall   = NA,
                                  mic       = NA,
                                  bias      = NA,
                                  mape      = NA,
                                  maape     = NA)
            }
            return(mtrcs)
          }) |>
          dplyr::bind_rows()
        metrics$source <- stll
        return(metrics)
      }) |>
      dplyr::bind_rows()
    all_metrics$analysis <- analysis
    all_metrics$filter <- 'all'
    
  } else {
    if(analysis == 'per_mth' | analysis == 'per_qrt'){
      ifelse(analysis == 'per_mth',
             {Data_lst <- Data |> dplyr::group_by(month) |> dplyr::group_split(); fltr <- 'month'},
             {Data_lst <- Data |> dplyr::group_by(quarter) |> dplyr::group_split(); fltr <- 'quarter'})
      
      all_metrics <- 1:length(Data_lst) |>
        purrr::map(.f = function(j){ # Looping through sub-dataset
          
          all_metrics <- 1:length(stlls) |>
            purrr::map(.f = function(s){
              stll <- stlls[s]
              cat('Evaluating',stll,'data source. Getting',analysis,'metrics.\n')
              # Filter dataset by observed and satellite data of interest
              dfm <- Data_lst[[j]] |> dplyr::select('station','valor_observado',stll) |> tidyr::drop_na()
              # Get metrics for all stations
              metrics <- 1:length(stations) |>
                purrr::map(.f = function(i){
                  cat('Comparing against',stations[i],'\n')
                  aux <- dfm |>
                    dplyr::select('station','valor_observado',stll) |>
                    dplyr::filter(station == stations[i])
                  if(nrow(aux) > 2){
                    names(aux)[3] <- 'valor_registrado'
                    obs <- aux$valor_observado
                    prd <- aux$valor_registrado
                    mtrcs <- data.frame(station = stations[i],
                                        n           = nrow(aux),
                                        rsq         = r2(obs = obs, pred = prd),
                                        rsq_adj     = r2_adj(r2 = r2(obs = obs, pred = prd), n = length(obs)),
                                        rmse        = rMSE(obs = obs, pred = prd),
                                        mae         = Metrics::mae(obs, prd),
                                        std_error   = std_error(obs = obs, pred = prd),
                                        kge         = kling_gupta(obs = obs, pred = prd),
                                        Spearman    = cor(x = obs, y = prd, method = 'spearman'),
                                        Kendall     = cor(x = obs, y = prd, method = 'kendall'),
                                        mic         = minerva::mine(x = obs, y = prd)$MIC,
                                        bias        = bias(obs = obs, pred = prd),
                                        mape        = mape(obs = obs, pred = prd),
                                        maape       = maape(obs = obs, pred = prd))
                  } else {
                    mtrcs <- data.frame(station = stations[i],
                                        n         = nrow(aux),
                                        rsq       = NA,
                                        rsq_adj   = NA,
                                        rmse      = NA,
                                        mae       = NA,
                                        std_error = NA,
                                        kge       = NA,
                                        Spearman  = NA,
                                        Kendall   = NA,
                                        mic       = NA,
                                        bias      = NA,
                                        mape      = NA,
                                        maape     = NA)
                  }
                  return(mtrcs)
                }) |>
                dplyr::bind_rows()
              metrics$source <- stll
              return(metrics)
            }) |>
            dplyr::bind_rows()
          all_metrics$analysis <- analysis
          all_metrics$filter <- (Data_lst[[j]] |> dplyr::pull(fltr) |> as.character() |> unique())
          return(all_metrics)
      }) |>
        dplyr::bind_rows()
    }
  }
  return(all_metrics)
}

# IDEAM
if(!file.exists(file.path(outd,'results/ideam/ideam_evaluation_year-month_general.parquet'))){
  ideam_general_mth <- stllts_evaluation(Data = ideam_mrg, analysis = 'general_mth')
  arrow::write_parquet(x = ideam_general_mth, sink = file.path(outd,'results/ideam/ideam_evaluation_year-month_general.parquet'), version = 'latest')
} else {
  ideam_general_mth <- arrow::read_parquet(file = file.path(outd,'results/ideam/ideam_evaluation_year-month_general.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideam_evaluation_year-quarter_general.parquet'))){
  ideam_general_qrt <- stllts_evaluation(Data = ideam_mrg, analysis = 'general_qrt')
  arrow::write_parquet(x = ideam_general_qrt, sink = file.path(outd,'results/ideam/ideam_evaluation_year-quarter_general.parquet'), version = 'latest')
} else {
  ideam_general_qrt <- arrow::read_parquet(file.path(outd,'results/ideam/ideam_evaluation_year-quarter_general.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideam_evaluation_per-month.parquet'))){
  ideam_per_mth <- stllts_evaluation(Data = ideam_mrg, analysis = 'per_mth')
  ideam_per_mth$filter <- factor(x = ideam_per_mth$filter, levels = unique(ideam_per_mth$filter))
  arrow::write_parquet(x = ideam_per_mth, sink = file.path(outd,'results/ideam/ideam_evaluation_per-month.parquet'), version = 'latest')
} else {
  ideam_per_mth <- arrow::read_parquet(file.path(outd,'results/ideam/ideam_evaluation_per-month.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideam_evaluation_per-quarter.parquet'))){
  ideam_per_qrt <- stllts_evaluation(Data = ideam_mrg, analysis = 'per_qrt')
  ideam_per_qrt$filter <- factor(x = ideam_per_qrt$filter, levels = c('EFM','AMJ','JAS','OND'))
  arrow::write_parquet(x = ideam_per_qrt, sink = file.path(outd,'results/ideam/ideam_evaluation_per-quarter.parquet'), version = 'latest')
} else {
  ideam_per_qrt <- arrow::read_parquet(file.path(outd,'results/ideam/ideam_evaluation_per-quarter.parquet'))
}

# Fedearroz
if(!file.exists(file.path(outd,'results/fedearroz/fedearroz_evaluation_year-month_general.parquet'))){
  fdrrz_general_mth <- stllts_evaluation(Data = fdrrz_mrg, analysis = 'general_mth')
  arrow::write_parquet(x = fdrrz_general_mth, sink = file.path(outd,'results/fedearroz/fedearroz_evaluation_year-month_general.parquet'), version = 'latest')
} else {
  fdrrz_general_mth <- arrow::read_parquet(file = file.path(outd,'results/fedearroz/fedearroz_evaluation_year-month_general.parquet'))
}
if(!file.exists(file.path(outd,'results/fedearroz/fedearroz_evaluation_year-quarter_general.parquet'))){
  fdrrz_general_qrt <- stllts_evaluation(Data = fdrrz_mrg, analysis = 'general_qrt')
  arrow::write_parquet(x = fdrrz_general_qrt, sink = file.path(outd,'results/fedearroz/fedearroz_evaluation_year-quarter_general.parquet'), version = 'latest')
} else {
  fdrrz_general_qrt <- arrow::read_parquet(file.path(outd,'results/fedearroz/fedearroz_evaluation_year-quarter_general.parquet'))
}
if(!file.exists(file.path(outd,'results/fedearroz/fedearroz_evaluation_per-month.parquet'))){
  fdrrz_per_mth <- stllts_evaluation(Data = fdrrz_mrg, analysis = 'per_mth')
  fdrrz_per_mth$filter <- factor(x = fdrrz_per_mth$filter, levels = unique(fdrrz_per_mth$filter))
  arrow::write_parquet(x = fdrrz_per_mth, sink = file.path(outd,'results/fedearroz/fedearroz_evaluation_per-month.parquet'), version = 'latest')
} else {
  fdrrz_per_mth <- arrow::read_parquet(file.path(outd,'results/fedearroz/fedearroz_evaluation_per-month.parquet'))
}
if(!file.exists(file.path(outd,'results/fedearroz/fedearroz_evaluation_per-quarter.parquet'))){
  fdrrz_per_qrt <- stllts_evaluation(Data = fdrrz_mrg, analysis = 'per_qrt')
  fdrrz_per_qrt$filter <- factor(x = fdrrz_per_qrt$filter, levels = c('EFM','AMJ','JAS','OND'))
  arrow::write_parquet(x = fdrrz_per_qrt, sink = file.path(outd,'results/fedearroz/fedearroz_evaluation_per-quarter.parquet'), version = 'latest')
} else {
  fdrrz_per_qrt <- arrow::read_parquet(file.path(outd,'results/fedearroz/fedearroz_evaluation_per-quarter.parquet'))
}

# Cenicafe
if(!file.exists(file.path(outd,'results/cenicafe/cenicafe_evaluation_year-month_general.parquet'))){
  cencf_general_mth <- stllts_evaluation(Data = cencf_mrg, analysis = 'general_mth')
  arrow::write_parquet(x = cencf_general_mth, sink = file.path(outd,'results/cenicafe/cenicafe_evaluation_year-month_general.parquet'), version = 'latest')
} else {
  cencf_general_mth <- arrow::read_parquet(file = file.path(outd,'results/cenicafe/cenicafe_evaluation_year-month_general.parquet'))
}
if(!file.exists(file.path(outd,'results/cenicafe/cenicafe_evaluation_year-quarter_general.parquet'))){
  cencf_general_qrt <- stllts_evaluation(Data = cencf_mrg, analysis = 'general_qrt')
  arrow::write_parquet(x = cencf_general_qrt, sink = file.path(outd,'results/cenicafe/cenicafe_evaluation_year-quarter_general.parquet'), version = 'latest')
} else {
  cencf_general_qrt <- arrow::read_parquet(file.path(outd,'results/cenicafe/cenicafe_evaluation_year-quarter_general.parquet'))
}
if(!file.exists(file.path(outd,'results/cenicafe/cenicafe_evaluation_per-month.parquet'))){
  cencf_per_mth <- stllts_evaluation(Data = cencf_mrg, analysis = 'per_mth')
  cencf_per_mth$filter <- factor(x = cencf_per_mth$filter, levels = unique(cencf_per_mth$filter))
  arrow::write_parquet(x = cencf_per_mth, sink = file.path(outd,'results/cenicafe/cenicafe_evaluation_per-month.parquet'), version = 'latest')
} else {
  cencf_per_mth <- arrow::read_parquet(file.path(outd,'results/cenicafe/cenicafe_evaluation_per-month.parquet'))
}
if(!file.exists(file.path(outd,'results/cenicafe/cenicafe_evaluation_per-quarter.parquet'))){
  cencf_per_qrt <- stllts_evaluation(Data = cencf_mrg, analysis = 'per_qrt')
  cencf_per_qrt$filter <- factor(x = cencf_per_qrt$filter, levels = c('EFM','AMJ','JAS','OND'))
  arrow::write_parquet(x = cencf_per_qrt, sink = file.path(outd,'results/cenicafe/cenicafe_evaluation_per-quarter.parquet'), version = 'latest')
} else {
  cencf_per_qrt <- arrow::read_parquet(file.path(outd,'results/cenicafe/cenicafe_evaluation_per-quarter.parquet'))
}

# IDEAM new
if(!file.exists(file.path(outd,'results/ideam/ideamnew_evaluation_year-month_general.parquet'))){
  ideamnew_general_mth <- stllts_evaluation(Data = ideamnew_mrg, analysis = 'general_mth')
  arrow::write_parquet(x = ideamnew_general_mth, sink = file.path(outd,'results/ideam/ideamnew_evaluation_year-month_general.parquet'), version = 'latest')
} else {
  ideamnew_general_mth <- arrow::read_parquet(file = file.path(outd,'results/ideam/ideamnew_evaluation_year-month_general.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideamnew_evaluation_year-quarter_general.parquet'))){
  ideamnew_general_qrt <- stllts_evaluation(Data = ideamnew_mrg, analysis = 'general_qrt')
  arrow::write_parquet(x = ideamnew_general_qrt, sink = file.path(outd,'results/ideam/ideamnew_evaluation_year-quarter_general.parquet'), version = 'latest')
} else {
  ideamnew_general_qrt <- arrow::read_parquet(file.path(outd,'results/ideam/ideamnew_evaluation_year-quarter_general.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideamnew_evaluation_per-month.parquet'))){
  ideamnew_per_mth <- stllts_evaluation(Data = ideamnew_mrg, analysis = 'per_mth')
  ideamnew_per_mth$filter <- factor(x = ideamnew_per_mth$filter, levels = unique(ideamnew_per_mth$filter))
  arrow::write_parquet(x = ideamnew_per_mth, sink = file.path(outd,'results/ideam/ideamnew_evaluation_per-month.parquet'), version = 'latest')
} else {
  ideamnew_per_mth <- arrow::read_parquet(file.path(outd,'results/ideam/ideamnew_evaluation_per-month.parquet'))
}
if(!file.exists(file.path(outd,'results/ideam/ideamnew_evaluation_per-quarter.parquet'))){
  ideamnew_per_qrt <- stllts_evaluation(Data = ideamnew_mrg, analysis = 'per_qrt')
  ideamnew_per_qrt$filter <- factor(x = ideamnew_per_qrt$filter, levels = c('EFM','AMJ','JAS','OND'))
  arrow::write_parquet(x = ideamnew_per_qrt, sink = file.path(outd,'results/ideam/ideamnew_evaluation_per-quarter.parquet'), version = 'latest')
} else {
  ideamnew_per_qrt <- arrow::read_parquet(file.path(outd,'results/ideam/ideamnew_evaluation_per-quarter.parquet'))
}

## Get elevation data ----
col_dem <- geodata::elevation_30s(country = 'COL', path = tempdir()) # Get Digital Elevation Model

ideam_general_mth <- ideam_general_mth |> dplyr::left_join(y = ideam_unq[,c('longitude','latitude','altitud','station')], by = 'station')
names(ideam_general_mth)[ncol(ideam_general_mth)] <- 'altitude'
ideam_general_qrt <- ideam_general_qrt |> dplyr::left_join(y = ideam_unq[,c('longitude','latitude','altitud','station')], by = 'station')
names(ideam_general_qrt)[ncol(ideam_general_qrt)] <- 'altitude'

fdrrz_general_mth <- fdrrz_general_mth |> dplyr::left_join(y = fdrrz_unq[,c('longitud','latitud','station')], by = 'station')
names(fdrrz_general_mth)[(ncol(fdrrz_general_mth)-1):ncol(fdrrz_general_mth)] <- c('longitude','latitude')
fdrrz_general_mth$altitude <- terra::extract(x = col_dem, y = fdrrz_general_mth[,c('longitude','latitude')]) |> dplyr::pull(COL_elv_msk); rm(col_dem)

cencf_general_mth <- cencf_general_mth |> dplyr::left_join(y = cencf_unq[,c('longitud','latitud','altitud','station')], by = 'station')
names(cencf_general_mth)[(ncol(cencf_general_mth)-2):ncol(cencf_general_mth)] <- c('longitude','latitude','altitude')

ideamnew_general_mth <- ideamnew_general_mth |> dplyr::left_join(y = ideamnew_unq[,c('longitude','latitude','station','altitude')], by = 'station')

## Composite index by using PCA ----
get_composite_index <- function(metrics = ideam_general_mth,
                                supplementary = fdrrz_general_mth,
                                outdir = 'results/fedearroz/figures'){
  if(!is.null(supplementary)){
    # Merge both data frames (active and supplementary)
    metrics_lst <- rbind(metrics, supplementary) |>
      dplyr::select(-analysis, -filter) |>
      dplyr::group_by(source) |>
      dplyr::group_split()
    stations_lst <- metrics_lst |> purrr::map(.f = function(dfm) dfm |> dplyr::pull(station))
    stations_lst <- stations_lst[[which.max(unlist(purrr::map(stations_lst, length)))[1]]]
    # Metrics index
    metrics_pca_idx <- metrics_lst |>
      purrr::map(.f = function(metrics_dfm){
        src <- unique(metrics_dfm$source)
        act_ids <- unique(metrics$station)
        wgs <- metrics_dfm$n[metrics_dfm$station %in% act_ids]/sum(metrics_dfm$n[metrics_dfm$station %in% act_ids])
        spl_ids <- which(!(metrics_dfm$station %in% act_ids))
        pca_fll <- metrics_dfm |> dplyr::select(-station,-source,-n) |> FactoMineR::PCA(scale.unit = T, ind.sup = spl_ids, row.w = wgs, graph = F)
        # Select variables with more than 50% of explained variability per component
        vrs <- apply(X = pca_fll$var$cos2[,1:2], MARGIN = 1, FUN = function(x){max(x) > 0.5})
        vrs <- names(vrs[vrs == T])
        metrics_dfm_aux <- metrics_dfm |> dplyr::mutate(rmse = -1*rmse, mae = -1*mae, bias = -1*bias, std_error = -1*std_error)
        metrics_dfm_vrs <- metrics_dfm_aux |> dplyr::select(vrs); rm(metrics_dfm_aux)
        pca_res <- metrics_dfm_vrs |> FactoMineR::PCA(scale.unit = T, ind.sup = spl_ids, row.w = wgs, graph = F)
        gg_pca_fll <<- factoextra::fviz_pca_var(pca_fll, col.var = 'cos2', gradient.cols = c('#00AFBB','#E7B800','#FC4E07'), repel = T)
        gg_pca_res <<- factoextra::fviz_pca_var(pca_res, col.var = 'cos2', gradient.cols = c('#00AFBB','#E7B800','#FC4E07'), repel = T)
        dir.create(path = file.path(outd,outdir), F, T)
        ggplot2::ggsave(filename = file.path(outd,paste0(outdir,'/Figure_annexes_',src,'_pca_full.png')), plot = gg_pca_fll, device = 'png', width = 6, height = 6, units = 'in', dpi = 350)
        ggplot2::ggsave(filename = file.path(outd,paste0(outdir,'/Figure_annexes_',src,'_pca_final.png')), plot = gg_pca_res, device = 'png', width = 6, height = 6, units = 'in', dpi = 350)
        res <- data.frame(index = as.numeric(index_cal(pca_res)))
        res$station <- metrics_dfm$station
        if (nrow(res) < length(stations_lst)) {
          aux <- data.frame(index = NA, station = base::setdiff(stations_lst, res$station))
          res <- rbind(res, aux); rm(aux)
        }
        names(res)[1] <- src
        return(res)
      }) |>
      purrr::reduce(dplyr::left_join, by = 'station')
    metrics_pca_idx$longitude <- metrics_lst[[1]] |> dplyr::pull(longitude)
    metrics_pca_idx$latitude  <- metrics_lst[[1]] |> dplyr::pull(latitude)
    metrics_pca_idx$altitude  <- metrics_lst[[1]] |> dplyr::pull(altitude)
    metrics_pca_idx[metrics_pca_idx == 50] <- NA
    return(metrics_pca_idx)
  } else {
    # Split data.frame per satellite source
    metrics_lst <- metrics |>
      dplyr::select(-analysis, -filter) |>
      dplyr::group_by(source) |>
      dplyr::group_split()
    # Metrics index
    metrics_pca_idx <- metrics_lst |>
      purrr::map(.f = function(metrics_dfm){
        src <- unique(metrics_dfm$source)
        wgs <- metrics_dfm$n/sum(metrics_dfm$n)
        pca_fll <- metrics_dfm |> dplyr::select(-station,-source,-n) |> FactoMineR::PCA(scale.unit = T, row.w = wgs, graph = F)
        # Select variables with more than 50% of explained variability per component
        vrs <- apply(X = pca_fll$var$cos2[,1:2], MARGIN = 1, FUN = function(x){max(x) > 0.5})
        vrs <- names(vrs[vrs == T])
        metrics_dfm_aux <- metrics_dfm |> dplyr::mutate(rmse = -1*rmse, mae = -1*mae, bias = -1*bias, std_error = -1*std_error)
        metrics_dfm_vrs <- metrics_dfm_aux |> dplyr::select(vrs); rm(metrics_dfm_aux)
        pca_res <- metrics_dfm_vrs |> FactoMineR::PCA(scale.unit = T, row.w = wgs, graph = F)
        gg_pca_fll <- factoextra::fviz_pca_var(pca_fll, col.var = 'cos2', gradient.cols = c('#00AFBB','#E7B800','#FC4E07'), repel = T)
        gg_pca_res <- factoextra::fviz_pca_var(pca_res, col.var = 'cos2', gradient.cols = c('#00AFBB','#E7B800','#FC4E07'), repel = T)
        dir.create(path = file.path(outd,outdir), F, T)
        ggplot2::ggsave(filename = file.path(outd,paste0(outdir,'/Figure_annexes_',src,'_pca_full.png')), plot = gg_pca_fll, device = 'png', width = 6, height = 6, units = 'in', dpi = 350)
        ggplot2::ggsave(filename = file.path(outd,paste0(outdir,'/Figure_annexes_',src,'_pca_final.png')), plot = gg_pca_res, device = 'png', width = 6, height = 6, units = 'in', dpi = 350)
        res <- data.frame(index = as.numeric(index_cal(pca_res)))
        names(res) <- src
        return(res)
      }) |>
      dplyr::bind_cols()
    metrics_pca_idx$station <- metrics_lst[[1]] |> dplyr::pull(station)
    metrics_pca_idx$longitude <- metrics_lst[[1]] |> dplyr::pull(longitude)
    metrics_pca_idx$latitude  <- metrics_lst[[1]] |> dplyr::pull(latitude)
    metrics_pca_idx$altitude  <- metrics_lst[[1]] |> dplyr::pull(altitude)
    metrics_pca_idx[metrics_pca_idx == 50] <- NA
    return(metrics_pca_idx)
  }
}

ideam_fdrrz_evaluation_mth <- get_composite_index(metrics = ideam_general_mth, supplementary = fdrrz_general_mth, outdir = 'results/fedearroz/figures')
ideam_fdrrz_evaluation_mth <- ideam_fdrrz_evaluation_mth[,c('station',base::setdiff(names(ideam_fdrrz_evaluation_mth),'station'))]
ideam_fdrrz_evaluation_mth$best_source <- apply(X = ideam_fdrrz_evaluation_mth[,2:5], MARGIN = 1, FUN = function(x){which(x == max(x, na.rm = T))})
ideam_fdrrz_evaluation_mth$best_svalue <- apply(X = ideam_fdrrz_evaluation_mth[,2:5], MARGIN = 1, FUN = max, na.rm = T)
ideam_fdrrz_evaluation_mth$best_source[ideam_fdrrz_evaluation_mth$best_source == 1] <- 'AgERA5'
ideam_fdrrz_evaluation_mth$best_source[ideam_fdrrz_evaluation_mth$best_source == 2] <- 'CHIRPS'
ideam_fdrrz_evaluation_mth$best_source[ideam_fdrrz_evaluation_mth$best_source == 3] <- 'IMERG'
ideam_fdrrz_evaluation_mth$best_source[ideam_fdrrz_evaluation_mth$best_source == 4] <- 'MSWEP'
ideam_fdrrz_evaluation_mth$best_source <- factor(ideam_fdrrz_evaluation_mth$best_source)

ideam_cencf_evaluation_mth <- get_composite_index(metrics = ideam_general_mth, supplementary = tidyr::drop_na(cencf_general_mth), outdir = 'results/cenicafe/figures')
ideam_cencf_evaluation_mth <- ideam_cencf_evaluation_mth[,c('station',base::setdiff(names(ideam_cencf_evaluation_mth),'station'))]
ideam_cencf_evaluation_mth$best_source <- apply(X = ideam_cencf_evaluation_mth[,2:5], MARGIN = 1, FUN = function(x){which(x == max(x, na.rm = T))})
ideam_cencf_evaluation_mth$best_svalue <- apply(X = ideam_cencf_evaluation_mth[,2:5], MARGIN = 1, max, na.rm = T)
ideam_cencf_evaluation_mth$best_source[ideam_cencf_evaluation_mth$best_source == 1] <- 'AgERA5'
ideam_cencf_evaluation_mth$best_source[ideam_cencf_evaluation_mth$best_source == 2] <- 'CHIRPS'
ideam_cencf_evaluation_mth$best_source[ideam_cencf_evaluation_mth$best_source == 3] <- 'IMERG'
ideam_cencf_evaluation_mth$best_source[ideam_cencf_evaluation_mth$best_source == 4] <- 'MSWEP'
ideam_cencf_evaluation_mth$best_source <- factor(ideam_cencf_evaluation_mth$best_source)

ideam_ideamnew_evaluation_mth <- get_composite_index(metrics = ideam_general_mth, supplementary = ideamnew_general_mth, outdir = 'results/ideamnew/figures')
ideam_ideamnew_evaluation_mth <- ideam_ideamnew_evaluation_mth[,c('station',base::setdiff(names(ideam_ideamnew_evaluation_mth),'station'))]
ideam_ideamnew_evaluation_mth$best_source <- apply(X = ideam_ideamnew_evaluation_mth[,2:5], MARGIN = 1, FUN = function(x){which(x == max(x, na.rm = T))})
ideam_ideamnew_evaluation_mth$best_svalue <- apply(X = ideam_ideamnew_evaluation_mth[,2:5], MARGIN = 1, max, na.rm = T)
ideam_ideamnew_evaluation_mth$best_source[ideam_ideamnew_evaluation_mth$best_source == 1] <- 'AgERA5'
ideam_ideamnew_evaluation_mth$best_source[ideam_ideamnew_evaluation_mth$best_source == 2] <- 'CHIRPS'
ideam_ideamnew_evaluation_mth$best_source[ideam_ideamnew_evaluation_mth$best_source == 3] <- 'IMERG'
ideam_ideamnew_evaluation_mth$best_source[ideam_ideamnew_evaluation_mth$best_source == 4] <- 'MSWEP'
ideam_ideamnew_evaluation_mth$best_source <- factor(ideam_ideamnew_evaluation_mth$best_source)

col_shp2 <- geodata::gadm(country = 'COL', level = 2, path = tempdir(), version = 'latest')
col_shp2_sf <- sf::st_as_sf(col_shp2)

ideam_fdrrz_evaluation_mth <- cbind(ideam_fdrrz_evaluation_mth,
                                    terra::extract(x = col_shp2, y = ideam_fdrrz_evaluation_mth[,c('longitude','latitude')])[,c('NAME_1','NAME_2')])
ideam_fdrrz_evaluation_mth <- ideam_fdrrz_evaluation_mth[,c('station','longitude','latitude','NAME_1','NAME_2','altitude','AgERA5','CHIRPS','IMERG','MSWEP','best_source','best_svalue')]
names(ideam_fdrrz_evaluation_mth)[4:5] <- c('departamento','municipio')
write.csv(x = ideam_fdrrz_evaluation_mth, file.path(outd,'results/fedearroz/fedearroz_calidad_estaciones.csv'), row.names = F, fileEncoding = 'latin1')

ideam_cencf_evaluation_mth <- cbind(ideam_cencf_evaluation_mth,
                                    terra::extract(x = col_shp2, y = ideam_cencf_evaluation_mth[,c('longitude','latitude')])[,c('NAME_1','NAME_2')])
ideam_cencf_evaluation_mth <- ideam_cencf_evaluation_mth[,c('station','longitude','latitude','NAME_1','NAME_2','altitude','AgERA5','CHIRPS','IMERG','MSWEP','best_source','best_svalue')]
names(ideam_cencf_evaluation_mth)[4:5] <- c('departamento','municipio')
write.csv(x = ideam_cencf_evaluation_mth, file.path(outd,'results/cenicafe/cenicafe_calidad_estaciones.csv'), row.names = F, fileEncoding = 'latin1')

ideam_ideamnew_evaluation_mth <- cbind(ideam_ideamnew_evaluation_mth,
                                       terra::extract(x = col_shp2, y = ideam_ideamnew_evaluation_mth[,c('longitude','latitude')])[,c('NAME_1','NAME_2')])
ideam_ideamnew_evaluation_mth <- ideam_ideamnew_evaluation_mth[,c('station','longitude','latitude','NAME_1','NAME_2','altitude','AgERA5','CHIRPS','IMERG','MSWEP','best_source','best_svalue')]
names(ideam_ideamnew_evaluation_mth)[4:5] <- c('departamento','municipio')
write.csv(x = ideam_ideamnew_evaluation_mth, file.path(outd,'results/ideam/ideam_full_calidad_estaciones.csv'), row.names = F, fileEncoding = 'latin1')









# List of departments in which Rice and Plantain are harvested
dpts <- c('La Guajira','Magdalena','Cesar','Sucre','Córdoba','Antioquia',
'Norte de Santander','Casanare','Meta','Tolima','Huila','Valle del Cauca')
departments <- factor(x = dpts, levels = dpts)

col_shp <- geodata::gadm(country = 'COL', level = 1, path = tempdir(), version = 'latest')
col_shp_sf <- sf::st_as_sf(col_shp)
col_flt <- col_shp[col_shp$NAME_1 %in% as.character(departments),]
col_flt_sf <- sf::st_as_sf(col_flt)

merged_unq <- rbind(ideam_unq[,c('longitude','latitude','station')], fdrrz_unq |> dplyr::rename(longitude = longitud, latitude = latitud) |> dplyr::select(longitude, latitude, station))
stations_in_departments_dfm <- terra::intersect(x = col_flt, y = terra::vect(merged_unq, c('longitude','latitude'), crs = 'EPSG:4326')) |> base::as.data.frame()
stations_in_departments_dfm <- stations_in_departments_dfm[,c('NAME_1','station')]
names(stations_in_departments_dfm)[1] <- 'department'
stations_in_departments <- terra::intersect(x = col_flt, y = terra::vect(merged_unq, c('longitude','latitude'), crs = 'EPSG:4326')) |> base::as.data.frame() |> dplyr::pull(station)

# ideam_general_mth_dpts <- ideam_general_mth |> dplyr::filter(station %in% stations_in_departments)
# ideam_general_qrt_dpts <- ideam_general_qrt |> dplyr::filter(station %in% stations_in_departments)
# ideam_per_mth_dpts <- ideam_per_mth |> dplyr::filter(station %in% stations_in_departments)
# ideam_per_qrt_dpts <- ideam_per_qrt |> dplyr::filter(station %in% stations_in_departments)
# arrow::write_parquet(x = ideam_general_mth_dpts, sink = file.path(outd,'ideam_evaluation_year-month_general_rice_dpts.parquet'), version = 'latest')
# arrow::write_parquet(x = ideam_general_qrt_dpts, sink = file.path(outd,'ideam_evaluation_year-quarter_general_rice_dpts.parquet'), version = 'latest')
# arrow::write_parquet(x = ideam_per_mth_dpts, sink = file.path(outd,'ideam_evaluation_per-month_rice_dpts.parquet'), version = 'latest')
# arrow::write_parquet(x = ideam_per_qrt_dpts, sink = file.path(outd,'ideam_evaluation_per-quarter_rice_dpts.parquet'), version = 'latest')

# Figure A1. General climatotology
gg_a1 <- rbind(ideam_mrg, fdrrz_mrg) |>
  dplyr::filter(station %in% stations_in_departments) |>
  tidyr::drop_na() |>
  dplyr::select(fecha,station,valor_observado) |>
  dplyr::mutate(month = lubridate::month(fecha)) |>
  dplyr::group_by(month) |>
  dplyr::summarise(valor_observado = median(valor_observado)) |>
  ggplot2::ggplot(aes(x = factor(month), y = valor_observado)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('Mes') +
  ggplot2::ylab('Precipitación (mm/mes)') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_A1.png'), plot = gg_a1, device = 'png', width = 5, height = 4, units = 'in', dpi = 350)

# Figure A2. Climatotology per department
gg_a2 <- rbind(ideam_mrg, fdrrz_mrg) |>
  dplyr::filter(station %in% stations_in_departments) |>
  tidyr::drop_na() |>
  dplyr::left_join(y = stations_in_departments_dfm, by = 'station') |>
  dplyr::select(fecha,station,valor_observado,department) |>
  dplyr::mutate(month = lubridate::month(fecha)) |>
  dplyr::group_by(month,department) |>
  dplyr::summarise(valor_observado = median(valor_observado)) |>
  dplyr::mutate(department = factor(department, levels = dpts)) |>
  ggplot2::ggplot(aes(x = factor(month), y = valor_observado)) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::xlab('Mes') +
  ggplot2::ylab('Precipitación (mm/mes)') +
  ggplot2::facet_wrap(~department) +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_A2.png'), plot = gg_a2, device = 'png', width = 8, height = 5, units = 'in', dpi = 350)

# -------------------------------------- #
# Temporal
# -------------------------------------- #

gg_tp1 <- rbind(ideam_per_mth, fdrrz_per_mth) |>
  dplyr::filter(station %in% stations_in_departments) |>
  dplyr::select(rsq:source,filter) |>
  tidyr::pivot_longer(cols = rsq:maape, names_to = 'metric', values_to = 'value') |>
  dplyr::group_by(source, filter, metric) |>
  dplyr::summarise(value = median(value, na.rm = T)) |>
  dplyr::ungroup() |>
  # dplyr::filter(metric == 'rmse') |>
  ggplot2::ggplot(aes(x = as.factor(as.numeric(filter)), y = value, group = source, color = source)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::facet_wrap(~metric, scales = 'free') +
  ggplot2::xlab('Mes') +
  ggplot2::ylab('') +
  ggplot2::labs(color = '') +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = 'bottom')
ggplot2::ggsave(filename = file.path(outd,'Figure_TP1.png'), plot = gg_tp1, device = 'png', width = 12.5, height = 6, units = 'in', dpi = 350)

gg_tp2 <- rbind(ideam_per_mth, fdrrz_per_mth) |>
  dplyr::filter(station %in% stations_in_departments) |>
  dplyr::left_join(y = stations_in_departments_dfm, by = 'station') |>
  dplyr::select(rsq:source,filter,department) |>
  tidyr::pivot_longer(cols = rsq:maape, names_to = 'metric', values_to = 'value') |>
  dplyr::group_by(source, filter, metric, department) |>
  dplyr::summarise(value = median(value, na.rm = T)) |>
  dplyr::ungroup() |>
  dplyr::mutate(department = factor(department, levels = dpts)) |>
  dplyr::filter(metric == 'Spearman') |>
  ggplot2::ggplot(aes(x = as.factor(as.numeric(filter)), y = value, group = source, color = source)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::facet_wrap(~department, scales = 'fixed') +
  ggplot2::xlab('Mes') +
  ggplot2::ylab('') +
  ggplot2::labs(color = '') +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = 'bottom')
ggplot2::ggsave(filename = file.path(outd,'Figure_TP2.png'), plot = gg_tp2, device = 'png', width = 12.5, height = 6, units = 'in', dpi = 350)

# -------------------------------------- #
# Spatial
# -------------------------------------- #
gg_sp0 <- ideam_general_mth |>
  dplyr::filter(station %in% stations_in_departments) |>
  ggplot2::ggplot(aes(x = reorder(source, -Spearman, FUN = median), y = Spearman, color = source)) +
  ggplot2::stat_summary(fun = median, geom = 'point', size = 5) +
  ggplot2::geom_jitter(position = position_jitter(), alpha = 0.1) +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::theme_minimal() +
  ggplot2::coord_flip() +
  ggplot2::xlab('') +
  ggplot2::ylab('Coeficiente de correlación de Spearman') +
  ggplot2::theme(legend.position = 'none')
ggplot2::ggsave(filename = file.path(outd,'Figure_SP0.png'), plot = gg_sp0, device = 'png', width = 6, height = 4, units = 'in', dpi = 350)

gg_sp1 <- ideam_evaluation_mth |>
  dplyr::filter(station %in% stations_in_departments) |>
  tidyr::pivot_longer(cols = AgERA5:MSWEP, names_to = 'source', values_to = 'quality') |>
  ggplot2::ggplot(aes(x = reorder(source, -quality, FUN = median), y = quality, color = source)) +
  ggplot2::stat_summary(fun = median, geom = 'point', size = 5) +
  ggplot2::geom_jitter(position = position_jitter(), alpha = 0.1) +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::ylim(0,100) +
  ggplot2::theme_minimal() +
  ggplot2::coord_flip() +
  ggplot2::xlab('') +
  ggplot2::ylab('Índice de calidad') +
  ggplot2::theme(legend.position = 'none')
ggplot2::ggsave(filename = file.path(outd,'Figure_SP1.png'), plot = gg_sp1, device = 'png', width = 6, height = 4, units = 'in', dpi = 350)

aux <- ideam_evaluation_mth[ideam_evaluation_mth$station %in% stations_in_departments,] |>
  dplyr::select(AgERA5:latitude) |>
  tidyr::pivot_longer(cols = AgERA5:MSWEP, names_to = 'source', values_to = 'quality')
gg_sp2 <- ggplot2::ggplot(data = col_flt_sf) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf(data = col_shp_sf, fill = NA) +
  ggplot2::geom_point(data = aux, aes(x = longitude, y = latitude, color = quality), size = 0.4) +
  ggplot2::scale_colour_gradientn(colours = c('#d7191c','#fdae61','#ffffbf','#a6d96a','#1a9641'), limits = c(0,100)) + # c('#D5AC0F','#F8DC72','#86C2E9','#2980B9','#1B5276')
  ggplot2::xlim(-77.6, -69.7) +
  ggplot2::ylim(1.4, 12.5) +
  ggplot2::xlab('Longitud') +
  ggplot2::ylab('Latitud') +
  ggplot2::labs(color = 'Índice de calidad') +
  ggplot2::facet_wrap(~source, scales = 'fixed') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_SP2.png'), plot = gg_sp2, device = 'png', width = 6, height = 6.5, units = 'in', dpi = 350)

gg_sp3 <- ggplot2::ggplot(data = col_flt_sf) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf(data = col_shp_sf, fill = NA) +
  ggplot2::geom_point(data = ideam_evaluation_mth[ideam_evaluation_mth$station %in% stations_in_departments,], aes(x = longitude, y = latitude, color = best_source), size = 0.5) +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::xlim(-77.6, -69.7) +
  ggplot2::ylim(1.4, 12.5) +
  ggplot2::xlab('Longitud') +
  ggplot2::ylab('Latitud') +
  ggplot2::labs(color = 'Mejor fuente\npor estación') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_SP3.png'), plot = gg_sp3, device = 'png', width = 4, height = 3.5, units = 'in', dpi = 350)

gg_sp4 <- ideam_evaluation_mth |>
  dplyr::filter(station %in% stations_in_departments) |>
  dplyr::left_join(stations_in_departments_dfm, by = 'station') |>
  dplyr::select(best_source, department) |>
  table() |>
  base::as.data.frame() |>
  dplyr::mutate(department = factor(department, levels = dpts)) |>
  ggplot2::ggplot(aes(x = department, y = Freq, fill = best_source)) + 
  ggplot2::geom_bar(position = 'fill', stat = 'identity') +
  ggplot2::scale_fill_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::coord_flip() +
  ggplot2::xlab('') +
  ggplot2::ylab('Porcentaje de estaciones (%)') +
  ggplot2::labs(fill = 'Fuente') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_SP4.png'), plot = gg_sp4, device = 'png', width = 6, height = 4, units = 'in', dpi = 350)

# High quality stations
ideam_evaluation_mth |>
  dplyr::filter(station %in% stations_in_departments) |>
  dplyr::left_join(stations_in_departments_dfm, by = 'station') |>
  dplyr::filter(best_svalue > 50) |>
  dplyr::select(best_source, department) |>
  table() |>
  base::as.data.frame() |>
  dplyr::mutate(department = factor(department, levels = dpts)) |>
  ggplot2::ggplot(aes(x = department, y = Freq, fill = best_source)) + 
  ggplot2::geom_bar(stat="identity", position=position_dodge()) +
  ggplot2::scale_fill_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::facet_wrap(~best_source) +
  ggplot2::coord_flip() +
  ggplot2::xlab('') +
  ggplot2::ylab('Número de estaciones') +
  ggplot2::labs(fill = 'Fuente') +
  ggplot2::theme_bw()

aux <- ideam_evaluation_mth[ideam_evaluation_mth$station %in% stations_in_departments,] |>
  dplyr::select(AgERA5:altitude) |>
  tidyr::pivot_longer(cols = AgERA5:MSWEP, names_to = 'source', values_to = 'quality')
gg_sp5 <- aux |>
  ggplot2::ggplot(aes(x = altitude, y = quality, color = source)) +
  ggplot2::geom_point(alpha = 0.1) +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::geom_smooth(aes(color = factor(source)), method = "lm", formula = y ~ x, se = F) +
  #ggplot2::ylim(0, 100) +
  ggplot2::xlab('Altitud (MSNM)') +
  ggplot2::ylab('Índice de calidad') +
    ggplot2::labs(color = 'Fuente') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_SP5.png'), plot = gg_sp5, device = 'png', width = 6, height = 4, units = 'in', dpi = 350)

ideam_evaluation_mth$slope <- terra::extract(x = slope, y = ideam_evaluation_mth[,c('longitude','latitude')]) |> dplyr::pull(slope)
aux <- ideam_evaluation_mth[ideam_evaluation_mth$station %in% stations_in_departments,] |>
  dplyr::select(AgERA5:MSWEP,slope) |>
  tidyr::pivot_longer(cols = AgERA5:MSWEP, names_to = 'source', values_to = 'quality')
gg_sp6 <- aux |>
  ggplot2::ggplot(aes(x = slope, y = quality, color = source)) +
  ggplot2::geom_point(alpha = 0.1) +
  ggplot2::scale_color_manual(values = c('#2596be','#C79F40','#7FC343','#90BBC8')) +
  ggplot2::geom_smooth(aes(color = factor(source)), method = "lm", formula = y ~ x, se = F) +
  #ggplot2::ylim(0, 100) +
  ggplot2::xlab('Pendiente topográfica') +
  ggplot2::ylab('Índice de calidad') +
    ggplot2::labs(color = 'Fuente') +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(outd,'Figure_SP6.png'), plot = gg_sp6, device = 'png', width = 6, height = 4, units = 'in', dpi = 350)

# Disclaimer table
tbl2shr <- ideam_evaluation_mth[ideam_evaluation_mth$station %in% stations_in_departments,]
rownames(tbl2shr) <- 1:nrow(tbl2shr)

names(tbl2shr) <- c('AgERA5','CHIRPS','IMERG','MSWEP','Estacion','Longitud','Latitud','Altitud','Mejor_fuente','Mejor_valor')

col_shp2 <- geodata::gadm(country = 'COL', level = 2, path = tempdir(), version = 'latest')
col_shp2_sf <- sf::st_as_sf(col_shp2)

stations_in_municipalities_dfm <- terra::intersect(x = col_shp2, y = terra::vect(tbl2shr, c('Longitud','Latitud'), crs = 'EPSG:4326')) |> base::as.data.frame()
tbl2shr <- cbind(tbl2shr[,c('Longitud','Latitud')], stations_in_municipalities_dfm[,c('Estacion','NAME_1','NAME_2','Altitud','AgERA5','CHIRPS','IMERG','MSWEP','Mejor_fuente','Mejor_valor')])
tbl2shr <- tbl2shr[,c('Estacion','Longitud','Latitud','NAME_1','NAME_2','Altitud','AgERA5','CHIRPS','IMERG','MSWEP','Mejor_fuente','Mejor_valor')]
names(tbl2shr)[4:5] <- c('Departamento','Municipio')
write.csv(x = tbl2shr, file.path(outd,'indice_calidad_estacion.csv'), row.names = F, fileEncoding = 'latin1')

## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##

## Recommendation system ----
#(CHIRPS, AgERA5, MSWEP, IMERG quality) ~ f(station, orography, cloud cover, etc)

## Potential predictors of the quality ----
# Rice harvested areas from CROPGRIDS
cropgrids <- terra::rast('D:/Data/Maps/cropgrids/NC_maps/CROPGRIDSv1.06_rice.nc')[['harvarea']]
names(cropgrids) <- 'cropgrids'
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = cropgrids, y = ideam_mtr[,c('longitude','latitude')]))); rm(cropgrids)
ideam_mtr$ID <- NULL
plot(ideam_mtr$cropgrids, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$cropgrids, ideam_mtr$quality, method = 'spearman')
abline(lm(quality~cropgrids, data = ideam_mtr), col = 'red')
# Rice harvested areas from MapSPAM
mapspam <- terra::rast('D:/Data/Maps/spam2010/spam2010V2r0_global_H_RICE_A.tif')
names(mapspam) <- 'mapspam'
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = mapspam, y = ideam_mtr[,c('longitude','latitude')]))); rm(mapspam)
ideam_mtr$ID <- NULL
plot(ideam_mtr$mapspam, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$mapspam, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~mapspam, data = ideam_mtr), col = 'red')

plot(ideam_mtr$mapspam, ideam_mtr$cropgrids, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$mapspam, ideam_mtr$cropgrids, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(cropgrids~mapspam, data = ideam_mtr), col = 'red')
abline(0, 1, col = 'blue')

# Cloud cover - annual mean
cloud_cover <- terra::rast(file.path(root,'MODCF_meanannual.tif'))
cloud_cover <- cloud_cover * 0.01
names(cloud_cover) <- 'cloud_cover_mean'
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = cloud_cover, y = ideam_mtr[,c('longitude','latitude')]))); rm(cloud_cover)
ideam_mtr$ID <- NULL
plot(ideam_mtr$cloud_cover_mean, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$cloud_cover_mean, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~cloud_cover_mean, data = ideam_mtr), col = 'red')
# Cloud cover - inter annual SD
cloud_cover <- terra::rast(file.path(root,'MODCF_interannualSD.tif'))
names(cloud_cover) <- 'cloud_cover_interannualsd'
cloud_cover <- cloud_cover * 0.01
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = cloud_cover, y = ideam_mtr[,c('longitude','latitude')]))); rm(cloud_cover)
ideam_mtr$ID <- NULL
plot(ideam_mtr$cloud_cover_interannualsd, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$cloud_cover_interannualsd, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~cloud_cover_interannualsd, data = ideam_mtr), col = 'red')

# Accessibility to cities
access <- terra::rast('D:/Data/Maps/accessibility/accessibility_to_cities_2015_v1.0.tif')
names(access) <- 'accessibility'
access[access == -9999] <- NA
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = access, y = ideam_mtr[,c('longitude','latitude')]))); rm(access)
ideam_mtr$ID <- NULL
plot(ideam_mtr$accessibility, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$accessibility, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~accessibility, data = ideam_mtr), col = 'red')

# Population density
pop <- terra::rast('D:/Data/Maps/population/gpw_v4_population_density_rev11_2020_2pt5_min.tif')
names(pop) <- 'population_density'
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = pop, y = ideam_mtr[,c('longitude','latitude')]))); rm(pop)
ideam_mtr$ID <- NULL
plot(ideam_mtr$population_density, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$population_density, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~population_density, data = ideam_mtr), col = 'red')

# Hillshade
pacman::p_load(geodata)
dem <- geodata::elevation_30s(country = 'COL', path = tempdir())
slope <- terra::terrain(dem, v = 'slope', unit = 'radians')
aspect <- terra::terrain(dem, v = 'aspect', unit = 'radians')
hillshade <- terra::shade(slope, aspect, angle = 45, direction = 315)
rm(dem, slope, aspect)
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = hillshade, y = ideam_mtr[,c('longitude','latitude')]))); rm(hillshade)
ideam_mtr$ID <- NULL
plot(ideam_mtr$hillshade, ideam_mtr$quality, col = alpha('black', 0.4), pch = 20, cex = 0.6)
cor(ideam_mtr$hillshade, ideam_mtr$quality, method = 'spearman', use = 'pairwise.complete.obs')
abline(lm(quality~hillshade, data = ideam_mtr), col = 'red')

## Classifications ----
koppen <- terra::rast(file.path(root,'koppen_geiger_0p00833333.tif'))
names(koppen) <- 'koppen_geiger'
ideam_mtr <- cbind(ideam_mtr,base::as.data.frame(terra::extract(x = koppen, y = ideam_mtr[,c('longitude','latitude')]))); rm(koppen)
ideam_mtr$ID <- NULL
ideam_mtr$koppen_geiger <- as.character(ideam_mtr$koppen_geiger)

ideam_mtr |>
  dplyr::mutate(koppen_geiger = factor(koppen_geiger,
                                       levels = gtools::mixedsort(unique(ideam_mtr$koppen_geiger)),
                                       labels = c('Tropical, rainforest','Tropical, monsoon','Tropical, savannah',
                                                  'Arid, desert, hot','Arid, steppe, hot','Temperate, dry summer, warm summer',
                                                  'Temperate, dry winter, warm summer','Temperate, no dry season, warm summer','Polar, tundra'))) |>
  ggplot2::ggplot(aes(x = reorder(koppen_geiger, quality), y = quality, color = koppen_geiger)) +
  ggplot2::stat_summary(fun = median, geom = 'point', size = 8, show.legend = F) +
  ggplot2::geom_jitter(alpha = 0.3, show.legend = F) +
  ggplot2::theme_minimal() +
  ggplot2::coord_flip()
