## ------------------------------------------ ##
## CSICAP - Satellite products evaluation
## Extract satellite data
## By: Harold Achicanoy
## Alliance Bioversity-CIAT
## Nov. 2024
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)}
pacman::p_load(terra, arrow, googledrive, cloudml, data.table, tidyverse,
               lubridate, Metrics, FactoMineR, factoextra, scales, minerva,
               readxl)

## Define directories ----
root <- '//CATALOGUE/WFP_ClimateRiskPr1'
outd <- paste0(root,'/7.Results/CSICAP')
# root <- 'C:/Users/haachicanoy/Downloads'
# outd <- root

## Load observed data ----
### IDEAM data ----
if(file.exists(file.path(outd,'ideam_prec_monthly.parquet'))|file.exists(file.path(outd,'ideam_feat_stations.parquet'))){
  ideam_mnt <- arrow::read_parquet(file = file.path(outd,'ideam_prec_monthly.parquet')) |> base::as.data.frame()
  ideam_unq <- arrow::read_parquet(file = file.path(outd,'ideam_feat_stations.parquet')) |> base::as.data.frame()
} else {
  # Loading IDEAM data
  ideam <- arrow::read_parquet(file = file.path(outd,'df_final_georeferencias.parquet'))
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
  # # Plotting 10 random stations
  # ideam_mnt |>
  #   dplyr::filter(codigo_estacion %in% sample(x = ideam_unq$codigo_estacion, size = 10)) |>
  #   ggplot2::ggplot(aes(x = fecha, y = valor_observado, group = codigo_estacion)) +
  #   ggplot2::geom_line(alpha = 0.1, show.legend = F) +
  #   ggplot2::theme_minimal()
  # Save monthly data
  arrow::write_parquet(x = ideam_mnt, sink = file.path(outd,'ideam_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = ideam_unq, sink = file.path(outd,'ideam_feat_stations.parquet'), version = 'latest')
  rm(ideam, ideamDT); gc(T)
}
ideam_unq <- ideam_unq[,c(c('longitude','latitude'),base::setdiff(names(ideam_unq),c('longitude','latitude')))]
names(ideam_unq)[3] <- 'station'
names(ideam_mnt)[2] <- 'station'

### Fedearroz data ----
if(file.exists(file.path(outd,'fedearroz_prec_monthly.parquet'))|file.exists(file.path(outd,'fedearroz_feat_stations.parquet'))){
  fdrrz_mnt <- arrow::read_parquet(file = file.path(outd,'fedearroz_prec_monthly.parquet')) |> base::as.data.frame()
  fdrrz_unq <- arrow::read_parquet(file = file.path(outd,'fedearroz_feat_stations.parquet')) |> base::as.data.frame()
} else {
  fdrrz_mnt <- readxl::read_excel(path = file.path(outd,'prec_monthly_obs_sat.xlsx'), sheet = 1) |> base::as.data.frame()
  fdrrz_mnt <- fdrrz_mnt[,-1]
  fdrrz_mnt <- fdrrz_mnt[fdrrz_mnt$fuente == 'fedearroz',]; rownames(fdrrz_mnt) <- 1:nrow(fdrrz_mnt)
  fdrrz_mnt$fuente <- NULL
  fdrrz_unq <- fdrrz_mnt[,c('longitud','latitud','station','dpto','mun')] |> unique()
  fdrrz_mnt <- fdrrz_mnt[,c('month_year','station','prec_month')]
  names(fdrrz_mnt) <- c('fecha','station','valor_observado')
  fdrrz_mnt$fecha <- as.Date(paste0(fdrrz_mnt$fecha,'-01'))
  arrow::write_parquet(x = fdrrz_mnt, sink = file.path(outd,'fedearroz_prec_monthly.parquet'), version = 'latest')
  arrow::write_parquet(x = fdrrz_unq, sink = file.path(outd,'fedearroz_feat_stations.parquet'), version = 'latest')
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
if(!file.exists(file.path(outd,'ideam_prec_monthly_merged.parquet'))|
   !file.exists(file.path(outd,'fedearroz_prec_monthly_merged.parquet'))){
  ideam_mrg <- extract_satellite_data(crds = ideam_unq, Data = ideam_mnt)
  ideam_mrg$IMERG <- ideam_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  fdrrz_mrg <- extract_satellite_data(crds = fdrrz_unq, Data = fdrrz_mnt)
  fdrrz_mrg$IMERG <- fdrrz_mrg$IMERG * 730.5 # To get mm/month instead of mm/hr
  arrow::write_parquet(x = ideam_mrg, sink = file.path(outd,'ideam_prec_monthly_merged.parquet'), version = 'latest')
  arrow::write_parquet(x = fdrrz_mrg, sink = file.path(outd,'fedearroz_prec_monthly_merged.parquet'), version = 'latest')
} else {
  ideam_mrg <- arrow::read_parquet(file.path(outd,'ideam_prec_monthly_merged.parquet'))
  fdrrz_mrg <- arrow::read_parquet(file.path(outd,'fedearroz_prec_monthly_merged.parquet'))
}
