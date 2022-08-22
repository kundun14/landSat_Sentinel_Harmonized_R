# spatial co-registration Sentinel MSI + OLI (LS 8 & 9 + S2 1A + 2A)
#stackable" for time series analysis


# paquetes
packages <- c('leaflet','rgdal','raster','jsonlite','sp','httr','rasterVis','ggplot2','magrittr', 'RColorBrewer','xml2','dygraphs','xts','lubridate','DT','rmarkdown', 'rprojroot')
# paquetes no instalados aun
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
# installar paquetes no isntalados aun
if(length(new.packages)) install.packages(new.packages, repos='http://cran.rstudio.com/') else print('All required packages are installed.')
# cargar paquetes
invisible(lapply(packages, library, character.only = TRUE)) # apply la funcion library a vector que contine los paquetes
#directorio
wd <- rprojroot::find_rstudio_root_file()
#creams una carpeta para los outputs (CVS de estadisiticos y rasters)
outDir <- file.path(wd, "Data", "R_Output", fsep="/")
suppressWarnings(dir.create(outDir)) 

## parametros de busqueda NASA's CMR-STAC API.
search_URL <- "https://cmr.earthdata.nasa.gov/stac/LPCLOUD/search"
HLS_col <- list("HLSS30.v2.0", "HLSL30.v2.0") # colecciones landsat 8 9 + sentinel 2

#parcela
parcela <-  readOGR("./Data/parcela.shp")
#mapa
leaflet() %>% 
  addPolygons(data = parcela, fill = FALSE) %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addMiniMap(zoomLevelFixed = 5)
#parametros de busqueda
#bbox, crs +proj=longlat +datum=WGS84
roi <- extent(parcela)
parcela_bbox <- paste(roi[1], roi[3], roi[2], roi[4], sep = ',')
#"-71.9089581485665,-13.5953803532065,-71.8676564870626,-13.5657056500853"
#rango de tiempo
parcela_datetime <- '2021-05-01T00:00:00Z/2021-05-30T23:59:59Z'  # 2 meses data # YYYY-MM-DDTHH:MM:SSZ/YYYY-MM-DDTHH:MM:SSZ
#parametros
cropland_search_body <- list(limit=100,
                             datetime=parcela_datetime,
                             bbox= parcela_bbox,
                             collections= HLS_col)

cropland_search_req <- httr::POST(search_URL, body = cropland_search_body, encode = "json") %>% 
  httr::content(as = "text") %>% 
  fromJSON()
cat('There are',cropland_search_req$numberMatched, 'features matched our request.')
#There are 28 features matched our request.
names(cropland_search_req)
cropland_search_req$features # todos los tiffs estan aqui
# visualizar el primer feature

search_features <- cropland_search_req$features
Feature1 <- search_features[1,]
Feature1 %>% 
  jsonlite::toJSON(auto_unbox = TRUE, pretty = TRUE) 

#imagen del primer feature

browse_image_url <- Feature1$assets$browse$href
browse_req <- GET(browse_image_url) %>% 
  httr::content()  

plot(0:1, 0:1, type = "n",ann=FALSE,axes=FALSE)
rasterImage(browse_req,0, 0, 1,1) # esta imagen tiene muchas nubes

# 'eo:cloud_cover' property

for (i in row.names(search_features)){
  cc <- search_features[i,]$properties$`eo:cloud_cover`
  d <- search_features[i,]$properties$datetime
  cat('Cloud Coverage of Feature ',i, 'Captured on ', d , 'Is: ' , cc , 'Percent' ,'\n')
}

#clip (subset parcela) + selecccionar bandas 

#set gdal config para acceder a assets en la nube

rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_UNSAFESSL", value = "YES")
rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_COOKIEFILE", value = ".rcookies")
rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_COOKIEJAR", value = ".rcookies")
rgdal::setCPLConfigOption(ConfigOption = "GDAL_DISABLE_READDIR_ON_OPEN", value = "YES")
rgdal::setCPLConfigOption(ConfigOption = "CPL_VSIL_CURL_ALLOWED_EXTENSIONS", value = "TIF")

# nombre de las bandas



# All Landsat-8 OLI and Sentinel-2 MSI reflective spectral bands nomenclatures are retained in
# the HLS products

#landsat
# NIr B8a
# red B04
# quality Fmask
# #sentinel
# NIr B05
# red B04
# quality Fmask

#dataframe de  las imagenes y sus propiedades

granule_list <- list()
n <- 1
for (item in row.names(search_features)){                       # Get the NIR, Red, and Quality band layer names
  if (search_features[item,]$collection == 'HLSS30.v2.0'){
    ndvi_bands <- c('B8A','B04','Fmask')
  }
  else{
    ndvi_bands <- c('B05','B04','Fmask')
  }
  for(b in ndvi_bands){
    f <- search_features[item,]
    b_assets <- f$assets[[b]]$href
    
    df <- data.frame(Collection = f$collection,                    # Make a data frame including links and other info
                     Granule_ID = f$id,
                     Cloud_Cover = f$properties$`eo:cloud_cover`,
                     band = b,
                     Asset_Link = b_assets, stringsAsFactors=FALSE)
    granule_list[[n]] <- df
    n <- n + 1
  }
}

# Create a searchable datatable
search_df <- do.call(rbind, granule_list)
DT::datatable(search_df)

# imagenes con cob nubes > 30%
search_df <- search_df[search_df$Cloud_Cover < 30, ]

#subset and stacking

# se necesita ncredenciales de earthdata NASA
# earthdata_netrc_setup.R script.


cat('The netrc file can be found in:', Sys.getenv('HOME'))
#The netrc file can be found in: C:\Users\LENOVO

# reprojeccion del area de estudio hacia el CRS del producto HLS

crs(parcela) # lon lat "EPSG",4326]

coordinate_reference <- raster(paste0(search_df$Asset_Link[1])) # crs de los assets
parcela_utm <- spTransform(parcela, crs(coordinate_reference)) # Transfer CRS
crs(parcela_utm)

# lista de rasters para cada banda  (i.e., Red, NIR, and Fmask). 
# cada banda es leida en la memoria, luego cropped en el area de estudio 
# los raster cropped seran stacekados 
# en esos rasters se utilizara una mascara para nubes y se calculara el NDVI

red_stack <- nir_stack <- fmask_stack <- list() # listas para los stacks de cada banda
# Add progress bar
pb = txtProgressBar(min = 0, max = length(search_df$band), initial = 0, style = 3)
l <- m <- n  <- 0
for (row in seq(length(search_df$band))){
  setTxtProgressBar(pb,row)
  if (search_df$band[row] == 'B04'){
    l = l+1
    red <- raster(search_df$Asset_Link[row])
    red_crop <- raster::mask (raster::crop(red, extent(parcela_utm)), parcela_utm)
    red_stack[[l]] <- red_crop
    rm (red, red_crop)
  }else if (search_df$band[row] == 'Fmask'){
    m = m+1
    fmask <- raster(search_df$Asset_Link[row])
    fmask_crop <- raster::mask (raster::crop(fmask, extent(parcela_utm)), parcela_utm)
    fmask_stack[[m]] <-  fmask_crop
    rm(fmask, fmask_crop)
  }else{
    n = n+1
    nir <- raster(search_df$Asset_Link[row])
    nir_crop <- raster::mask (raster::crop(nir, extent(parcela_utm)), parcela_utm)
    nir_stack [[n]] <- nir_crop
    rm(nir, nir_crop)
  }
}
close(pb)

# los stacks cropeados estan en la memoria?
inMemory(nir_stack[[1]]) # true

# calculo de NDVI 

calculate_NDVI <- function(nir, red) {
  ndvi <- (nir - red)/(nir + red)
  return(ndvi)
}

ndvi_stack <- list()
for (i in 1:length(nir_stack)){ 
  # Calcular NDVI 
  ndvi_stack[[i]] <- calculate_NDVI (nir_stack[[i]], red_stack[[i]])  
  # excluir Inf and -Inf del NDVI
  ndvi_stack[[i]][ndvi_stack[[i]] == Inf] <- NA
  ndvi_stack[[i]][ndvi_stack[[i]] == -Inf] <- NA
  layer_name <- names(nir_stack[[i]]) # fecha de cada raster layer
  date <- strsplit(layer_name, "[.]")[[1]][4]
  date <- substr(date[1], 1, 7) # fecha del producto ejemplo "2021129" == año 2021 dia 129
  date_time <- as.Date(as.integer(substr(date,5,7)), origin = paste0(substr(date,1,4), "-01-01"))
  
  if (strsplit(layer_name, ".T")[[1]][1] == 'HLS.S30'){           # fechas de los rasters NDVI
    names(ndvi_stack[[i]]) <- paste0('S',as.character(date_time)) # S de sentinel + fecha
    
  }else{
    
    names(ndvi_stack[[i]]) <- paste0('L',as.character(date_time)) # L de landsat + fecha
  }
}
ndvi_stacks <- stack(ndvi_stack)  # stack  NDVI

# hay imagenes para el mismo dia S2021.05.10.1, S2021.05.10.2 -------> por que???

# graficar
## paleta
pal <- colorNumeric(terrain.colors(n = 100), c(0,1) ,na.color = "transparent", reverse = TRUE)

leaflet() %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addRasterImage(ndvi_stacks[[3]], color = pal, opacity = 1) %>%
  addPolygons(data = parcela_utm, fill = FALSE) %>%
  leaflet::addLegend(pal = pal, values = c(0,1), title = "NDVI")


#Quality Filtering

# In HLS, both value of 0 and 64 in the Fmask layer indicate the pixel without cloud,
# cloud shadow, water, or snow/ice. A value of 0 also shows climatology aerosol level 
# and 64 shows low aerosol level.

mask_raster <- list()
ndvi_filtered <- list()

for (i in 1:length(fmask_stack)){
  mask_raster[[i]] <- fmask_stack[[i]]
  mask_raster[[i]][values(mask_raster[[i]])!= 0 || values(mask_raster[[i]])!= 64] <- NA
  ndvi_filtered[[i]] <- mask(ndvi_stacks[[i]], mask_raster[[i]], maskvalue=NA )
  names(ndvi_filtered[[i]]) <- names(ndvi_stacks[[i]])
}
ndvi_filtered_stacks <- stack(ndvi_filtered)


# series de tiempo por raster NDVI

base <- c("map<-leaflet()%>%
            addProviderTiles(providers$Esri.WorldImagery) %>%
            addMiniMap(zoomLevelFixed = 5) %>%")

# make a string including the addRasterImage function for every layer in the raster stack
X <- lapply(1:nlayers(ndvi_filtered_stacks), function(j){
  paste(paste("addRasterImage(ndvi_filtered_stacks[[",j,"]],
             colors = pal,
             opacity=1,
             group=names(ndvi_filtered_stacks[[",j,"]]))", sep=""),"%>% \n")})

X <- do.call(paste, X)

controls<-"addLayersControl(baseGroups=names(ndvi_stacks),
               options = layersControlOptions(collapsed=F), position = 'topleft')%>%"

legend <- "leaflet::addLegend(pal = pal, values = c(0,1), title = 'NDVI')"

final <- paste(base, X, controls, legend ,sep="\n")
eval(parse(text=final))
map


#graficos / plots
raster::boxplot(ndvi_filtered_stacks, col=c('olivedrab3'), 
                main='NDVI Time Series', ylab='NDVI', 
                names = names(ndvi_stacks), las=2)
#estadisticos
ndvi_mean <- cellStats(ndvi_filtered_stacks, stat='mean', na.rm=TRUE)
ndvi_max <- cellStats(ndvi_filtered_stacks, stat='max', na.rm=TRUE)
ndvi_min <- cellStats(ndvi_filtered_stacks, stat='min', na.rm=TRUE)
ndvi_sd <- cellStats(ndvi_filtered_stacks, stat='sd', na.rm=TRUE)
ndvi_median <- cellStats(ndvi_filtered_stacks, stat='median', na.rm=TRUE)

#graficas de los estadisticos 
stats <- data.frame(
  Date=substr(names(ndvi_filtered_stacks), 2,11),
  NDVI_Max = ndvi_max,
  NDVI_Min = ndvi_min,
  NDVI_Median = ndvi_median,
  NDVI_mean = ndvi_mean,
  NDVI_SD = ndvi_sd
)

stats$Date <- lubridate::ymd(stats$Date) # reformat the date
variables = xts::xts(x=stats[,-1], order.by=stats$Date) # Elegir variables, x sin dates
dygraph(variables) 

# exportar en csv
stats_name <- file.path(outDir, "NDVI_Statistics.csv")
write.csv(stats,stats_name)

# exportar geotiffs

for (i in 1:length(ndvi_filtered)){
  output_name <- file.path(outDir, paste0(names(ndvi_filtered[[i]]), "_NDVI.tif"))
  writeRaster(ndvi_filtered[[i]], output_name, overwrite = TRUE)
}

