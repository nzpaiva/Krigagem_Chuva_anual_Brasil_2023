
#inicio 05/02/2026
#
###########################################
# 1. Dados meteorologicos e contorno do Brasil:
##########################################
# Install devtools if not already installed
install.packages("devtools")
# Install BrazilMet
devtools::install_github("FilgueirasR/BrazilMet")
# Load the package for station data:
library(BrazilMet)

#Pacotes utilizados ao longo do programa:
if (!require("pacman")) install.packages("pacman")
pacman::p_load(geobr, sf, ggplot2, dplyr, gstat, gt, stars, 
               patchwork, viridis, sp, ggspatial, automap)


# Coleta os dados de todas as estações disponiveis no periodo:
df <- download_AWS_INMET_daily(stations = see_stations_info()$station_code,
                               start_date = "2024-01-01",
                               end_date = "2024-12-31")


# Cria um data frame da estação e a precipitação média no periodo selecionado
dat1 <- aggregate(rainfall_mm ~ station_code +
                    longitude_degrees +
                    latitude_degrees, data = df, mean)

# Transformar para sistema metrico 5880 ("SIRGAS 2000 / Brazil Polyconic")
dat1_sf <- st_as_sf(dat1, coords = c("longitude_degrees", "latitude_degrees"),
  crs = 4326)
dat1_m <- st_transform(dat1_sf, crs = 5880)
st_crs(dat1_sf)

# Coleta o contorno do brasil e faz a transformação para "SIRGAS 2000 / Brazil Polyconic"
contorno_brasil <- read_country(year=2020)
st_crs(contorno_brasil)
contorno_brasil_m <- st_transform(contorno_brasil, crs = 5880)
st_crs(contorno_brasil_m)

estados <- read_state(code_state = "all", year=2010)
st_crs(estados)
estados_m <- st_transform(estados, crs = 5880)

ggplot() +
  geom_sf(data = estados_m, fill = "white", color = "gray50") +
  geom_sf(data = dat1_m) +
  theme_minimal() + 
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Estações meteorológicas no Brasil")


#############################################################
# 2. Fazer a analise de variogramas neste ponto? REML Robusto ou automap
#############################################################

variograma_auto_ord <- autofitVariogram(rainfall_mm ~ 1, 
                                        input_data = as(dat1_m, "Spatial"))

preds <- variogramLine(variograma_auto_ord$var_model, maxdist = max(variograma_auto_ord$exp_var$dist))

ggplot() +
  geom_point(data = variograma_auto_ord$exp_var, aes(x = dist, y = gamma), size = 2.5) +
  geom_line(data = preds, aes(x = dist, y = gamma),
            color = "red", linewidth =.5) +
  labs(
    title = "Semivariograma empírico e modelo ajustado (Krig. ordinária)",
    subtitle = paste0(
      "Modelo: Ste",
      " | Pepita (Nugget): ", round(variograma_auto_ord$var_model$psill[1], 2),
      " | Patamar (Sill): ", round(variograma_auto_ord$var_model$psill[2], 2),
      " | Alcance (Range): ", round(variograma_auto_ord$var_model$range[2], 0)
    ),
    x = "Distância (m)",
    y = "Semivariância"
  ) +
  theme_classic()

ggsave("Escrita\\Gráficos\\semivariograma_ord.png", scale = 0.8)


##########################################
# 3. Aplicação da Krigagem ordinaria (KO)
##########################################

# Criar Grid de Predição, o cellsize determina a precisão
# Precisamos definir ONDE queremos estimar:
res = 20000 #tamanho de cada celula em metros
grid <- st_make_grid(contorno_brasil_m, cellsize = res, what = "polygons") |>
  st_as_sf() |>
  st_intersection( contorno_brasil_m)

krigagem_ord <- predict(
  gstat(
    formula = rainfall_mm ~ 1,
    data = dat1_m,
    model = variograma_auto_ord$var_model,
    nmax = 40
  ),
  newdata = grid
)

# var1.pred = Valor Estimado 
# var1.var = Variância de Krigagem (Erro)
# Mapa 1: Predição
ggplot() +
  geom_sf(data = krigagem_ord, aes(fill = var1.pred)) +
  geom_sf(data = contorno_brasil_m, fill = NA, color = "black", size = 0.5) + # Contorno
  scale_fill_gradientn(
    colours = c("#e0f3ff", "#a6d8ff", "#4ea5ff", "#0066cc", "#003f7f"),
    name = "Valor"
  ) +
  labs(title = "Predição de precipitação média anual (mm)") +
  theme_minimal() + #para colocar coordenadas
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank()) +
        annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black", "white")) +
        annotation_north_arrow(location = "tl", which_north = "true", 
                               pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                               style = north_arrow_fancy_orienteering) +
  
# Mapa 2: Variância (Erro)
ggplot() +
  geom_sf(data = krigagem_ord, aes(fill = var1.var)) +
  geom_sf(data = contorno_brasil_m, fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Erro") +
  labs(title = "Variância da predição") +
  theme_minimal() + #para colocar coordenadas
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank()) +
        annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black", "white")) +
        annotation_north_arrow(location = "tl", which_north = "true", 
                               pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                               style = north_arrow_fancy_orienteering)


ggsave("Gráficos\\ajuste_krig_ord.png")


###############################################
# 4. Aplicação da Krigagem com deriva externa (KED)
###############################################

###################################################
# Inter_step. Inicialização do rgee:
###################################################
# general libraries
install.packages("pacman")
library(pacman)

install.packages(c("remotes", "devtools", "googledrive"))
p_load(rgee, geojsonio, remotes, reticulate, devtools, googledrive)

rgee_environment_dir = "C:\\Users\\nzpai\\miniconda3\\envs\\rgee_py"

install_github("r-spatial/rgee")

## sometimes at this point you are required to restart R or the computer before proceeding
## try restarting if the installation do not finish properly and run the installation again after restart

# set python
reticulate::use_python(rgee_environment_dir, required=T)
rgee::ee_install_set_pyenv(
  py_path = rgee_environment_dir, # Change it for your own Python PATH
  py_env = "rgee_py" # Change it for your own Python ENV
)
Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

# Initialize the Python Environment
# to clean credentials: ee_clean_credentials()
rgee::ee_Initialize(drive = T)

## It worked if some text about google drive credentials appeared, and asked you to log in your GEE account.
## Congrats.

########### Extração dos dados de altitude: NASA/NASADEM_HGT/001

brasil_ee <- sf_as_ee(contorno_brasil_gps)  

modis <- ee$ImageCollection("NASA/NASADEM_HGT/001")$
  filterDate("2023-01-01", "2023-12-31")$
  select("Percent_Tree_Cover")

#extrai a media do periodo selecionado: imagino eu?
ndvi_mean <- modis$mean()

#extrai o ndvi para o terrotiorio brasileiro:
task <- ee_image_to_drive(
  image = ndvi_mean,
  description = "NDVI_BRASIL_2023",
  folder = "GEE",
  fileFormat = "GeoTIFF",
  region = brasil_ee$geometry(),
  scale = 250,
  maxPixels = 400000000
)
#isso gera um arquivo no google drive, preciso baixa-lo para esse ambiente
task$start()

#########################################################################

########## Extração dos dados a vegetação: ##############################

brasil_ee <- sf_as_ee(contorno_brasil_gps)  

modis <- ee$ImageCollection("MODIS/061/MOD44B")$
  filterDate("2023-01-01", "2023-12-31")$
  select("Percent_Tree_Cover")

#extrai a media do periodo selecionado: imagino eu?
ndvi_mean <- modis$mean()

#extrai o ndvi para o terrotiorio brasileiro:
task <- ee_image_to_drive(
  image = ndvi_mean,
  description = "NDVI_BRASIL_2023",
  folder = "GEE",
  fileFormat = "GeoTIFF",
  region = brasil_ee$geometry(),
  scale = 250,
  maxPixels = 400000000
)
#isso gera um arquivo no google drive, preciso baixa-lo para esse ambiente
task$start()


library(terra)

#carrega o raster no R:
ndvi <- rast("NDVI_BRASIL_2023_2026_02_17_14_52_02.tif")

#reprojeta em um sistema metrico de metros:
ndvi_m   <- terra::project(ndvi, "EPSG:5880")

# buffer de 10km metros ao redor dos pontos
buf <- buffer(vect(dat1_m), width = res/2)

# extrair média do NDVI dentro de cada buffer
ndvi_mean <- terra::extract(ndvi_m, buf, fun = mean, na.rm = TRUE)

# adicionar ao data.frame
dat1_m$ndvi_mean <- ndvi_mean[,2]
#####################################################################

# Ajusta o modelo de variograma com deriva externa:

variograma_auto_kde <- autofitVariogram(rainfall_mm ~ ndvi_mean, 
                                        input_data = as(dat1_m, "Spatial"))


preds <- variogramLine(variograma_auto_kde$var_model, maxdist = max(variograma_auto_kde$exp_var$dist))

ggplot() +
  geom_point(data = variograma_auto_kde$exp_var, aes(x = dist, y = gamma), size = 2.5) +
  geom_line(data = preds, aes(x = dist, y = gamma),
            color = "red", linewidth =.5) +
  labs(
    title = "Semivariograma empírico e modelo ajustado (Krig. deriva externa)",
    subtitle = paste0(
      "Modelo: Ste",
      " | Pepita (Nugget): ", round(variograma_auto_kde$var_model$psill[1], 2),
      " | Patamar (Sill): ", round(variograma_auto_kde$var_model$psill[2], 2),
      " | Alcance (Range): ", round(variograma_auto_kde$var_model$range[2], 0)
    ),
    x = "Distância (m)",
    y = "Semivariância"
  ) +
  theme_classic()

ggsave("Escrita\\Gráficos\\semivariograma_kde.png", scale = 0.8)

# Questão da resolução do grid:
# os dados do ndvi tem uma resolução muito superior ao grid de predição criado
grid_pts <- st_centroid(grid)

# buffer de 10km metros ao redor dos centros das celulas
buf_malha <- buffer(vect(grid_pts), width = res/2) 

# Extrair média do NDVI dentro de cada buffer
grid$ndvi_mean <- terra::extract(ndvi_m, buf_malha, fun = mean, na.rm = TRUE)

grid$ndvi_mean <- grid$ndvi_mean$Percent_Tree_Cover 

# Realiza a krigagem com deriva externa: 
kde <- predict(
  gstat(
    formula = rainfall_mm ~ ndvi_mean,
    data = dat1_m,
    model = variograma_auto_kde$var_model,
    nmax = 40
  ),
  newdata = grid
)

# var1.pred = Valor Estimado 
# var1.var = Variância de Krigagem (Erro)
# Mapa 1: Predição
ggplot() +
  geom_sf(data = kde, aes(fill = var1.pred)) +
  geom_sf(data = contorno_brasil_m, fill = NA, color = "black", size = 0.5) + 
  scale_fill_gradientn(
    colours = c("#e0f3ff", "#a6d8ff", "#4ea5ff", "#0066cc", "#003f7f"),
    name = "Valor"
  ) +
  labs(title = "Predição de precipitação média anual (mm)") +
  theme_minimal() + #para colocar coordenadas
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank()) +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  
# Mapa 2: Variância (Erro)
ggplot() +
  geom_sf(data = kde, aes(fill = var1.var)) +
  geom_sf(data = contorno_brasil_m, fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Erro") +
  labs(title = "Variância da predição") +
  theme_minimal() + #para colocar coordenadas
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank()) +
  annotation_scale(location = "bl", width_hint = 0.3, bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)


ggsave("Gráficos\\ajuste_krig_deriva_externa.png")

########################################################
# Aplicar a validação cruzada aos modelos
# Os dados de vegetação melhoram a previsão?
########################################################

#Executar Validação Cruzada (Leave-One-Out)
# A função krige.cv faz o loop automaticamente.

  # Aplicada a krigagem ordinaria e com devira externa:
  
  cv_results_ord <- krige.cv(rainfall_mm ~ 1, 
                         locations = dat1_m, 
                         model = variograma_auto_ord$var_model, debug.level = 0)


  cv_results_kde <- krige.cv(rainfall_mm ~ ndvi_mean, 
                        locations = dat1_m, 
                         model = variograma_auto_kde$var_model, debug.level = 0)

  # Estatísticas de Diagnóstico:

  # MSDR (Razão de Desvio Quadrático Médio)
  msdr_ord <- mean(cv_results_ord$zscore^2)
  msdr_kde <- mean(cv_results_kde$zscore^2)
  # ME (Erro Médio) - Deve ser próximo de 0 (viés)
  me_ord <- mean(cv_results_ord$residual)
  me_kde <- mean(cv_results_kde$residual)
  # RMSE (Raiz do Erro Quadrático Médio)
  rmse_ord <- sqrt(mean(cv_results_ord$residual^2))
  rmse_kde <- sqrt(mean(cv_results_kde$residual^2))
  
  resultados <- tibble(
    Indicador = c("Mean Error (Viés)", 
                  "RMSE (Precisão)", 
                  "MSDR (Calibração)"),
    "Krigagem ordinária" = round(c(me_ord, rmse_ord, msdr_ord), 4),
    "Krigagem com deriva externa" = round(c(me_kde, rmse_kde, msdr_kde), 4)
  )

  resultados |>
    knitr::kable()
  
  tab <- resultados 
  
  saveRDS(tab, "Escrita\\tabela_1.rds")
  
  # Calculo do coeficiente de determinação R^2:
  r2_ord <- 1-(sum((cv_results_ord$observed - cv_results_ord$var1.pred)**2)/sum((cv_results_ord$observed - mean(cv_results_ord$observed))**2)) ; r2_ord
  
  r2_kde <- 1-(sum((cv_results_kde$observed - cv_results_kde$var1.pred)**2)/sum((cv_results_kde$observed - mean(cv_results_kde$observed))**2)) ; r2_kde
  
  p1 <- ggplot(cv_results_ord, aes(x = observed, y = var1.pred)) +
    geom_point(color = "black", alpha = 0.5) +
    geom_abline(color = "red", slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Acurácia: Observado vs Predito",
         subtitle = paste0(
           "R^2: ", round(r2_ord, 3)
         ),
         x = "Precipitação observada", y = "Precipitação predita") +
    theme_classic()
    
  # Histograma dos Z-Scores (Deve parecer uma Normal(0,1))
  p2 <- ggplot(cv_results_ord, aes(x = zscore)) +
    geom_histogram(aes(y = ..density..), fill = "steelblue", bins = 20, color = "white", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", linewidth = 1) +
     labs(title = "Calibração: Z-Scores",
         subtitle = paste("MSDR:", round(msdr_ord, 3), "(Ideal = 1.0)"),
         x = "Resíduo Padronizado", y = "Densidade") +
    theme_classic()
  
  p1 + p2

  ggsave("Escrita\\Gráficos\\validação_cruzada_ord.png", scale = 0.8)
  
  ########################################
  # Aplicada a krigagem com deriva externa:
  ######################################## 
  
  p1 <- ggplot(cv_results_kde, aes(x = observed, y = var1.pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Acurácia: Observado vs Predito",
         subtitle = paste0("R^2: ", round(r2_kde, 3)),
        x = "Precipitação observada", y = "Precipitação predita") +
    theme_classic()
  
  # Histograma dos Z-Scores (Deve parecer uma Normal(0,1))
  p2 <- ggplot(cv_results_kde, aes(x = zscore)) +
    geom_histogram(aes(y = ..density..), bins = 20, fill = "steelblue", color = "white", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", size = 1) +
    labs(title = "Calibração: Z-Scores",
         subtitle = paste("MSDR:", round(msdr_kde, 3), "(Ideal = 1.0)"),
         x = "Resíduo Padronizado", y = "Densidade") +
    theme_classic()
  
  p1 + p2
  
  ggsave("Escrita\\Gráficos\\validação_cruzada_kde.png", scale = 0.8)
  
  
  
  fim <- proc.time() - inicio
  print(fim)
  
  
  