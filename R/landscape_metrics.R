#' ---
#' title: 'Calculating landscape metrics for jaguars in the Neotropics'
#' author: Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
#' ---

#' This script aims at calculating landscape metrics for the landscapes inhabited by
#' jaguars in the Neotropics. These jaguars have been monitored with GPS collars.

# packages
library(install.load)
install.load::install_load("tidyverse")
install.load::install_load("landscapemetrics", "landscapetools")

# data
env.data <- read.csv2(file = "maps/table_akde.csv", header = T) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(kernel = as.numeric(kernel),
                id = as.character(id) %>% ifelse(grepl("J", .), paste(., "per", sep = "_"), .))
  
jaguar.data <- read.csv2(file = "HR_year_2d.csv", header = T) %>% 
  dplyr::as_tibble() %>% 
  tidyr::pivot_longer(cols = starts_with("area_sq_km"), 
                      names_to = "kernel", 
                      names_prefix = "area_sq_km_",
                      values_to = "akde_area_sq_km") %>% 
  dplyr::mutate(kernel = as.numeric(kernel),
                id = as.character(id), 
                ID = as.numeric(as.character(ID)))
                      
jaguar.data$ID[(length(jaguar.data$ID)-3):length(jaguar.data$ID)] <- c(118, 118, 119, 119)

# merge
jaguar.env.data <- jaguar.data %>% 
  dplyr::left_join(env.data, by = c("id", "kernel")) %>% 
  dplyr::mutate(ID = as.numeric(as.character(ID))) %>% 
  dplyr::arrange(ID, kernel)

# check
jaguar.env.data %>% print(n = 300)

unique(jaguar.env.data$id)

# landscape metrics
maps.names <- list.files("maps", pattern = "forest_year", full.names = T)

maps <- list()
 
for(i in 1:nrow(jaguar.env.data)) {
  
  print(i)
  
  id <- str_remove(jaguar.env.data$id[i], pattern = "_per")
  kern <- jaguar.env.data$kernel[i]
  id_kern <- paste0("ind", id, "_", kern, sep = "")
  
  map.open <- maps.names %>% grep(pattern = id_kern, value = T)
  
  if(length(map.open) == 0) {
    maps[[i]] <- NULL
  } else {
    maps[[i]] <- raster::raster(map.open)
  }
  
}
names(maps) <- paste(jaguar.env.data$id, jaguar.env.data$kernel, sep = "_")

maps

la01 <- landscapetools::show_landscape(maps$`01_95`, discrete = TRUE) +
  scale_fill_manual(values = c("blue", "orange", "forestgreen")) +
  labs(title = "Paisagem 01") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = .5))
la01

# function to convert latlong to UTM zone
long2UTM <- function(rast) {
  center.long <- sum(raster::extent(rast)[1:2])/2
  zone <- (floor((center.long + 180)/6) %% 60) + 1
  
  ifelse(zone > 19, 
         paste0("+proj=utm +zone=", zone, " +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
         paste0("+proj=utm +zone=", zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
}

# test to find the utm zone
long2UTM(maps[[97]])

# for one landscape
metric <- raster::projectRaster(maps[[2]], crs = long2UTM(maps[[2]]))
metric[] <- ifelse(metric[] != 0, 1, 0)
plot(metric)

landscapemetrics::check_landscape(metric)

class_metrics <- landscapemetrics::list_lsm() %>%
  dplyr::filter(level == "class") %>% 
  dplyr::arrange(type)
class_metrics %>%  
  print(n = 100)

landscapemetrics::calculate_lsm(landscape = metric, 
                                what = c("lsm_c_area_mn", "lsm_c_area_sd", "lsm_c_area_cv",
                                         "lsm_c_pland", "lsm_c_te", "lsm_c_ed",
                                         "lsm_c_para_mn", "lsm_c_para_sd", "lsm_c_para_cv")) %>% 
  dplyr::filter(class == 1)

# for all landscapes
jaguar.env.data <- jaguar.env.data %>% 
  dplyr::mutate(patch_area_hectares_CV = NA,
                patch_area_hectares_avg = NA,
                patch_area_hectares_sd = NA,
                edge_density_m_sqm = NA,
                perimeter_area_ratio_m_sqm_CV = NA,
                perimeter_area_ratio_m_sqm_avg = NA,
                perimeter_area_ratio_m_sqm_sd = NA,
                p_landscape = NA,
                total_edge_m = NA)

for(i in 1:length(maps)) {
  
  print(i)
  
  if(!is.null(maps[[i]])) {
    metric <- raster::projectRaster(maps[[i]], crs = long2UTM(maps[[i]]))
    metric[] <- ifelse(metric[] != 0, 1, 0)
    
    vals <- landscapemetrics::calculate_lsm(landscape = metric, 
                                    what = c("lsm_c_area_mn", "lsm_c_area_sd", "lsm_c_area_cv",
                                             "lsm_c_pland", "lsm_c_te", "lsm_c_ed",
                                             "lsm_c_para_mn", "lsm_c_para_sd", "lsm_c_para_cv")) %>% 
      dplyr::filter(class == 1) %>% 
      pull(value)
    
    jaguar.env.data[i,20:ncol(jaguar.env.data)] <- vals
  } 
  
}

jaguar.env.data

# export
jaguar.env.data %>% 
  readr::write_csv2(path = "jaguar_env_data_complete.csv")
