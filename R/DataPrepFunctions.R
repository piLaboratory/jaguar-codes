## Load packages (check if all those are necessary for this script)
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load("maptools",'move',"circular","RCurl","dplyr","readr","caTools","adehabitatLT","rgl",
                           "lubridate","raster","amt","tibble","ezknitr","lattice","rgdal","sp")

## Add 2000 to years
get.year <- function(time.stamp) {
  init <- gregexpr('/', time.stamp, fixed = T)[[1]][2] + 1
  end <- gregexpr(' ', time.stamp, fixed = T)[[1]][1] - 1
  substr(time.stamp, init, end)
}

## New dates
set.year <- function(time.stamp, year) {
  init <- gregexpr('/', time.stamp, fixed = T)[[1]][2]
  end <- gregexpr(' ', time.stamp, fixed = T)[[1]][1]
  paste(substr(time.stamp, 1, init), year,
        substr(time.stamp, end, nchar(time.stamp)), sep = "")
}

#' Converts coordinate systems (CRC)
#' @param data a dataframe with track points, where x and y coordinates must be names 'x' and 'y'
#' @param crs.input string with the descripition of the CRS used in data
#' @param crs.output string of the desired CRS output
#' @param point.names string vector of length 4; names of for the columns with the original and transformed x and y coordinates.
#' @return a dataframe like data with the converted coordinates added as the last two columns.
crs.convert <- function(data, crs.input="+proj=longlat +datum=WGS84", crs.output,
                        point.names = c("converted_x","converted_y","original_x","original_y")) {
    ## Projection in the input CRS
    coord.latlong <- SpatialPoints(cbind(data$x,data$y), proj4string = CRS(crs.input))
    ## Transforming coordinates to the output CRS 
    coord.UTM <- spTransform(coord.latlong , CRS(crs.output))
    coord.latlong.df <- as.data.frame(coord.latlong)
    coord.utm.df <- as.data.frame(coord.UTM)
    locsj_matx <- cbind(coordinates(coord.UTM), coordinates(coord.latlong))
    locsj_df <- as.data.frame(locsj_matx)
    colnames(locsj_df) <- point.names
    return(cbind(data, locsj_df))
}

#' Converts a dataframe with movement data to trk object
#' @param data dataframe with movement data, each line a location point
#' @param \ldots further arguments to be passed to 'mk_track' function (see amt::mk_track).
#' @note current version follows the commands of Fienberg's code that adds movement and time data to the output trk object; check if this is necessary.
#' @return a trk objetc with dir_abs, dir_rel, step length, nsd, week, month year and hour of each track location.
trk.convert <- function(data, ...){
    trk <- mk_track(data, ...)
    trk <- trk %>% arrange(id)
    trk <- trk %>% time_of_day(solar.dep = 18,
                               include.crepuscule = TRUE)  # this considers the Astronomical Twilight
    nesttrk<-trk%>%nest(-id)
    ## We can add a columns to each nested column of data using purrr::map
    trk<-trk %>% nest(-id) %>% 
        mutate(sl = map(data, step_lengths),
               nsd_=map(data, nsd),
               dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
               dir_rel = map(data, direction_rel))%>%unnest()
    trk.class<-class(trk)
    #' Calculate month, year, hour, week of each observation and append these to the dataset
    #' Unlike the movement charactersitics, these calculations can be done all at once, 
    #' since they do not utilize successive observations (like step lengths and turn angles do).
    trk<-trk%>% 
        mutate(
            week=week(t_),
            month = month(t_, label=TRUE), 
            year=year(t_),
            hour = hour(t_)
        )
    trk <- trk %>%dplyr::select(x_,y_,t_,id,tod_, everything())
    trk <- trk%>%dplyr::select(-sl,-project_region,-nsd_, -dir_abs,-dir_rel,
                               -week, -month, -year,-hour,-long_x,-lat_y, everything())
    ## Now, we need to again tell R that this is a track (rather than just a data frame)
    class(trk)<-trk.class
    return(trk)
}

