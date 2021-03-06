---
title: "MS-GWR Analyses with the Italian Data"
author:
  - Nicolas Kuehn^[Unversity of California, Los Angeles, kuehn@ucla.edu]
date: "14 September, 2021"
output:
  html_document:
  #pdf_document:
    keep_md: true
    toc: true
    toc_depth: 2
    number_sections: true
    highlight: tango
link-citations: yes
linkcolor: blue
citecolor: blue
urlcolor: blue
bibliography: /Users/nico/BIBLIOGRAPHY/BIBTEX/references.bib
---
  
  
  

# Introduction

This code contains the code to calculate the MS-GWR analyses.
The code mostly follows mostly @Caramenti2020, available at <{https://github.com/lucaramenti/ms-gwr/>.

# Getting Started


```r
# load required packages
library(lme4)
library(sf)
library(raster)
library(rgdal)
library(GWmodel)
library(psych)
library(ggplot2)
library(tictoc)

source(file.path("R_Caramenti", "functions.R"))
source(file.path("R_Caramenti", "predict_new.R"))
```


```r
# read in data file with event and station ids
# corrected for Vs30 differences between stations of same location
data <- read.csv(file.path('DATA', 'italian_data_pga_id_utm_stat.csv'))
print(dim(data))
```

```
## [1] 4784   16
```

```r
# Set linear predictors
mh = 5.5
mref = 5.324
h = 6.924
attach(data)
b1 = (mag-mh)*(mag<=mh)
b2 = (mag-mh)*(mag>mh)
c1 = (mag-mref)*log10(sqrt(JB_complete^2+h^2))
c2 = log10(sqrt(JB_complete^2+h^2))
c3 = sqrt(JB_complete^2+h^2)
f1 = as.numeric(fm_type_code == "SS")
f2 = as.numeric(fm_type_code == "TF")
k = log10(vs30/800)*(vs30<=1500)+log10(1500/800)*(vs30>1500)
y = log10(rotD50_pga)
detach(data)

eq <- data$EQID
stat <- data$STATID
n_rec <- length(b1)
n_eq <- max(eq)
n_stat <- max(stat)

# spatial coordinates
shape_utm = st_read(file.path("DATA",'confini_ut33.shp'))
```

```
## Reading layer `confini_ut33' from data source 
##   `/Users/nico/GROUNDMOTION/PROJECTS/NONERGODIC/ITALY/Git/NonErgGMM_Italy/DATA/confini_ut33.shp' 
##   using driver `ESRI Shapefile'
## Simple feature collection with 1 feature and 2 fields
## Geometry type: MULTIPOLYGON
## Dimension:     XY
## Bounding box:  xmin: -158718.7 ymin: 3930217 xmax: 800075.7 ymax: 5219220
## Projected CRS: WGS 84 / UTM zone 33N
```

```r
shape_utm_no_lamp = shape_utm
shape_utm_no_lamp$geometry[[1]][[12]] = NULL

shape_utm_no_sardinia = shape_utm
shape_utm_no_sardinia$geometry[[1]][[3]] = NULL
for (i in 4:10){
  shape_utm_no_sardinia$geometry[[1]][[4]] = NULL
}
for (i in 12:57){
  shape_utm_no_sardinia$geometry[[1]][[5]] = NULL
}

#Project coordinates of the events to work with UTM33:
latitude_ev = data$ev_latitude
longitude_ev = data$ev_longitude
long_lat_ev = cbind(longitude_ev, latitude_ev)
utm_ev = project(long_lat_ev, "+proj=utm +zone=33 ellps=WGS84") 
utm_ev = as.data.frame(utm_ev)
long_lat_ev = as.data.frame(long_lat_ev)

#Project coordinates of the stations to work with UTM33:
latitude_st = data$st_latitude
longitude_st = data$st_longitude
long_lat_st = cbind(longitude_st, latitude_st)
utm_st = project(long_lat_st, "+proj=utm +zone=33 ellps=WGS84") 
utm_st = as.data.frame(utm_st)
long_lat_st = as.data.frame(long_lat_st)

#Build spatial points data frame for UTM33 coordinates:
utm_ev_sp = SpatialPointsDataFrame(utm_ev, data[,1:7])
utm_st_sp = SpatialPointsDataFrame(utm_st, data[,1:7])

#Transform shapefile into spatial dataframe:
shape_utm_spatial = as_Spatial(shape_utm_no_sardinia)
grid_utm = makegrid(shape_utm_spatial, cellsize = 10000) # cellsize in map units!
grid_utm = SpatialPoints(grid_utm, proj4string = CRS(proj4string(shape_utm_spatial)))
```

```
## Warning in proj4string(shape_utm_spatial): CRS object has comment, which is lost
## in output
```

```r
grid_inside_utm = grid_utm[shape_utm_spatial,]
coords_utm = grid_inside_utm@coords

#Transform shapefile into spatial dataframe:
shape_utm_spatial = as_Spatial(shape_utm_no_sardinia)
grid_utm = makegrid(shape_utm_spatial, cellsize = 10000) # cellsize in map units!
grid_utm = SpatialPoints(grid_utm, proj4string = CRS(proj4string(shape_utm_spatial)))
```

```
## Warning in proj4string(shape_utm_spatial): CRS object has comment, which is lost
## in output
```

Set design matrices for constant, event, and site related coefficients.


```r
Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

# best bandwidth is set
bwe = 25000
bws = 75000
coords_df_utm = as.data.frame(coords_utm)
```


```r
tic()
sec = SEC_only_calibration(Xc, Xe, Xs, y, "c", bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
toc()
```

```
## 11434.796 sec elapsed
```

```r
#Now we can compute all the regression coefficients:
result = SEC_grid_creation(Xc, Xe, Xs, y,"c", bwe, bws, coordinates(utm_ev_sp),
                           coordinates(utm_st_sp), coords_utm, sec)
beta_const = result$beta_c

beta_k = t(result$beta_s)
beta_k_coord = cbind(coords_utm, beta_k)
beta_k_coord = as.data.frame(beta_k_coord)

beta_c2 = result$beta_e[1,]
beta_c2_coord = cbind(coords_utm, beta_c2)
beta_c2_coord = as.data.frame(beta_c2_coord)

beta_c3 = result$beta_e[2,]
beta_c3_coord = cbind(coords_utm, beta_c3)
beta_c3_coord = as.data.frame(beta_c3_coord)


save(sec, result, file = file.path('RESULTS', 'results_caramenti.Rdata'))
```

Calculate prediction for the data set.


```r
pcoords = as.data.frame(cbind(coordinates(utm_ev_sp), coordinates(utm_st_sp)))
tic()
prediction_msgwr <- emgwr_prediction_points_no_param2(Xc, Xe, Xs, y, sec, "sec", bwe, bws,
                                                      coordinates(utm_ev_sp), coordinates(utm_st_sp), Xc, Xe, Xs,
                                                      pcoords, 0.05)
toc()
```

```
## 16630.961 sec elapsed
```

```r
save(prediction_msgwr, file = file.path('RESULTS', 'residuals_caramenti.Rdata'))
```

Calculate prediction for one event in Middle Italy, with directions to North and South.


```r
data_pred <- read.csv(file.path("DATA", "data_prediction_line.csv"))

attach(data_pred)
b1p = (mag-mh)*(mag<=mh)
b2p = (mag-mh)*(mag>mh)
c1p = (mag-mref)*log10(sqrt(JB_complete^2+h^2))
c2p = log10(sqrt(JB_complete^2+h^2))
c3p = sqrt(JB_complete^2+h^2)
f1p = as.numeric(fm_type_code == "SS")
f2p = as.numeric(fm_type_code == "TF")
kp = log10(vs30/800)*(vs30<=1500)+log10(1500/800)*(vs30>1500)
yp = log10(rotD50_pga)
detach(data_pred)

data_pred_reg <- data.frame(Y = yp,
                            M1 = b1p,
                            M2 = b2p,
                            MlnR = c1p,
                            lnR = c2p,
                            R = c3p,
                            Fss = f1p,
                            Frv = f2p,
                            lnVS = kp,
                            #eq = data_pred$EQID,
                            eq = NA,
                            stat = data_pred$STATID,
                            intercept = 1,
                            resid = NA
)
data_pred_reg$idx_cell <- (n_rec + 1):(n_rec + length(data_pred_reg$Y))
data_pred_reg$idx <- (n_rec + 1):(n_rec + length(data_pred_reg$Y))

co_eq_pred <- data_pred[,c(13,14)]
co_stat_pred <- data_pred[,c(15,16)]

co_eq_stat_pred <- cbind(rowMeans(cbind(co_eq_pred$X_ev, co_stat_pred$X_stat)),
                         rowMeans(cbind(co_eq_pred$Y_ev, co_stat_pred$Y_stat)))


###############################################################################
co_eq_pred_msgwr <- co_eq_pred * 1000
co_stat_pred_msgwr <- co_stat_pred * 1000

#Build spatial points data frame for UTM33 coordinates:
utm_ev_pred_sp = SpatialPointsDataFrame(co_eq_pred_msgwr, data_pred[,1:6])
utm_st_pred_sp = SpatialPointsDataFrame(co_stat_pred_msgwr, data_pred[,1:6])

pcoords = cbind(co_eq_pred_msgwr, co_stat_pred_msgwr)

Xc_pred = cbind(b1p,b2p,f1p,f2p,c1p)
Xe_pred = cbind(c2p,c3p)
Xs_pred = kp

#load(file = file.path('RESULTS', 'results_caramenti.Rdata'))
prediction = emgwr_prediction_points_no_param2(Xc, Xe, Xs, y, sec, "sec", bwe, bws,
                                               coordinates(utm_ev_sp), coordinates(utm_st_sp), Xc_pred, Xe_pred, Xs_pred,
                                               pcoords, 0.05)
save(prediction,
     file = file.path('RESULTS', 'results_pred_mswgr_line.Rdata'))
```

Calculate predictions for three events, with four directions (East, West, South, North) each.


```r
data_pred <- read.csv(file.path("DATA", "data_prediction.csv"))

attach(data_pred)
b1p = (mag-mh)*(mag<=mh)
b2p = (mag-mh)*(mag>mh)
c1p = (mag-mref)*log10(sqrt(JB_complete^2+h^2))
c2p = log10(sqrt(JB_complete^2+h^2))
c3p = sqrt(JB_complete^2+h^2)
f1p = as.numeric(fm_type_code == "SS")
f2p = as.numeric(fm_type_code == "TF")
kp = log10(vs30/800)*(vs30<=1500)+log10(1500/800)*(vs30>1500)
yp = log10(rotD50_pga)
detach(data_pred)

data_pred_reg <- data.frame(Y = yp,
                            M1 = b1p,
                            M2 = b2p,
                            MlnR = c1p,
                            lnR = c2p,
                            R = c3p,
                            Fss = f1p,
                            Frv = f2p,
                            lnVS = kp,
                            #eq = data_pred$EQID,
                            eq = NA,
                            stat = data_pred$STATID,
                            intercept = 1,
                            resid = NA
)
data_pred_reg$idx_cell <- (n_rec + 1):(n_rec + length(data_pred_reg$Y))
data_pred_reg$idx <- (n_rec + 1):(n_rec + length(data_pred_reg$Y))

co_eq_pred <- data_pred[,c(13,14)]
co_stat_pred <- data_pred[,c(15,16)]

co_eq_stat_pred <- cbind(rowMeans(cbind(co_eq_pred$X_ev, co_stat_pred$X_stat)),
                         rowMeans(cbind(co_eq_pred$Y_ev, co_stat_pred$Y_stat)))


###############################################################################
co_eq_pred_msgwr <- co_eq_pred * 1000
co_stat_pred_msgwr <- co_stat_pred * 1000

#Build spatial points data frame for UTM33 coordinates:
utm_ev_pred_sp = SpatialPointsDataFrame(co_eq_pred_msgwr, data_pred[,1:6])
utm_st_pred_sp = SpatialPointsDataFrame(co_stat_pred_msgwr, data_pred[,1:6])

pcoords = cbind(co_eq_pred_msgwr, co_stat_pred_msgwr)

Xc_pred = cbind(b1p,b2p,f1p,f2p,c1p)
Xe_pred = cbind(c2p,c3p)
Xs_pred = kp

#load(file = file.path('RESULTS', 'results_caramenti.Rdata'))
prediction = emgwr_prediction_points_no_param2(Xc, Xe, Xs, y, sec, "sec", bwe, bws,
                                               coordinates(utm_ev_sp), coordinates(utm_st_sp), Xc_pred, Xe_pred, Xs_pred,
                                               pcoords, 0.05)
save(prediction,
     file = file.path('RESULTS', 'results_pred_mswgr.Rdata'))
```

# References
