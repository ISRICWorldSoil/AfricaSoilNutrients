## Prediction of soil nutrients N,P,K,S,Ca,Mg,B,Cu,Fe,Mn,Zn (Mehlich-3 extractable)
## by: Tom.Hengl@isric.org (with contributions by G.B.M. Heuvelink and M.G. Walsh)

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "ROCR", "randomForest", "R.utils", "plyr", "parallel", "psych", "mda", "dismo", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret", "plotKML", "compositions", "h2o", "bartMachine", "spatstat", "maptools", "maps", "grDevices", "ggbiplot", "factoextra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/AFMicroNutrients")
load(".RData")
library(plyr)
library(stringr)
library(sp)
library(foreign)
library(rgdal)
library(tools)
#library(devtools)
#devtools::install_github('dmlc/xgboost')
library(xgboost)
#devtools::install_github("imbs-hl/ranger", subdir = "ranger-r-package/ranger")
## https://github.com/imbs-hl/ranger/issues/121
library(ranger)
library(caret)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
#library(snowfall)
library(utils)
library(plotKML)
library(R.utils)
library(GSIF)
library(parallel)
library(doParallel)
library(maps)
library(scales)
library(compositions)
library(snowfall)
library(raster)
library(spatstat)
library(maptools)
library(stringr)
library(lattice)
library(ggplot2)
library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
#install_github("kassambara/factoextra")
library("factoextra")

plotKML.env(convert="convert", show.env=FALSE)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "gdal_translate"
  gdalwarp =  "gdalwarp"
  gdalbuildvrt = "gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
  gdaladdo = "gdaladdo"
}
system("/usr/local/bin/gdal-config --version")
## rJava install requires rJava locations
#Sys.setenv(JAVA_HOME='/usr/lib/jvm/default-java')
#system('sudo apt-get install r-cran-rjava')
#system('export LD_LIBRARY_PATH=/usr/lib/jvm/default-java/jre/lib/amd64:/usr/lib/jvm/default-java/jre/lib/amd64/default')
#system('sudo /usr/lib64/microsoft-r/3.3/lib64/R/bin/R CMD javareconf')

#source("/data/models/extract.equi7t3.R")
source("/data/models/wrapper.predict_cs.R")
source("/data/models/saveRDS_functions.R")
source("Africa_soilnutrients_functions.R")

##----------- IMPORT AND MERGING OF POINTS  -----------

## To convert Cmol/kg to ppm multiply by 1000 x atomic weight[4]. To convert % to ppm divide by 10,000.
#xAl <- Al/90
#xCa <- Ca/200
#xK <- K/391
#xMg <- Mg/121
#xNa <- Na/230
#xN <- N/14

#tvars <- c("NTOM3S", "PHOM3S", "KALM3S", "SULM3S", "CALM3S", "MAGM3S", "BORM3S", "CUTM3S", "FETM3S", "MNTM3S", "ZNTM3S")
t.vars <- c("C","N","P","P.B","P.O","P.T","K","S","Ca","Mg","Na","B","Cu","Fe","Mn","Zn","Al")

## IMPORT POINT DATA:
PBL_AfSP01302Qry <- read.csv("PBL_AfSP01302Qry.csv")
## "AvailP_B" - Available P based on Bray (most comparable to Mehlic3)
AfSP01302Qry <- rename(PBL_AfSP01302Qry, replace=c("OrgC"="C", "TotalN"="N", "AvailP_B"="P.B", "AvailP_O"="P.O", "TotalP"="P.T", "oFe"="Fe", "oMn"="Mn", "oZn"="Zn", "oCu"="Cu", "oB"="B", "oS"="S", "ExtrCa"="Ca", "ExtrMg"="Mg", "ExtrNa"="Na", "ExtrK"="K"))
AfSP01302Qry <- AfSP01302Qry[,c("C","N","P.B","P.O","P.T","K","S","Ca","Mg","Na","B","Cu","Fe","Mn","Zn","X_LonDD","Y_LatDD","LayerID","T_Year","UpDpth","LowDpth")]
for(j in names(AfSP01302Qry)){
  if(is.numeric(AfSP01302Qry[,j])){
    AfSP01302Qry[,j] <- ifelse(AfSP01302Qry[,j]==-9999, NA, AfSP01302Qry[,j])
  }
}
AfSP01302Qry$X_LonDD <- ifelse(AfSP01302Qry$X_LonDD==0&AfSP01302Qry$Y_LatDD==0, NA, AfSP01302Qry$X_LonDD)
AfSP01302Qry$Y_LatDD <- ifelse(AfSP01302Qry$X_LonDD==0&AfSP01302Qry$Y_LatDD==0, NA, AfSP01302Qry$Y_LatDD)
## Convert to PPM:
for(j in c("C","N")){  AfSP01302Qry[,j] <- AfSP01302Qry[,j] * 1000 }

summary(!is.na(AfSP01302Qry$P.O)); summary(!is.na(AfSP01302Qry$P.B)); summary(!is.na(AfSP01302Qry$P.T))
xyplot(AfSP01302Qry$P.O~AfSP01302Qry$P.B, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="P (Bray)", ylab="P (Olsen)", main=paste("Total points:", sum(!is.na(AfSP01302Qry$P))))
xyplot(AfSP01302Qry$P.T~AfSP01302Qry$P.B, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="P (Bray)", ylab="P (total)", main=paste("Total points:", sum(!is.na(AfSP01302Qry$P))))
## P.T is most reliable, P shows some groupings around number 1 and 5 --> looks to be OK
## Test "P" with model fitting

P.T.xy <- AfSP01302Qry[!is.na(AfSP01302Qry$P.T)&!is.na(AfSP01302Qry$X_LonDD),c("X_LonDD","Y_LatDD","LayerID","P.T")]
coordinates(P.T.xy) <- ~ X_LonDD+Y_LatDD
proj4string(P.T.xy) <- CRS("+proj=longlat +datum=WGS84")
#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#kml(P.T.xy, colour=log1p(P.T), size=P.T, balloon=TRUE, colour_scale=SAGA_pal[[1]], shape=shape, points_names=P.T.xy$LayerID)
## remove some suspicious points?
AfSP01302Qry <- AfSP01302Qry[!AfSP01302Qry$LayerID %in% c("ET_5322_5-083_2","ET_5322_5-083_3"),]
## Take out points from Guinea? (suspiciosly too high values)
#sel.GA <- grep(pattern=glob2rx("^GN_*"), paste(AfSP01302Qry$LayerID))
#AfSP01302Qry[sel.GA,c("P.B","P.O","P.T")] <- NA
AfSP01302Qry$SOURCEDB = "AfSPDB"

PBLTopsoilData <- read.csv("PBLTopsoilData.csv")									 	
PBLTopsoil <- rename(PBLTopsoilData, replace=c("OrgC_g.kg"="C", "TotalN_g.kg"="N", "ExtrP_ppm"="P", "ExtrFe_ppm"="Fe", "ExtrMn_ppm"="Mn", "ExtrZn_ppm"="Zn", "ExtrCu_ppm"="Cu", "ExtrB_ppm"="B", "ExtrS_ppm"="S", "ExtrCa_ppm"="Ca", "ExtrMg_ppm"="Mg", "ExtrNa_ppm"="Na", "ExtrK_ppm"="K", "ExtrAl_ppm"="Al"))
PBLTopsoil <- PBLTopsoil[,c("C","N","P","K","S","Ca","Mg","Na","B","Cu","Fe","Mn","Zn","Al","Longitude","Latitude","UpDpth","LowDpth","Year","ProfileID")]
for(j in names(PBLTopsoil)){
  if(is.numeric(PBLTopsoil[,j])){
    PBLTopsoil[,j] <- ifelse(PBLTopsoil[,j]==-9999, NA, PBLTopsoil[,j])
  }
}
PBLTopsoil$Longitude <- ifelse(PBLTopsoil$Longitude==0&PBLTopsoil$Latitude==0, NA, PBLTopsoil$Longitude)
PBLTopsoil$Latitude <- ifelse(PBLTopsoil$Longitude==0&PBLTopsoil$Latitude==0, NA, PBLTopsoil$Latitude)
## Convert to PPM:
for(j in c("C","N")){  PBLTopsoil[,j] <- PBLTopsoil[,j] * 1000 }
PBLTopsoil$SOURCEDB = "PBL_IFDC"
## plot in GE:
P.xy <- PBLTopsoil[!is.na(PBLTopsoil$P)&!is.na(PBLTopsoil$Longitude),c("Longitude","Latitude","ProfileID","P")]
coordinates(P.xy) <- ~ Longitude+Latitude
proj4string(P.xy) <- CRS("+proj=longlat +datum=WGS84")
#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#kml(P.xy, colour=log1p(P), size=P, balloon=TRUE, colour_scale=SAGA_pal[[1]], shape=shape, points_names=P.xy$ProfileID)
x = lapply(PBLTopsoil, function(x){sum(!is.na(x))})

## ----------- AfSIS point data -----------
afD <- read.csv("AfSIS_Micro60.csv")
afD2 <- read.csv("AfSIS_Samples.csv")
summary(afD2$N)
afL <- read.csv("AfSIS_Profiles.csv")
afALL <- read.csv("AfSIS_pr.csv")
afALL$Depth <- -100 * afALL$Depth
summary(afALL$xN)
afALL$N <- afALL$xN*10000
afALL$C <- afALL$C*10000
## cmolc/kg (xCa = Ca/200, xK = K/391, xMg = Mg/121, xNa = Na/230).
afALL$K <- afALL$xK*391
afALL$Ca <- afALL$xCa*200
afALL$Mg <- afALL$xMg*121
afALL$Na <- afALL$xNa*230
afALL$PID <- paste(afALL$Site, afALL$SSN, sep="_")
## remove calibration data:
afALL <- afALL[which(!afALL$SSN %in% afD2$SSN),]

## ----------- EthioSIS points -----------
etD <- read.csv("ET_Samples.csv")
etL <- read.csv("ET_Profiles.csv")
## plot next to each other:
summary(etD$K); summary(afD2$K)
histbackback(log1p(afALL$N), log1p(afD2$N), xlab=c("SS","AfSS"))
histbackback(log1p(afALL$C), log1p(afD2$C), xlab=c("SS","AfSS"))
histbackback(log1p(etD$K), log1p(afD2$K), xlab=c("ET","AfSS"))
histbackback(log1p(etD$S), log1p(afD2$S), xlab=c("ET","AfSS"))
histbackback(log1p(etD$B), log1p(afD$B), xlab=c("ET","AfSS"))
histbackback(log1p(etD$Ca), log1p(afD2$Ca), xlab=c("ET","AfSS"))
histbackback(log1p(etD$Zn), log1p(afD$Zn), xlab=c("ET","AfSS"))
histbackback(log1p(etD$Mg), log1p(afD2$Mg), xlab=c("ET","AfSS"))
histbackback(log1p(afALL$Mg), log1p(etD$Mg), xlab=c("AfSS-sp","ET"))
histbackback(log1p(etD$Cu), log1p(afD$Cu), xlab=c("ET","AfSS"))
histbackback(log1p(etD$P), log1p(afD2$P), xlab=c("ET","AfSS"))
## should all match
af <- join(join(afD, afD2, by="SSN"), afL, by="PID", type="left")
af$SOURCEDB = "AfSIS_I"
x = lapply(af, function(x){sum(!is.na(x))})
et <- join(etD, etL, by="PID")
et$SOURCEDB = "EthioSIS"
x = lapply(et, function(x){sum(!is.na(x))})

## Simulated points:
load("deserts.pnt.rda")
nut.sim <- as.data.frame(spTransform(deserts.pnt, CRS("+proj=longlat +datum=WGS84")))
#load("nut.sim.rda")
nut.sim[,1] <- NULL
nut.sim <- nut.sim[nut.sim$x > -20 & nut.sim$x < 53 & nut.sim$y > -37 & nut.sim$x < 37,]
nut.sim <- plyr::rename(nut.sim, c("x"="Lon", "y"="Lat"))
nut.sim$PID <- paste("Simulated", 1:nrow(nut.sim), sep="_")
## Keith: I agree that inserting pseudo-points based on expert knowledge is likely to improve predictions.
## insert zeros for all nutrients except for the once we are not sure:
## http://www.decodedscience.org/chemistry-sahara-sand-elements-dunes/45828
nut.sim[,t.vars[!t.vars %in% c("Ca","Mg","Na","Al")]] <- 0
str(nut.sim)
nut.sim$SOURCEDB = "Simulated"

##----------- Vital signs data-----------
vs_UG = read.csv("Eplot_UGA_with_Lat_Lon.csv")
## 1374 points
vs_UG.s = plyr::rename(vs_UG, c("latitude"="Lat", "longitude"="Lon", "SSN"="PID", "sample_depth_top..cm."="UpDpth", "sample_depth_bottom..cm."="LowDpth", "Nitrogen.content.for.acid.treated.sample.to.remove.carbonates....by.weight."="N", "Exchangeable.aluminium.concentration.by.Mehlich.3.extraction..m"="Al", "Copper.concentration.by.Mehlich.3.extraction..mg.kg..1."="Cu", "Iron.concentration.by.Mehlich.3.extraction..mg.kg..1."="Fe", "Exchangeable.Manganese.concentration.by.Mehlich.3.extraction..m"="Mn", "Phosphorus.by.Mehlich.3.extraction..mg.kg..1."="P", "Boron.concentration.by.Mehlich.3.extraction..mg.kg..1."="B", "Phosphorus.total.element.concentration.TXRF..mg.kg..1."="P.T", "Sulphur.by.Mehlich.3.extraction..mg.kg..1."="S", "Zinc.by.Mehlich.3.extraction..mg.kg..1."="Zn", "Exchangeable.Magnesium.by.wet.method..mg.kg..1."="Mg", "Exchangeable.Sodium.concentration.by.Mehlich.3.extraction..mg.k"="Na", "Potassium.concentration.by.Mehlich.3.extraction..mg.kg..1."="K", "Exchangeable.calcium.concentration.by.Mehlich.3.extraction..mg"="Ca"))
vs_UG.s$Year = as.numeric(format(as.Date(vs_UG.s$Sampling.date, format="%m/%d/%Y"),'%Y'))
vs_UG.s = vs_UG.s[,c("PID","Lat","Lon","Year","UpDpth","LowDpth","N","P","P.T","S","B","Cu","Fe","Mn","Zn","Al","Ca","Na","Mg","K")]
## convert values
vs_UG.s$N = vs_UG.s$N * 1e2
## Ca, Na, Mg and K have values in cmolc/kg:
vs_UG.s$K = vs_UG.s$K * 391
vs_UG.s$Ca = vs_UG.s$Ca * 200
vs_UG.s$Mg = vs_UG.s$Mg * 121
vs_UG.s$Na = vs_UG.s$Na * 230
#xCa <- Ca/200
#xK <- K/391
#xMg <- Mg/121
#xNa <- Na/230
vs_UG.s$SOURCEDB = "VitalSigns"
xC = lapply(vs_UG.s[,c("N","P","P.T","S","B","Cu","Fe","Mn","Zn","Al","Ca","Na","Mg","K")], summary)
names(xC) = c("N","P","P.T","S","B","Cu","Fe","Mn","Zn","Al","Ca","Na","Mg","K")

## Nigeria data AfSIS II:
## Markus Walsh: All the M3 values are in ppm. C&N in %, EC is in mS. All the samples are taken from 0-20 cm, compositing 4 samples per 100m^2 (see attached).
system('7za x top_MIR_pred_long.csv.zip')
nigeria = read.csv("top_MIR_pred_long.csv")
nigeria.s = plyr::rename(nigeria, c("lat"="Lat", "lon"="Lon", "ssid"="PID"))
nigeria.s$UpDpth = 0
nigeria.s$LowDpth = 20
nigeria.s$Depth = 10
nigeria.s$Year = as.numeric(format(as.Date(nigeria.s$date, format="%m/%d/%y"),'%Y'))
nigeria.s$Hp = NULL
nigeria.s$restrict = NULL
nigeria.s$plevel = NULL
nigeria.s$SSN = NULL
nigeria.s$depth = NULL
nigeria.s$date = NULL
nigeria.s$SOURCEDB = "AfSIS_II"
# shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
# nigeria.sp = nigeria[,c("ssid","lon","lat","pH","B")]
# coordinates(nigeria.sp) = ~lon+lat
# proj4string(nigeria.sp) = "+proj=longlat +datum=WGS84"
# kml(nigeria.sp, file="samples_nigeria_B.kml", folder.name="B", colour=B, shape=shape, points_names=signif(nigeria.sp$B,3), kmz=TRUE)
# kml(nigeria.sp, file="samples_nigeria_pH.kml", colour_scale=R_pal[["pH_pal"]], folder.name="pH", colour=pH, shape=shape, points_names=signif(nigeria.sp$pH,3), kmz=TRUE)

## Add/rename some columns:
AfSP01302Qry$Depth <- (AfSP01302Qry$LowDpth - AfSP01302Qry$UpDpth)/2 + AfSP01302Qry$UpDpth
PBLTopsoil$Depth <- (PBLTopsoil$LowDpth - PBLTopsoil$UpDpth)/2 + PBLTopsoil$UpDpth
summary(PBLTopsoil$Depth)
AfSP01302Qry <- rename(AfSP01302Qry, replace=c("Y_LatDD"="Lat","X_LonDD"="Lon","LayerID"="PID","T_Year"="Year"))
PBLTopsoil <- rename(PBLTopsoil, replace=c("Latitude"="Lat","Longitude"="Lon","ProfileID"="PID"))

## merge EVERYTHING together:
nut <- plyr::rbind.fill(list(AfSP01302Qry, PBLTopsoil, af[,-which(names(af) %in% c("SSN","ICRAFPID","Site","CC"))], afALL[,-which(names(afALL) %in% c("SSN","PHIHOX","Site","xCa","xK","xNa","EXB","xN","xMg"))], et[,-which(names(et) %in% c("Woreda"))], nut.sim, nigeria.s, vs_UG.s))
sel.rm.xy = is.na(nut$Lat)|is.na(nut$Lon)
nut = nut[-which(sel.rm.xy),]
str(nut)
## 119,329 obs. of  27 variables
#nut <- nut[,c("PID","Lat","Lon","Depth",t.vars)]
nut$Depth <- ifelse(is.na(nut$Depth), 15, nut$Depth)
summary(nut$Depth)
nut <- nut[nut$Depth>=0,]
nut$LOC_ID <- as.factor(paste(nut$Lon, nut$Lat, sep="_"))
summary(!duplicated(nut$LOC_ID))
## Total of 58,463 points
summary(!duplicated(nut$PID))
nut$PID <- make.unique(nut$PID)
## there should be no PID duplicates
nut <- plyr::rename(nut, c("Lon"="LONWGS84", "Lat"="LATWGS84", "Depth"="DEPTH"))
coordinates(nut) <- ~ LONWGS84+LATWGS84
proj4string(nut) <- CRS("+proj=longlat +datum=WGS84")
#plot(nut, pch="+")
save(nut, file="nut.rda")
#writeOGR(nut, "Africa_soil_nutrients.shp", "Africa_soil_nutrients", "ESRI Shapefile")
unlink("/data/GEOG/Nutrients/Africa_soil_nutrients.gpkg")
writeOGR(nut, "/data/GEOG/Nutrients/Africa_soil_nutrients.gpkg", "Africa_soil_nutrients", "GPKG")
## write per property:
for(j in t.vars){
  x <- nut[,c("PID","Year","UpDpth","LowDpth","DEPTH","SOURCEDB",j)]
  x <- x[!is.na(x@data[,j]),]
  unlink(paste0("Africa_soil_nutrients_", j, ".shp"))
  writeOGR(x, paste0("Africa_soil_nutrients_", j, ".shp"), paste0("Africa_soil_nutrients_", j), "ESRI Shapefile")
}

#shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#kml(nut, file="samples_af_K.kml", folder.name="Potassium", colour=log1p(K), shape=shape, points_names=nut$K, altitude=100-DEPTH, kmz=TRUE)
#kml(nut, file="samples_af_P.kml", folder.name="Phosphorus", colour=log1p(P), shape=shape, points_names=nut$P, altitude=100-DEPTH, kmz=TRUE)
#kml(nut, file="samples_af_N.kml", folder.name="Nitrogen", colour=log1p(N), shape=shape, points_names=nut$N, altitude=100-DEPTH, kmz=TRUE)
## scales:
sum.df <- data.frame(t(rbind(data.frame(lapply(nut@data[,t.vars], quantile, probs=c(0,0.01,0.5,0.99,1), na.rm=TRUE)), data.frame(lapply(nut@data[,t.vars], function(x){sum(!is.na(x))})))))
names(sum.df) = c("Min","Q_01","Q_50","Q_99","Max","N")
sum.df$N = unlist(lapply(nut@data[,t.vars], function(x){sum(!is.na(x))}))
write.csv(sum.df, file="Summary_values_nutrients.csv")
## P, P.O, P.T, Na, S, B, Cu, Zn = multiply by 100!
## Filter out some properties:
nut$C <- ifelse(nut$C>550000, NA, nut$C)
nut$B <- ifelse(nut$B<0, NA, nut$B)

## plot maps:
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
## c(-32.42, -42.90, 64.92, 41.08)
## Nitrogen
plot(country, col="darkgrey", ylim=c(-42.90, 41.08), xlim=c(-32.42, 64.92)) 
points(nut[!is.na(nut$N),], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
## Aluminum
plot(country, col="darkgrey", ylim=c(-42.90, 41.08), xlim=c(-32.42, 64.92)) 
points(nut[!is.na(nut$Al),], pch=21, bg=alpha("blue", 0.6), cex=.8, col="black")

## ------------- Grid definition -----------

af.prj = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
xllcorner = -4729500
yllcorner = -4328500
xurcorner = 3406000
yurcorner = 2602500
cellsize = 250
ncols = (xurcorner-xllcorner)/cellsize
## 32542
nrows = (yurcorner-yllcorner)/cellsize
## 27724
NODATA_value = -32768

## PREPARE MAPS
## Geology map (http://pubs.usgs.gov/of/1997/ofr-97-470/OF97-470A/):
surfgeo = readOGR("./geology/geo7_2ag.shp", "geo7_2ag")
summary(surfgeo$GLG)
surfgeo$GLG_INT = as.integer(surfgeo$GLG)
surfgeo = spTransform(surfgeo, CRS(af.prj))
writeOGR(surfgeo, "./geology/geo7_2ag_af.shp", "geo7_2ag_af", "ESRI Shapefile")
surfgeo.sum = data.frame(NAMES=levels(surfgeo$GLG), Value=1:length(levels(surfgeo$GLG)))
write.csv(surfgeo.sum, "surfacegeology_legend.csv")
system(paste0('/usr/local/bin/saga_cmd -c=56 grid_gridding 0 -INPUT \"./geology/geo7_2ag_af.shp\" -FIELD \"GLG_INT\" -GRID \"surfgeo_250m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
system(paste(gdal_translate, 'surfgeo_250m.sdat ./stacked250m/GESUSG6.tif -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\" -overwrite'))
unlink('surfgeo_250m.sdat')
## Coal fields:
coalpoly = readOGR("./geology/Africa_coal.shp", "Africa_coal")
proj4string(coalpoly) = "+proj=longlat +datum=WGS84"
summary(coalpoly$age)
coalpoly = spTransform(coalpoly, CRS(af.prj))
writeOGR(coalpoly, "./geology/Africa_coal_af.shp", "Africa_coal_af", "ESRI Shapefile")
system(paste0('/usr/local/bin/saga_cmd -c=56 grid_gridding 0 -INPUT \"./geology/Africa_coal_af.shp\" -FIELD \"age\" -OUTPUT 0 -GRID \"coalpoly_250m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
system(paste(gdal_translate, 'coalpoly_250m.sdat ./stacked250m/COAUSG6.tif -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\" -overwrite'))
unlink('coalpoly_250m.sdat')
## Mangroves:
system(paste0(gdalwarp, ' /data/mangroves/MANGPRf_500m.sdat ./stacked250m/MNGMAP.tif -ot \"Byte\" -dstnodata \"255\" -co \"COMPRESS=DEFLATE\" -t_srs \"',af.prj,'\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner, ' -tr 250 250 -overwrite'))

## Ground water maps (http://www.bgs.ac.uk/research/groundwater/international/africangroundwater/mapsDownload.html):
gwd.nlst = c("xyzASCII_dtwmap_v1.txt","xyzASCII_gwprod_v1.txt")
gwd.lst <- lapply(paste0("./geology/", gwd.nlst), read.table, header=TRUE)
summary(gwd.lst[[1]]$DTWAFRICA)
gwd.lst[[1]]$int = ifelse(gwd.lst[[1]]$DTWAFRICA=="VS", 3.5, ifelse(gwd.lst[[1]]$DTWAFRICA=="S", 16, ifelse(gwd.lst[[1]]$DTWAFRICA=="SM",37.5,ifelse(gwd.lst[[1]]$DTWAFRICA=="M",75,ifelse(gwd.lst[[1]]$DTWAFRICA=="D", 175,ifelse(gwd.lst[[1]]$DTWAFRICA=="VD",400,NA))))))
summary(gwd.lst[[2]]$GWPROD_V2)
gwd.lst[[2]]$int = ifelse(gwd.lst[[2]]$GWPROD_V2=="VH", 40, ifelse(gwd.lst[[2]]$GWPROD_V2=="H", 12.5, ifelse(gwd.lst[[2]]$GWPROD_V2=="M", 3, ifelse(gwd.lst[[2]]$GWPROD_V2=="LM", 0.75, ifelse(gwd.lst[[2]]$GWPROD_V2=="L", .3, ifelse(gwd.lst[[2]]$GWPROD_V2=="VL", .05, 0))))))
for(i in 1:length(gwd.nlst)){
  coordinates(gwd.lst[[i]]) = ~ X + Y
  gridded(gwd.lst[[i]]) = TRUE
  proj4string(gwd.lst[[i]]) = "+proj=longlat +datum=WGS84"
  writeGDAL(gwd.lst[[i]]["int"], paste0("./geology/", gsub("xyzASCII_", "", gsub(".txt", ".tif", gwd.nlst[i]))), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
}

## ----------- Density of mineral expolation sites-----------
#system('wget https://mrdata.usgs.gov/mrds/mrds-trim.zip /data/Geochem/mrds-trim.zip')
#system('7za x /data/Geochem/mrds-trim.zip')
mrds = readOGR("/data/Geochem/mrds-trim.shp", "mrds-trim")
## 304389 features
mrds = spTransform(mrds, CRS(af.prj))
## subset to Africa:
sel.mrds = mrds@coords[,1]>xllcorner & mrds@coords[,1]<xurcorner & mrds@coords[,2]>yllcorner & mrds@coords[,2]<yurcorner
## 2296 points
mrds = mrds[sel.mrds,]
mrds@data$CODE_LIST = droplevels(mrds@data$CODE_LIST)
mrds$mag = 1
summary(mrds@data$CODE_LIST)
## area of interest:
system(paste0(gdalwarp, ' /data/LDN/GAUL_COUNTRIES_300m_ll.tif GAUL_COUNTRIES_5km.tif -ot \"Int16\" -r \"near\" -co \"COMPRESS=DEFLATE\" -t_srs \"', af.prj, '\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner, ' -tr 5000 5000'))
mask = readGDAL("GAUL_COUNTRIES_5km.tif")
summary(mask$band1)
mask$band1 = ifelse(mask$band1<0, NA, mask$band1)
wowin <- as(mask, "owin")
CODE.lst = c("FE","CU","MN","MG","AL","ZN")

## run in parallel:
sfInit(parallel=TRUE, cpus=length(CODE.lst))
sfExport("kern.density", "mrds", "wowin", "CODE.lst")
sfLibrary(rgdal)
sfLibrary(spatstat)
sfLibrary(maptools)
x <- sfClusterApplyLB(CODE.lst, function(i){try( kern.density(i, input.points=mrds, wowin=wowin, COLUMN="CODE_LIST", sigma=80000, out.dir="./mrds", ext="5km") )})
sfStop()
## Phosphates only:
phosph = readOGR("/data/Geochem/phosphate/phosphate.shp", "phosphate")
phosph = spTransform(phosph, CRS(af.prj))
phosph$CODE_LIST = "phosphates"
kern.density(i="phosphates", input.points=phosph, wowin=wowin, COLUMN="CODE_LIST", sigma=80000, out.dir="./mrds", ext="5km")

## ----------- List of covariates -----------
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
ms.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
DEM.lst <- paste0("/data/AFsoilmapper/", c("VDPSRT6.tif", "TPISRT6.tif", "SLPSRT6.tif", "POSSRT6.tif", "NEGSRT6.tif", "MRNSRT6.tif", "DVMSRT6.tif", "DEM100m.tif"))
landsat.lst <- c(list.files("/data/Landsat", pattern=glob2rx("Water_*.tif"), full.names = TRUE), list.files("/data/Landsat", pattern=glob2rx("Landsat????_*.tif"), full.names = TRUE))
## MODIS EVI images:
EVI.lst <- paste0("/data/MOD13Q1/M_liste", ms.lst, ".vrt")
EVI2.lst <- paste0("/data/MOD13Q1/SD_liste", ms.lst, ".vrt")
## MODIS LST images:
LST.lst = c(paste0("/data/MOD11A2/M_liste", m.lst, "_D.vrt"), paste0("/data/MOD11A2/M_liste", m.lst, "_N.vrt"))
LST2.lst = c(paste0("/data/MOD11A2/SD_liste", m.lst, "_D.vrt"), paste0("/data/MOD11A2/SD_liste", m.lst, "_N.vrt"))
NEO.lst <- paste0("/data/NEO/SKY_WV_M_", ms.lst, "_10km.tif")
geog.lst <- c(paste0("./geology/", c("dtwmap_v1.tif", "gwprod_v1.tif")), paste0("./stacked250m/", c("COAUSG6.tif", "GESUSG6.tif", "MNGMAP.tif")))
MRDS.lst = list.files("./mrds", pattern=".tif", full.names=TRUE)
SoilGrids.lst = paste0("/data/GEOG/PHIHOX_M_sl", 1:7, "_250m_ll.tif")
## Total list all covs:
covs250m.lst = c(list.files(path="/data/EarthEnv", pattern=glob2rx("MODCF_*.tif$"), full.names = TRUE), list.files(path="/data/Climate", pattern=glob2rx("CHELSA_prec_*.tif$"), full.names = TRUE), geog.lst, "/data/DAAC/average_soil_and_sedimentary-deposit_thickness.tif", "/data/EarthEnv/giems_d15_v10.tif", "/data/ESA_global/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif", "/data/MCD43A4/desertPR_sin.tif", "/data/EcoTapestry/EF_LF_Desc_250m.tif", "/data/MDEM/MTWI_AF_250m.sdat", NEO.lst, DEM.lst, landsat.lst, EVI.lst, EVI2.lst, LST.lst, LST2.lst, MRDS.lst, SoilGrids.lst)
raster("/data/Groundwater/World_wtd_v2_f.tif")
raster("/data/Climate/CHELSA_prec_1979-2013_land.tif")
str(covs250m.lst)
## 135
## print out all problematic rasters
r.check = sapply(covs250m.lst, function(i){proj4string(raster(i))})
r.check[is.na(r.check)]
## /data/Climate/CHELSA_prec_1979-2013_land.tif

## ------------- OVERLAY / CREATE REGRESSION MATRIX -----------

## overlay in parallel TAKES 20 MINS:
sfInit(parallel=TRUE, cpus=56)
sfExport("nut", "covs250m.lst", "extract.tif")
sfLibrary(rgdal)
sfLibrary(raster)
m.out <- data.frame(sfClusterApplyLB(covs250m.lst, function(i){try( extract.tif(i, nut) )}))
sfStop()

names(m.out) = file_path_sans_ext(basename(covs250m.lst))
#head(m.out)
summary(as.factor(m.out$`ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1`))
summary(m.out$`CHELSA_prec_7_1979-2013`)
summary(m.out$`CHELSA_prec_2_1979-2013`)
m.out$desertPR_sin = NULL
summary(m.out$dtwmap_v1) ## misses Madagascar, hence better take out
m.out$dtwmap_v1 = NULL
## Filter out missing values in some covs:
#for(i in c(grep("CHELSA_prec", names(m.out)))){
#  m.out[,i] = ifelse(m.out[,i]<0, 0, m.out[,i])
#}
## remove layers with too many missing values?
c.check = sapply(m.out, function(i){sum(is.na(i))})
names(m.out)[!c.check<500]
#m.out = m.out[,c.check<500]
## "gwprod_v1"
summary(m.out$gwprod_v1)
names(m.out) = gsub("-", "_", names(m.out))
m.out$PID = nut@data$PID
## Convert land cover and geology to indicators
m.out$GESUSG6 = as.factor(m.out$GESUSG6)
summary(m.out$GESUSG6)
m.out$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1 = as.factor(m.out$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1)
summary(m.out$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1)
m.out$EF_LF_Desc_250m = as.factor(m.out$EF_LF_Desc_250m)
summary(m.out$EF_LF_Desc_250m)
geo.mat = data.frame(model.matrix(~GESUSG6-1, m.out))
xsum.geo = sapply(geo.mat, sum, na.rm=TRUE)
sel.rm.geo = names(geo.mat)[which(xsum.geo<6)]
lc.mat = data.frame(model.matrix(~ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1-1, m.out))
xsum.lc = sapply(lc.mat, sum, na.rm=TRUE)
sel.rm.lc = names(lc.mat)[which(xsum.lc<6)]
lf.mat = data.frame(model.matrix(~EF_LF_Desc_250m-1, m.out))
xsum.lf = sapply(lf.mat, sum, na.rm=TRUE)
sel.rm.lf = names(geo.mat)[which(xsum.lf<6)]

## merge everything together:
geo.mat = cbind(data.frame(PID=m.out$PID[as.numeric(rownames(geo.mat))]), geo.mat[,names(geo.mat)[which(xsum.geo>5)]])
lc.mat = cbind(data.frame(PID=m.out$PID[as.numeric(rownames(lc.mat))]), lc.mat[,names(lc.mat)[which(xsum.lc>5)]])
lf.mat = cbind(data.frame(PID=m.out$PID[as.numeric(rownames(lf.mat))]), lf.mat[,names(lf.mat)[which(xsum.lf>5)]])
m.outN = cbind(as.data.frame(nut), m.out[-which(names(m.out)=="PID")])
m.outN$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1 = NULL
m.outN$GESUSG6 = NULL
m.outN$EF_LF_Desc_250m = NULL

## assign SoilGrids values to corresponding depths:
## 7 standard depths
breaks = c(0, rowMeans(data.frame(c(0,5,15,30,60,100), c(5,15,30,60,100,200))), 200, 4500)
m.outN$DEPTH_c = cut(m.outN$DEPTH, breaks, labels = paste0("sl", 1:8))
summary(m.outN$DEPTH_c)
#sl1   sl2   sl3   sl4   sl5   sl6   sl7   sl8 
#912 33217 35425 20172 12129 13992  2759   684
hist(m.outN$DEPTH)
## Estimate SG variables per depth:
m.outN$PHIHOX_M = ifelse(m.outN$DEPTH_c=="sl1", m.outN$PHIHOX_M_sl1_250m_ll, ifelse(m.outN$DEPTH_c=="sl2", m.outN$PHIHOX_M_sl2_250m_ll, ifelse(m.outN$DEPTH_c=="sl3", m.outN$PHIHOX_M_sl3_250m_ll, ifelse(m.outN$DEPTH_c=="sl4", m.outN$PHIHOX_M_sl4_250m_ll, ifelse(m.outN$DEPTH_c=="sl5", m.outN$PHIHOX_M_sl5_250m_ll, ifelse(m.outN$DEPTH_c=="sl6", m.outN$PHIHOX_M_sl6_250m_ll, ifelse(m.outN$DEPTH_c=="sl7", m.outN$PHIHOX_M_sl7_250m_ll, NA)))))))
## Bind everything together:
ov = plyr::join_all(list(m.outN[-grep("PHIHOX_M_sl", names(m.outN))], geo.mat, lc.mat, lf.mat), by="PID")
## histogram of depths:
hist(ov$DEPTH[ov$DEPTH<300], breaks=45, xlim=c(0,300), col="gray")

## Regression matrix:
write.csv(ov, file="ov.Nutrients_SoilGrids250m.csv")
unlink("ov.Nutrients_SoilGrids250m.csv.gz")
gzip("ov.Nutrients_SoilGrids250m.csv")
saveRDS.gz(ov, file="regression_matrix_AFNutrients.rds")
save.image()

## ------------- MODEL FITTING -----------

## FIT MODELS:
pr.lst <- names(ov)[-which(names(ov) %in% c(names(as.data.frame(nut)), "DEPTH_c"))]
## 193 vars
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ DEPTH +', paste(pr.lst, collapse="+")))})
names(formulaString.lst) = t.vars
#all.vars(formulaString.lst[[1]])

## Takes ca 45 mins to fit all models (fully parallelized mode)
## Caret training settings (reduce number of combinations to speed up):
Nsub <- 8e3 
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta=c(0.3,0.4), nrounds=c(50,100), max_depth=2:3, gamma=0, colsample_bytree=0.8, min_child_weight=1)
rf.tuneGrid <- expand.grid(mtry=seq(5,55,by=10))
#bm.tuneGrid <- expand.grid(num_trees=c(20,50), k=2, alpha=0.95, beta=2, nu=3)

## Initiate cluster
cl <- makeCluster(56)
registerDoParallel(cl)
## Takes 1 hour to fit all models:
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="AfNutrients_resultsFit.txt")
for(j in t.vars){
  cat("\n", file="AfNutrients_resultsFit.txt", append=TRUE)
  cat(paste("Variable:", j), file="AfNutrients_resultsFit.txt", append=TRUE)
  cat("\n", file="AfNutrients_resultsFit.txt", append=TRUE)
  PID <- ov$PID
  out.rf <- paste0("./models/mrf.",j,".rds")
  if(!file.exists(out.rf)){
    dfs <- ov[,all.vars(formulaString.lst[[j]])]
    sel <- complete.cases(dfs)
    dfs <- dfs[sel,]
    if(nrow(dfs)<Nsub){Nsub=nrow(dfs)}
    ## optimize mtry parameter:
    if(!file.exists(gsub("mrf","t.mrf",out.rf))){
      t.mrfX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
      saveRDS.gz(t.mrfX, file=gsub("mrf","t.mrf",out.rf))
    } else {
      t.mrfX <- readRDS.gz(gsub("mrf","t.mrf",out.rf))
    }
    ## fit RF model using 'ranger' (fully parallelized)
    ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
    mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
    saveRDS.gz(mrfX, file=out.rf)
    ## Top 15 covariates:
    sink(file="AfNutrients_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="AfNutrients_resultsFit.txt", append=TRUE)
    xl <- as.list(ranger::importance(mrfX))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
    ## save fitting success vectors:
    fit.df <- data.frame(PID=PID[sel], observed=dfs[,1], predicted=predictions(mrfX))
    unlink(paste0("RF_fit_", j, ".csv.gz"))
    write.csv(fit.df, paste0("./models/RF_fit_", j, ".csv"))
    gzip(paste0("./models/RF_fit_", j, ".csv"))
    ## XGBoost
    mb.out = paste0("./models/mgb.",j,".rds")
    if(!file.exists(mb.out)){
      ## fit XGBoost model (uses all points):
      mgbX <- caret::train(formulaString.lst[[j]], data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid) 
      saveRDS.gz(mgbX, file=mb.out)
      ## save also binary model for prediction purposes:
      xgb.save(mgbX$finalModel, paste0("./models/Xgb.",j))
    } else {
      mgbX <- readRDS.gz(mb.out)
    }
    importance_matrix <- xgb.importance(mgbX$coefnames, model = mgbX$finalModel)
    cat("\n", file="AfNutrients_resultsFit.txt", append=TRUE)
    print(mgbX)
    cat("\n XGBoost variable importance:\n", file="AfNutrients_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:25,])
    ## bartMachine (needs a lot of RAM and takes at the order of magnitude more time)
    #try( detach("package:bartMachine", unload=TRUE), silent = TRUE)
    #options(java.parameters="-Xmx8g")
    #library(bartMachine)
    #set_bart_machine_num_cores(56)
    #bm.out = paste0("./models/mbm.",j,".rds")
    #if(!file.exists(gsub("mbm","t.mbm",bm.out))){
    #  t.mbMX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="bartMachine", trControl=ctrl, tuneGrid=bm.tuneGrid)
    #  saveRDS.gz(t.mbMX, file=gsub("mbm","t.mbm",bm.out))
    #} else {
    #  t.mbMX <- readRDS.gz(gsub("mbm","t.mbm",bm.out))
    #}
    #mBMX <- train(formulaString.lst[[j]], data=dfs, method="bartMachine", tuneGrid=t.mbMX$finalModel$tuneValue, trControl=trainControl(method="none"))
    #print(summary(mBMX))
    #saveRDS.gz(mBMX, file=bm.out)
    cat("--------------------------------------\n", file="AfNutrients_resultsFit.txt", append=TRUE)
    sink()
  }
}
rm(mrfX); rm(mgbX); #rm(mBMX)
stopCluster(cl); closeAllConnections()

mrfX_lst <- list.files(path="./models", pattern=glob2rx("^mrf.*.rds$"), full.names = TRUE)
mgbX_lst <- list.files(path="./models", pattern=glob2rx("^mgb.*.rds$"), full.names = TRUE)
#mbmX_lst <- list.files(path="./models", pattern=glob2rx("^mBM.*.rds$"), full.names = TRUE)

names(mrfX_lst) <- paste(sapply(basename(mrfX_lst), function(x){ strsplit(strsplit(x, "mrf.")[[1]][2],  ".rds")[[1]][1]}))
names(mgbX_lst) <- paste(sapply(basename(mgbX_lst), function(x){ strsplit(strsplit(x, "mgb.")[[1]][2],  ".rds")[[1]][1]}))
#names(mbmX_lst) <- paste(sapply(basename(mbmX_lst), function(x){ strsplit(strsplit(x, "mBM.")[[1]][2],  ".rds")[[1]][1]}))
save.image()

## ------------- PREPARE PREDICTION LOCATIONS -----------

## Countries:
system(paste0(gdalwarp, ' /data/ESA_global/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif ./stacked250m/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif -ot \"Byte\" -r \"near\" -co \"COMPRESS=DEFLATE\" -t_srs \"', af.prj, '\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner, ' -tr 250 250'))
system(paste0(gdalwarp, ' /data/EcoTapestry/EF_LF_Desc_250m.tif ./stacked250m/EF_LF_Desc_250m.tif -ot \"Byte\" -r \"near\" -co \"COMPRESS=DEFLATE\" -t_srs \"', af.prj, '\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner, ' -tr 250 250'))
system(paste0(gdalwarp, ' /data/LDN/GAUL_COUNTRIES_300m_ll.tif ./stacked250m/GAUL_COUNTRIES.tif -ot \"Int16\" -r \"near\" -co \"COMPRESS=DEFLATE\" -t_srs \"', af.prj, '\" -te ', xllcorner, ' ', yllcorner, ' ', xurcorner, ' ', yurcorner, ' -tr 250 250'))

## Tiling SSA at 250m:
obj <- GDALinfo("./stacked250m/GAUL_COUNTRIES.tif")
tile.lst <- getSpatialTiles(obj, block.x=1e5, return.SpatialPolygons=TRUE)
tile.tbl <- getSpatialTiles(obj, block.x=1e5, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## 5740 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
unlink(c("tiles_100km.shp","ov_ADMIN_tiles.shp"))
writeOGR(tile.pol, "tiles_100km.shp", "tiles_100km", "ESRI Shapefile")
system(paste('gdal_translate ./stacked250m/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif LC2010_250m.sdat -of \"SAGA\" -ot \"Int16\"'))
system(paste0(saga_cmd, ' -c=56 shapes_grid 2 -GRIDS=\"LC2010_250m.sgrd\" -POLYGONS=\"tiles_100km.shp\" -PARALLELIZED=1 -RESULT=\"ov_ADMIN_tiles.shp\"'))
ov_ADMIN = readOGR("ov_ADMIN_tiles.shp", "ov_ADMIN_tiles")
summary(sel.t <- !ov_ADMIN$LC2010_250m.5==210)
## 3227 tiles with values
ov_ADMIN = ov_ADMIN[sel.t,]
plot(ov_ADMIN)
writeOGR(ov_ADMIN, "tiles_100km_sel.shp", "tiles_100km_sel", "ESRI Shapefile")

t.sel = as.character(ov_ADMIN$ID)
new.dirs <- paste0("/data/AFMicroNutrients/tiled/T", t.sel)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

## Some covariates do not match exactly the same grid and need to be stacked:
r250m = raster("./stacked250m/GAUL_COUNTRIES.tif")
r.sel = covs250m.lst[unlist(sapply(names(m.out), function(x){grep(pattern=x, gsub("-", "_", covs250m.lst))}))]
#r.sel = c(r.sel, "/data/ESA_global/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif", "/data/MCD43A4/desertPR_sin.tif")

## TAKES CA 30 mins
sfInit(parallel=TRUE, cpus=38)
sfLibrary(raster)
sfExport("gdalwarp", "gdalwarp_250m", "r250m", "r.sel", "cellsize", "af.prj")
out <- sfClusterApplyLB(r.sel, gdalwarp_250m)
sfStop()

## Countries of interest:
GAUL_names = read.csv("Africa_soil_nutrients_GAUL_names.csv")
sel.country = GAUL_names$Value[GAUL_names$SSA_select==1]
## 47 countries (50 ADMIN areas)
## Existing tifs:
covs.tif = list.files(path="./stacked250m", pattern = ".tif", full.names = TRUE)
covs.tif = covs.tif[-grep("GAUL_COUNTRIES.tif", covs.tif)]
## levels
GESUSG6.lev = levels(m.out$GESUSG6)
LC.lev = levels(m.out$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1)
LF.lev = levels(m.out$EF_LF_Desc_250m)
ov.quant = lapply(ov[,pr.lst], function(i){quantile(i, probs=.5, na.rm=TRUE)})
names(ov.quant) = pr.lst

## ----------- Total clean-up  -----------
#del.lst = list.files(path="/data/AFMicroNutrients/tiled", pattern = ".tif", recursive = TRUE, full.names = TRUE)
#unlink(del.lst)
#del.lst = list.files(path="/data/AFMicroNutrients/tiled", pattern = ".rds", recursive = TRUE, full.names = TRUE)
#unlink(del.lst)
countries = read.csv("/data/LDN/countries/GAUL_names.csv")
## overlay countries and tiles (3-4 minutes):
detach("package:snowfall", unload=TRUE)
detach("package:snow", unload=TRUE)
cl <- makeCluster(56, type="FORK")
ov_ADMIN.lst <- parLapply(cl, as.numeric(t.sel), function(x){ try( summary_GAUL_tiles(x, tile.tbl, countries=countries) )})
stopCluster(cl)
closeAllConnections()
ov_ADMIN.tbl = do.call(rbind, ov_ADMIN.lst)
saveRDS.gz(ov_ADMIN.tbl, file="summary_Africa_GAUL_tiles.rds")
str(ov_ADMIN.tbl) ## 4153

## Clean up some countries:
#del.lst = as.vector(ov_ADMIN.tbl[ov_ADMIN.tbl$Value %in% c("Democratic Republic of the Congo", "Côte d'Ivoire", "South Africa", "United Republic of Tanzania"),"ID"])
#del.lst = as.vector(ov_ADMIN.tbl[ov_ADMIN.tbl$Value %in% c("Abyei", "Hala'ib triangle", "Ilemi triangle"),"ID"])
#x = lapply(del.lst, function(x){file.remove(list.files(path=paste0("/data/AFMicroNutrients/tiled/T", x), pattern=glob2rx(paste0("*_M_sl*_T", x, ".tif")), full.names = TRUE))})
#x = lapply(del.lst, function(x){file.remove(paste0("/data/AFMicroNutrients/tiled/T", x, "/T", x, ".rds"))})
#x = lapply(del.lst, function(x){file.remove(list.files(path=paste0("/data/AFMicroNutrients/tiled/T", x), pattern=glob2rx(paste0("*_M_agg30cm_T", x, ".tif")), full.names = TRUE))})

## ----------- Function to make predictions locations -----------
## Test:
#make_newdata(i=4594, tile.tbl, in.path="/data/AFMicroNutrients/stacked250m", out.path="/data/AFMicroNutrients/tiled", pr.lst=pr.lst, covs.tif=covs.tif, mask="/data/AFMicroNutrients/stacked250m/GAUL_COUNTRIES.tif", sel.country=sel.country, ov.quant=ov.quant, GESUSG6.lev=GESUSG6.lev, LF.lev=LF.lev, LC.lev=LC.lev)

## run in parallel (TAKES ca. 1 hr):
library(snowfall)
sfInit(parallel=TRUE, cpus=56)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(tools)
sfExport("make_newdata", "t.sel", "tile.tbl", "sel.country", "ov.quant", "pr.lst", "covs.tif", "GESUSG6.lev", "LC.lev", "LF.lev")
out <- sfClusterApplyLB(t.sel, function(i){ make_newdata(i, tile.tbl, in.path="/data/AFMicroNutrients/stacked250m", out.path="/data/AFMicroNutrients/tiled", pr.lst=pr.lst, covs.tif=covs.tif, mask="/data/AFMicroNutrients/stacked250m/GAUL_COUNTRIES.tif", sel.country=sel.country, ov.quant=ov.quant, GESUSG6.lev=GESUSG6.lev, LF.lev=LF.lev, LC.lev=LC.lev) })
sfStop()

## ------------- PREDICTIONS -----------

#del.lst = list.files(path="/data/AFMicroNutrients/tiled", pattern=glob2rx("*_T*_rf.rds"), full.names = TRUE, recursive = TRUE)
#unlink(del.lst)

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="./tiled", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)
## 2634 dirs
## Clean-up empty directories
pr.dirs.c <- list.dirs("./tiled")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})

## limit predictions to the succesful models:
sel.vars = c("P","P.B","P.T","N","K","Ca","Mg","Na","B","Fe","Mn","Zn","Cu","Al")
z.min <- as.list(rep(0, length(sel.vars)))
names(z.min) = sel.vars
z.max <- as.list(rep(6e5, length(sel.vars)))
names(z.max) = sel.vars

## Predictions:
type.lst <- as.list(rep("Int16", length(sel.vars)))
names(type.lst) = sel.vars
type.lst[[which(sel.vars=="Ca")]] = "Int32"
type.lst[[which(sel.vars=="Na")]] = "Int32"
mvFlag.lst <- as.list(rep(-32768, length(sel.vars)))
names(mvFlag.lst) = sel.vars

library(ranger)
library(xgboost)
library(tools)
library(parallel)
library(doParallel)
library(rgdal)
library(plyr)
## Run per property (RF = 20 tiles per minute)
for(j in sel.vars){
  try( detach("package:snowfall", unload=TRUE), silent=TRUE)
  try( detach("package:snow", unload=TRUE), silent=TRUE)
  if(j %in% c("P","P.B","S","B","Cu","Zn")){ multiplier = 100 }
  if(j %in% c("C","P.T","N","K","Ca","Mg","Na","Fe","Mn","Al")){ multiplier = 1 }
  ## Random forest predictions:
  gm = readRDS.gz(paste0("./models/mrf.", j,".rds"))
  gm1.w = 1/gm$prediction.error
  ## Estimate amount of RAM needed per core
  cpus = unclass(round((256-50)/(3.5*(object.size(gm)/1e9))))
  cl <- makeCluster(ifelse(cpus>46, 46, cpus), type="FORK")
  x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(any(!file.exists(paste0("/data/AFMicroNutrients/tiled/", x, "/", j, "_M_sl", 1:5, "_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="/data/AFMicroNutrients/tiled", out.path="/data/AFMicroNutrients/tiled", split_no=NULL, varn=j, method="ranger", sd=c(0, 5, 15, 30, 60), DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/AFMicroNutrients/tiled/", x, "/", x,".rds"), SG.col=paste0("PHIHOX_M_sl", 1:5, "_250m_ll"), SG.col.name="PHIHOX_M" ) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## XGBoost:
  gm = readRDS.gz(paste0("./models/mgb.", j,".rds"))
  gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
  cpus = unclass(round((256-30)/(3.5*(object.size(gm)/1e9))))
  cl <- makeCluster(ifelse(cpus>46, 46, cpus), type="FORK")
  x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(any(!file.exists(paste0("/data/AFMicroNutrients/tiled/", x, "/", j, "_M_sl", 1:5, "_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="/data/AFMicroNutrients/tiled", out.path="/data/AFMicroNutrients/tiled", split_no=NULL, varn=j, method="xgboost", sd=c(0, 5, 15, 30, 60), DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/AFMicroNutrients/tiled/", x, "/", x,".rds"), SG.col=paste0("PHIHOX_M_sl", 1:5, "_250m_ll"), SG.col.name="PHIHOX_M" ) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## sum up predictions:
  library(snowfall)
  sfInit(parallel=TRUE, cpus=46)
  sfExport("t.sel", "sum_predict_ensemble", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(paste0("T", t.sel), fun=function(x){ try( sum_predict_ensemble(x, in.path="/data/AFMicroNutrients/tiled", out.path="/data/AFMicroNutrients/tiled", varn=j, num_splits=NULL, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]], rds.file=paste0("/data/AFMicroNutrients/tiled/", x, "/", x,".rds")) ) } )
  sfStop()
  gc(); gc()
}

## clean-up:
#for(i in sel.vars){ ## c("N","K","Ca","Mg","B","Fe","Mn","Zn","Cu") 
# del.lst <- list.files(path="/data/AFMicroNutrients/tiled", pattern=glob2rx(paste0("^", i, "_M_sl?_*.tif$")), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)
#}

## Aggregate values to get 0-30 cm estimates:
## test it:
#agg_layers(i="T3713", varn="Al", ot=type.lst["Al"], dstnodata=mvFlag.lst["Al"])
for(j in sel.vars){
  sfInit(parallel=TRUE, cpus=46)
  sfExport("pr.dirs", "agg_layers", "sel.vars", "type.lst", "mvFlag.lst", "j")
  sfLibrary(rgdal)
  sfLibrary(raster)
  x <- sfClusterApplyLB(pr.dirs, function(x){try( agg_layers(i=x, varn=j, ot=type.lst[j], dstnodata=mvFlag.lst[j]) )})
  sfStop()
}

## Organic carbon, BLD and CRF from SoilGrids250m:
for(k in c("ORCDRC","BLDFIE","CRFVOL")){
  sfInit(parallel=TRUE, cpus=56)
  sfExport("t.sel", "tile_SoilGrids", "gdalwarp", "af.prj", "tile.tbl", "k")
  sfLibrary(rgdal)
  x <- sfClusterApplyLB(paste(t.sel), function(x){try( tile_SoilGrids(i=x, tile.tbl=tile.tbl, varn=k, in.tifs = paste0("/data/GEOG/", k, "_M_sl", 1:4, "_250m_ll.tif"), out.path="/data/AFMicroNutrients/tiled", t_srs=af.prj) ) })
  sfStop()
}

## ----------- Aggregate SoilGrids -----------
for(k in c("ORCDRC","BLDFIE","CRFVOL")){
  sfInit(parallel=TRUE, cpus=56)
  sfExport("pr.dirs", "agg_layers", "k")
  sfLibrary(rgdal)
  sfLibrary(raster)
  x <- sfClusterApplyLB(pr.dirs, function(x){try( agg_layers(i=x, varn=k, d=c(0,5,15,30)) )})
  sfStop()
}

## Create mosaicks:
sfInit(parallel=TRUE, cpus=ifelse(length(sel.vars)>45, 45, length(sel.vars)))
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "sel.vars", "mvFlag.lst", "make_mosaick", "type.lst")
out <- sfClusterApplyLB(1:length(sel.vars), function(x){ try( make_mosaick("M_agg30cm", varn=sel.vars[x], in.path="/data/AFMicroNutrients/tiled", ot=type.lst[x], dstnodata=mvFlag.lst[x], tile.name="AF") )})
sfStop()

sel.sg = c("ORCDRC","BLDFIE","CRFVOL")
sfInit(parallel=TRUE, cpus=3)
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "make_mosaick", "sel.sg")
out <- sfClusterApplyLB(1:length(sel.sg), function(x){ try( make_mosaick("M_agg30cm", varn=sel.sg[x], in.path="/data/AFMicroNutrients/tiled", tile.name="AF") )})
sfStop()

## ----- Cluster analysis ----

library(compositions)
library(leaflet)
library(ranger)

tMV.vars <- c("C","N","P","P.T","K","S","Ca","Mg","Na","B","Cu","Fe","Mn","Zn","Al")
#load("nut.rda")
#nut$P.a <- rowMeans(nut@data[,c("P","P.O")], na.rm = TRUE)
#nut$P.a <- ifelse(is.nan(nut$P.a), NA, nut$P.a)
#str(nut@data)

## Fig_AfNutrients_histograms_macro.pdf
## Macro nutrients
tMV.varsH <- c("C","N","K","P","P.T","Ca","Mg","Na","S")
tMV.varsHnames <- c("org. C", "org. N", "ext. K", "ext. P","tot. P","ext. Ca","ext. Mg","ext. Na", "ext. S")
df = data.frame(ppm=unlist(nut@data[,tMV.varsH]), var=as.vector(sapply(tMV.varsHnames, function(x){rep(x, nrow(nut@data))})))
df = df[!is.na(df$ppm)&df$ppm<5e4,]
#str(df)
hist(log1p(df$ppm))
#brks = seq(range(log1p(df$ppm))[1], range(log1p(df$ppm))[2], length.out = 25)
#df$values <- cut(df$ppm, breaks=brks, labels = round(expm1(brks))[-1])
p <- ggplot(data = df, aes(x=ppm))
p <- p + scale_x_log10()
p <- p + geom_histogram(aes(fill=var))
p <- p + scale_fill_brewer(palette="Set3")
p <- p + facet_wrap( ~ var, ncol=1)
p

## Fig_AfNutrients_histograms_micro.pdf
## Micro nutrients:
tMV.varsL <- c("Fe","Mn","Zn","Cu","B")
tMV.varsLnames <- c("ext. Fe", "ext. Mn", "ext. Zn", "ext. Cu","ext. B")
df2 = data.frame(ppm=unlist(nut@data[,tMV.varsL]), var=as.vector(sapply(tMV.varsLnames, function(x){rep(x, nrow(nut@data))})))
df2 = df2[!is.na(df2$ppm)&df2$ppm<4e4,]
p2 <- ggplot(data = df2, aes(x=ppm))
p2 <- p2 + scale_x_log10()
p2 <- p2 + geom_histogram(aes(fill=var))
p2 <- p2 + scale_fill_brewer(palette="Set2")
p2 <- p2 + facet_wrap( ~ var, ncol=1)
p2

## Many missing values that need to be predicted
#load("ov.rda")
tMV.nuts = c("C","N","P","P.O","P.T","K","S","Ca","Mg","Na","B","Cu","Fe","Mn","Zn","Al")
nutMV.df = nut@data
## Replace missing values with predictions (this is tricky)
nutMV.df = nutMV.df[apply(nutMV.df[,tMV.nuts], 1, function(x){sum(is.na(x))})<length(tMV.nuts),]
for(i in tMV.nuts){
  gm = readRDS(paste0("/data/AFMicroNutrients/models/mrf.",i,".rds"))
  ## missing values
  sel = nutMV.df$PID[which(is.na(nutMV.df[,i]))]
  pr.df = ov[match(sel, ov$PID),]
  pr.df = pr.df[complete.cases(pr.df[,gm$forest$independent.variable.names]),]
  v = predict(gm, pr.df, na.action=na.pass)
  nutMV.df[match(pr.df$PID, nutMV.df$PID),i] = ifelse(v$predictions<0, 0, v$predictions)
}

#nutMV.df$P.a <- (0.8*nutMV.df$P + 0.2*nutMV.df$P.O)/2
nutMV.df$C = ifelse(nutMV.df$C>expm1(13), NA, nutMV.df$C)
sel2 = complete.cases(nutMV.df[,c(tMV.varsH,tMV.varsL)])
dfc0 = nutMV.df[sel2,c(tMV.varsH,tMV.varsL)]
## 110,668 observations with all values
## convert to compositions:
dfc0.c = acomp(dfc0, total=1e6)
sum(is.na(dfc0.c))
sum(is.infinite(dfc0.c))
#boxplot(dfc0.c[1:500,])
pc.nut <- princomp.acomp(dfc0.c) ## , center = TRUE, scale. = TRUE
plot(pc.nut)
#straight(mean(dfc0.c), pc.nut$Loadings)
fviz_pca_var(pc.nut, title=NULL, labelsize=6) ## Fig_prcomp_1_2_biplot.pdf
fviz_pca_var(pc.nut, axes=3:4, title=NULL, labelsize=6) ## Fig_prcomp_3_4_biplot.pdf

## Optimal number of clusters (http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters):
library(fpc)
library(cluster)
#pamk.best <- pamk(dfc0.c[sample.int(nrow(dfc0.c),6e3),], krange=4:25)
#pamk.best$nc
## 6 clusters -> probably too few?
#plot(pam(dfc0.c[sample.int(nrow(dfc0.c),1e3),], pamk.best$nc))
library(mclust)
# Run the function to see how many clusters (TAKES 10 MINS)
d_clust <- Mclust(dfc0.c[sample.int(nrow(dfc0.c),1.5e4),], G=4:28)
m.best <- dim(d_clust$z)[2]
m.best
# 20 clusters
#plot(d_clust)

## Fuzzy k-means clustering:
fkm.df = plyr::join(nutMV.df, ov, by="PID")
library(h2o)
h2o.init(nthreads = -1)
dfc0.c.hex <- as.h2o(dfc0.c, destination_frame = "dfc0.c")
km.nut <- h2o.kmeans(training_frame=dfc0.c.hex, k=m.best, keep_cross_validation_predictions = TRUE)
m.km <- as.data.frame(h2o.predict(km.nut, dfc0.c.hex, na.action=na.pass))
## Class centers:
class_dfc0.c = as.data.frame(h2o.centers(km.nut))
str(class_dfc0.c)
write.csv(class_dfc0.c, "NCluster_20_class_centers.csv")
## back-transform acomp?

fkm.df$Clusters = NA
fkm.df[sel2,"Clusters"] = m.km$predict
fkm.df$Clusters = as.factor(fkm.df$Clusters)
summary(fkm.df$Clusters)
## largest classes: "17", "21"and "0"
saveRDS(fkm.df, "ov_mcluster.rds")
class_dfc0.cL = data.frame(cluster=paste0("c", 1:m.best))
for(i in names(dfc0)){
  fm = as.formula(paste(i, "~ Clusters-1"))
  class_dfc0.cL[,i] = signif(coef(lm(fm, fkm.df)), 3)
}
write.csv(class_dfc0.cL, "NCluster_20_class_centers.csv")
h2o.shutdown()

## fit random forest model:
fm.class = as.formula(paste("Clusters ~ ", paste(gm$forest$independent.variable.names[-1], collapse="+")))
fkm.df.s = fkm.df[complete.cases(fkm.df[,all.vars(fm.class)]),all.vars(fm.class)]
mrf.class = ranger(fm.class, fkm.df.s, importance = "impurity", probability = TRUE, write.forest = TRUE, num.trees=120)
mrf.class
## 70% accuracy
saveRDS(mrf.class, "mrf_mcluster.rds")


## ----------- Predict Mclusters -----------
#mrf.class = readRDS("mrf_mcluster.rds")

library(ranger)
library(tools)
library(rgdal)
library(plyr)
library(parallel)
library(doParallel)
cl <- makeCluster(46, type="FORK")
x = parLapply(cl, paste0("T", t.sel), fun=function(x){ try( predict_mcluster_tiles(i=x, gm=mrf.class, tile.tbl=tile.tbl, in.path="/data/AFMicroNutrients/tiled", out.path="/data/AFMicroNutrients/tiled", nlevs=ncol(mrf.class$predictions)) ) } )
stopCluster(cl)

library(snowfall)
cl.lst = c("M", 1:ncol(mrf.class$predictions))
sfInit(parallel=TRUE, cpus=ifelse(length(cl.lst)>45, 45, length(cl.lst)))
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "cl.lst", "make_mosaick")
out <- sfClusterApplyLB(1:length(cl.lst), function(x){ try( make_mosaick(cl.lst[x], varn="NCluster", in.path="/data/AFMicroNutrients/tiled", ot="Byte", dstnodata="255", tile.name="AF") )})
sfStop()

## ----------- derive Scaled Shannon Entropy -----------
##(100 is a maximum error; 0 is perfect prediction)

varn = "NCluster"
levs = paste0(1:ncol(mrf.class$predictions))
sfInit(parallel=TRUE, cpus=56)
sfExport("entropy_tile", "levs", "varn")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(pr.dirs, function(i){try( entropy_tile(i, in.path="/data/AFMicroNutrients/tiled", varn, levs) )})
sfStop()

tmp.lst <- list.files(path="/data/AFMicroNutrients/tiled", pattern=glob2rx(paste0("SSI_NCluster_*.tif$")), full.names=TRUE, recursive=TRUE)
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/Nutrients/SSI_NCluster_AF_250m.tif -ot \"Byte\" -dstnodata \"255\" -r \"near\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/Nutrients/SSI_NCluster_AF_1km.tif -ot \"Byte\" -dstnodata \"255\" -tr 1000 1000 -r \"average\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' ', vrt.tmp, ' SSI_NCluster_AF_5km.tif -ot \"Byte\" -dstnodata \"255\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.05 0.05 -r \"average\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -overwrite'))

## 5 km maps (final output)
tif.lst <- list.files(path="/data/GEOG/Nutrients", pattern="1km.tif", full.names = TRUE)
sfInit(parallel=TRUE, cpus=length(tif.lst))
sfExport("gdalwarp", "tif.lst")
out <- sfClusterApplyLB(tif.lst, function(x){ try( system(paste0(gdalwarp, ' ', x, ' ', basename(gsub("1km", "5km", x)), ' -r \"average\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.05 0.05 -co \"COMPRESS=DEFLATE\" -overwrite')) ) })
sfStop()

## NCluster:
as.vector(extent(raster("N_M_agg30cm_AF_5km.tif")))[c(1,3,2,4)]
system(paste0(gdalwarp, ' /data/GEOG/Nutrients/NCluster_M_AF_250m.tif NCluster_M_AF_5km.tif  -r \"near\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.05 0.05 -co \"COMPRESS=DEFLATE\" -overwrite -te ', paste(as.vector(extent(raster("N_M_agg30cm_AF_5km.tif")))[c(1,3,2,4)], collapse = " ")))
## USDA classes:
system(paste0(gdalwarp, ' /data/GEOG/TAXOUSDA_Usterts_250m_ll.tif TAXOUSDA_Usterts_5km_ll.tif -r \"average\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.05 0.05 -co \"COMPRESS=DEFLATE\" -overwrite -te ', paste(as.vector(extent(raster("N_M_agg30cm_AF_5km.tif")))[c(1,3,2,4)], collapse = " ")))
system(paste0(gdalwarp, ' /data/GEOG/TAXOUSDA_Humults_250m_ll.tif TAXOUSDA_Humults_5km_ll.tif -r \"average\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.05 0.05 -co \"COMPRESS=DEFLATE\" -overwrite -te ', paste(as.vector(extent(raster("N_M_agg30cm_AF_5km.tif")))[c(1,3,2,4)], collapse = " ")))

#unlink(paste0("CV_", sel.vars, ".rda"))
## Cross-validation 10-fold (TH: this does not takes into account high spatial clustering):
source("/data/cv/cv_functions.R")
cat("Results of Cross-validation:\n\n", file="AfNutrients_resultsCV_5fold.txt")
cv_lst <- rep(list(NULL), length(sel.vars))
for(j in 1:length(sel.vars)){
  if(!file.exists(paste0("CV_", sel.vars[j], ".rda"))){
    cat(paste("Variable:", sel.vars[j]), file="AfNutrients_resultsCV_5fold.txt", append=TRUE)
    cat("\n", file="AfNutrients_resultsCV_5fold.txt", append=TRUE)
    rmatrix = ov[complete.cases(ov[,all.vars(formulaString.lst[[sel.vars[j]]])]),c("PID", all.vars(formulaString.lst[[sel.vars[j]]]))]
    cv_lst[[j]] <- cv_numeric(formulaString=formulaString.lst[[sel.vars[j]]], rmatrix=rmatrix, nfold=5, idcol="PID", Log=TRUE)
    sink(file="AfNutrients_resultsCV_5fold.txt", append=TRUE, type="output")
    print(cv_lst[[j]]$Summary)
    cat("\n", file="AfNutrients_resultsCV_5fold.txt", append=TRUE)
    sink()
    assign(paste0("CV_", sel.vars[j]), cv_lst[[j]])
    save(list=paste0("CV_", sel.vars[j]), file=paste0("CV_", sel.vars[j], ".rda"))
  }
}

## ----------- correlation plots -----------
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
source("/data/models/plot_hexbin.R")
plt.names <- c("ext. K", "ext. Mg", "ext. Ca", "ext. Na", "ext. Zn", "ext. B", "ext. Cu", "ext. Fe", "ext. Mn", "ext. Al", "org. N", "tot. P", "ext. P") ## "ext. S", "ext. P"
names(plt.names) = c("K", "Mg", "Ca", "Na", "Zn", "B", "Cu", "Fe", "Mn", "Al", "N", "P.T", "P") # "S"
breaks.lst <- list(
  c(0,5,10,seq(20,1500,length=47)), # K
  c(0,seq(5,2500,length=49)), # Mg
  c(0,seq(10,14500,length=49)), # Ca
  c(0,seq(10,2700,length=49)), # Na
  c(0,seq(0.5,26,length=49)), # Zn
  c(0,seq(0.2,2.1,length=49)), # B
  c(0,seq(0.5,11,length=49)), # Cu
  c(0,seq(2,580,length=49)), # Fe
  c(0,seq(1,480,length=49)), # Mn
  c(0,seq(5,2200,length=49)), # Al
  c(0,seq(20,4300,length=49)), # N
  c(0,seq(20,3100,length=49)), # P.T
  c(0,seq(5,440,length=49)) # P
) ## c(0,seq(0.5,190,length=49)),
names(breaks.lst) = names(plt.names)
plt.log <- rep(TRUE, length(plt.names))
names(plt.log) = names(plt.names)

for(j in 1:length(plt.names)){
  plot_hexbin(varn=names(plt.names)[j], breaks=breaks.lst[[j]], main=plt.names[j], in.file=paste0("CV_", names(plt.names)[j], ".rda"), log.plot=plt.log[j])
}

## Fig_results_CV.png
system(paste('montage -mode concatenate -tile 3x5 ', paste(paste0("/data/AFMicroNutrients/plot_CV_", c("N","P.T","K","Ca","Mg","Na","S","Al","P","B","Cu","Fe","Mn","Zn"),".png"), collapse=" "),' Fig_results_CV.png'))

## ----------- OFRA data  -----------
## available for download [http://ec2-54-93-187-255.eu-central-1.compute.amazonaws.com/]
GrainYield.OFRA <- readRDS("GrainYield.OFRA.rds")
summary(GrainYield.OFRA$grainyield)
## clean-up:
dim(GrainYield.OFRA)
levels(GrainYield.OFRA$crop1)
GrainYield.OFRA$crop1 = as.factor(ifelse(GrainYield.OFRA$crop1=="maize", "Maize", ifelse(GrainYield.OFRA$crop1=="Millet", "Pearl Millet", paste(GrainYield.OFRA$crop1))))
summary(GrainYield.OFRA$crop1)
summary(GrainYield.OFRA$variety1)
GrainYield.OFRA$variety1 = as.factor(ifelse(GrainYield.OFRA$variety1=="", "Unknown", paste(GrainYield.OFRA$variety1)))
summary(GrainYield.OFRA$applimethod)
GrainYield.OFRA$applimethod = as.factor(ifelse(GrainYield.OFRA$applimethod=="", "Unknown", paste(GrainYield.OFRA$applimethod)))
plot(admin.af)
points(GrainYield.OFRA)
crop1.mat = data.frame(model.matrix(~crop1-1, GrainYield.OFRA@data))
## Remove classes with too little observations
variety1.mat = data.frame(model.matrix(~variety1-1, GrainYield.OFRA@data))
xsum.variety1 = sapply(variety1.mat, sum, na.rm=TRUE)
sel.variety1 = names(variety1.mat)[which(xsum.variety1>20)]
applimethod.mat = data.frame(model.matrix(~applimethod-1, GrainYield.OFRA@data))
xsum.applimethod = sapply(applimethod.mat, sum, na.rm=TRUE)
sel.applimethod = names(applimethod.mat)[which(xsum.applimethod>20)]

## overlay with nutrients and climatic vars
Ncovs250m.lst = c(list.files(path="/data/Climate", pattern=glob2rx("CHELSA_prec_*.tif$"), full.names = TRUE), list.files(path="/data/Climate", pattern=glob2rx("CHELSA_temp_*.tif$"), full.names = TRUE), list.files(path="/data/GEOG/Nutrients", pattern=glob2rx("*_M_agg30cm_AF_250m.tif$"), full.names = TRUE))
#TAKES 20 MINS:
sfInit(parallel=TRUE, cpus=56)
sfExport("GrainYield.OFRA", "Ncovs250m.lst", "extract.tif")
sfLibrary(rgdal)
sfLibrary(raster)
Nm.out <- data.frame(sfClusterApplyLB(Ncovs250m.lst, function(i){try( extract.tif(i, GrainYield.OFRA["SOURCEID"]) )}))
sfStop()
names(Nm.out) = file_path_sans_ext(basename(Ncovs250m.lst))
names(Nm.out) = gsub("-", "_", names(Nm.out))

## Regression matrix:
Nov = do.call(cbind, list(GrainYield.OFRA@data[,c("SOURCEID","grainyield","crop1")], Nm.out, crop1.mat, variety1.mat[,sel.variety1], applimethod.mat[,sel.applimethod]))
summary(Nov$applimethodUnknown==1)
write.csv(Nov, file="ov.OFRAyields_AFNutrients.csv")
unlink("ov.OFRAyields_AFNutrients.csv.gz")
gzip("ov.OFRAyields_AFNutrients.csv")
saveRDS.gz(Nov, file="regression_matrix_OFRAyields_AFNutrients.rds")

## Yield as a function of "crop1", "variety1", "applimethod", nutrients and climatic data
N.fm = as.formula(paste(' grainyield ~ ', paste(names(Nm.out), collapse="+"), "+", paste(names(crop1.mat), collapse = "+"), "+", paste(sel.variety1, collapse = "+"), "+", paste(sel.applimethod, collapse = "+")))
m.grainyield.rfX <- ranger(N.fm, data=Nov, importance="impurity", write.forest=TRUE, num.trees=85)
m.grainyield.rfX
## 78%
xl <- as.list(ranger::importance(m.grainyield.rfX))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:35]])))
xl = t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:20]]))
save.image()

pdf(file = "Fig_variable_importance_grainyield.pdf", width = 7, height = 7.5)
par(mar=c(2.5,10,2.5,0.5), oma=c(1,1,1,1))
plot(x=rev(xl)/max(xl)*100, y=1:20, pch = 19, col="blue", xlab="Importance (%)", xlim=c(0,105), ylim=c(0,21), yaxp=c(0,20,20), xaxs="i", yaxs="i", cex=1.4, yaxt="n", ylab="", main="Model importance plot (grain yield)", cex.main=1)
abline(h=1:20, lty=2, col="grey")
#axis(2, at=1:20, labels=rev(attr(xl, "dimnames")[[1]]), las=2)
axis(2, at=1:20, labels=rev(c("Crop type - Cassava", "Crop type - Maize", "Crop type - Potato", expression(bold("ext. Mn")), expression(bold("ext. Zn")), "Variety - B53", expression(bold("ext. Cu")), "Precipitation June", "Precipitation October", "Precipitation total", "Precipitation September", "Variety - Asante", "Precipitation May", expression(bold("ext. Al")), "Precipitation July", "Temperature November", "Precipitation December", "Precipitation April", "App. - fertilizer applied", "Precipitation March")), las=2)
dev.off()

## Stack to 1 km resolution:
r1km = raster("/data/GEOG/Nutrients/Fe_M_agg30cm_AF_1km.tif")
r1km.sel = c(Ncovs250m.lst[1:26], paste0("/data/GlobCover30/L0", 1:9, "GLC3a.tif"))

sfInit(parallel=TRUE, cpus=length(r1km.sel))
sfLibrary(raster)
sfExport("gdalwarp", "r1km.sel", "gdalwarp_1km", "r1km")
out <- sfClusterApplyLB(r1km.sel, gdalwarp_1km, r1km=r1km)
sfStop()
## Error: Missing data in columns: CHELSA_prec_1979_2013_land, CHELSA_temp_1979_2013_land.
system('gdalinfo /data/Climate/CHELSA_prec_1_1979-2013.tif')
system(paste0('gdalwarp /data/Climate/CHELSA_prec_1979-2013_land.tif ./stacked1km/CHELSA_prec_1979-2013_land.tif -co \"COMPRESS=DEFLATE\" -tr 1000 1000 -co \"BIGTIFF=YES\" -multi -wm 2000 -overwrite -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"', proj4string(r1km), '\" -te ', paste0(extent(r1km)[c(1,3,2,4)], collapse=" ")))
system(paste0('gdalwarp /data/Climate/CHELSA_temp_1979-2013_land.tif ./stacked1km/CHELSA_temp_1979-2013_land.tif -co \"COMPRESS=DEFLATE\" -tr 1000 1000 -co \"BIGTIFF=YES\" -multi -wm 2000 -overwrite -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"', proj4string(r1km), '\" -te ', paste0(extent(r1km)[c(1,3,2,4)], collapse=" ")))

## Agricultural land mask:
grid1km = readGDAL("/data/GEOG/Nutrients/Fe_M_agg30cm_AF_1km.tif")
names(grid1km) = "Fe_M_agg30cm_AF_250m"
## 23 million pixels
grid1km = as(grid1km, "SpatialPixelsDataFrame")
glc.lst = list(NULL)
for(i in paste0("./stacked1km/L0", c(1,3,4,8), "GLC3a.tif")){
  glc.lst[[basename(i)]] = readGDAL(i)$band1[grid1km@grid.index]
}
sel.pix1km = glc.lst[["L01GLC3a.tif"]] > 0 | glc.lst[["L03GLC3a.tif"]] > 0 | glc.lst[["L04GLC3a.tif"]] > 0 & glc.lst[["L08GLC3a.tif"]] < 80 
summary(sel.pix1km)
## 17 million pixels
rm(glc.lst)
grid1km.s = grid1km[sel.pix1km,]
#plot(raster(grid1km.s), col=SAGA_pal[[1]])

## ----------- Predict yield using crop1, climate and nutrients  -----------
pred.1km = c(list.files(path="./stacked1km", pattern=glob2rx("CHELSA_*.tif$"), full.names = TRUE),  list.files(path="/data/GEOG/Nutrients", pattern=glob2rx("*_M_agg30cm_AF_1km.tif$"), full.names = TRUE))
for(i in 1:length(pred.1km)){
  grid1km.s@data[,file_path_sans_ext(basename(pred.1km[i]))] = readGDAL(pred.1km[i])$band1[grid1km.s@grid.index]
}
names(grid1km.s) = gsub("-", "_", names(grid1km.s))
names(grid1km.s) = gsub("1km", "250m", names(grid1km.s))
summary(grid1km.s$CHELSA_prec_1979_2013_land)
summary(grid1km.s$CHELSA_temp_1979_2013_land)
summary(grid1km.s$ORCDRC_M_agg30cm_AF_250m)
## Add all dummy vars:
for(i in c(paste(names(crop1.mat)), paste(sel.variety1), paste(sel.applimethod))){
  grid1km.s@data[,i] = 0
}
## filter out missing values (takes few minutes):
for(i in names(grid1km.s)){
  if(sum(is.na(grid1km.s@data[,i]))>0){
    grid1km.s@data[,i] = ifelse(is.na(grid1km.s@data[,i]), quantile(grid1km.s@data[,i], 0.5, na.rm=TRUE), grid1km.s@data[,i])
  }
}

## Predict and generate a geoTIFF (takes >30 minutes per crop1!):
grid1km.s@data[,"applimethodUnknown"] = 1
grid1km.s@data[,"variety1Unknown"] = 1
for(k in names(crop1.mat)){
  out.tif = paste0("/data/GEOG/Nutrients/PotYield_", k, "_1km.tif")
  if(!file.exists(out.tif)){
    for(i in names(crop1.mat)){ grid1km.s@data[,i] = 0 }
    grid1km.s@data[,k] = 1
    grid1km.s$pred = round(predict(m.grainyield.rfX, grid1km.s@data, num.threads=56)$predictions * 10)
    writeGDAL(grid1km.s["pred"], out.tif, type="Int16", mvFlag=-32678, options="COMPRESS=DEFLATE")
    gc(); gc()
  }
}
#plot(raster(grid1km.s["pred"]), col=SAGA_pal[[1]])

## ------------- Generate Inspire metadata files for each zipped GEOTIFF -----------

landmask = readGDAL("/data/GEOG/Nutrients/Zn_M_agg30cm_AF_1km.tif")
#landmask <- as(landmask["mask"], "SpatialPixelsDataFrame")

out.tif.lst = c(paste0("/data/GEOG/Nutrients/",c("Ca","Na","N","Zn","Mg","Al","Cu","K","B","Fe","P.T","Mn","P"), "_M_agg30cm_AF_250m.tif"), "/data/GEOG/Nutrients/NCluster_M_AF_250m.tif")
ATTRIBUTE_LABEL = c("Ca","Na","N","Zn","Mg","Al","Cu","K","B","Fe","P.T","Mn","P","Ncluster")
ATTRIBUTE_TITLE = c(paste0("Extractable ", c("Ca","Na","N","Zn","Mg","Al","Cu","K","B","Fe","P.T","Mn","P"), " for 0--30 cm depth in ppm"), "Nutrient clusters based on fuzzy k-means")
ATTRIBUTE_TITLE = gsub("Extractable N", "Organic N", ATTRIBUTE_TITLE)
ATTRIBUTE_TITLE = gsub("Extractable P.T", "Total P", ATTRIBUTE_TITLE)
ATTRIBUTE_TITLE = sapply(c("P","B","Cu","Zn"), function(i){ gsub(paste0(i," for 0--30 cm depth in ppm"), paste0(i," for 0--30 cm depth in 100 ppm"), ATTRIBUTE_TITLE)})
ATTRIBUTE_TITLE

for(i in 1:length(out.tif.lst)){
  out.xml.file = gsub(".tif", ".tif.xml", out.tif.lst[i]) 
  if(!file.exists(out.xml.file)){
    xmd <- spMetadata(landmask, out.xml.file=out.xml.file,
                      md.type="INSPIRE",
                      CI_Citation_title = paste("Africa SoilGrids:", ATTRIBUTE_LABEL[i], ":", ATTRIBUTE_TITLE[i]),
                      CI_Electronic_mail_address = "tom.hengl@isric.org",
                      CI_Organisation_name = "ISRIC / AfSIS / EthioSIS and other",
                      Date_stamp = as.Date("2017-03-14", format="%Y-%m-%d"),
                      CI_Citation_date = as.Date("2016-12-05", format="%Y-%m-%d"),
                      MD_Thesaurus_date = as.Date("2016-12-05", format="%Y-%m-%d"),
                      CI_Unique_name = ATTRIBUTE_LABEL[i],
                      MD_Abstract = ATTRIBUTE_TITLE[i],
                      MD_Organisation_name = "ISRIC - World Soil Information",
                      MD_Electronic_mail_address = "tom.hengl@isric.org",
                      MD_Keyword = c("soil", "nutrients"),
                      MD_Use_limitations = "Open Data Commons Open Database License (ODbL). This means that: You must attribute any public use of the database, or works produced from the database, in the manner specified in the ODbL. For any use or redistribution of the database, or works produced from it, you must make clear to others the license of the database and keep intact any notices on the original database. Share-Alike: If you publicly use any adapted version of this database, or works produced from an adapted database, you must also offer that adapted database under the ODbL. Keep open: If you redistribute SoilGrids geotifs, or an adapted version of it, then you may use technological measures that restrict the work (such as DRM) as long as you also redistribute a version without such measures.",
                      MD_Other_restrictions = "http://www.isric.org/content/disclaimer-soilgrids",
                      MD_Equivalent_scale = "200000",
                      MD_Resolution = 250,
                      Time_period_begin = "1950-01-01", ## as.Date("1950-01-01", format="%Y-%m-%d"),
                      Time_period_end = "2015-12-31",  ## as.Date("2005-12-31", format="%Y-%m-%d")
                      Extent_West_Longitude = -33,
                      Extent_East_Longitude = 59,
                      Extent_South_Latitude = -36,
                      Extent_North_Latitude = 30,
                      DQ_Lineage_statement = paste("Values M = mean value predicted using Machine Learning (ensemble between random forest and gradient boosting). Measurement unit: 250 m. The automated soil mapping system 'SoilGrids' is explained in detail in http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169748")
    )
    metadata = c(ATTRIBUTE_LABEL[i], "ISRIC / AfSIS / EthioSIS / PBL and other", ATTRIBUTE_TITLE[i], "http://geonode.isric.org")
    names(metadata) = c("VARIABLE_NAME", "CITATION_ORIGINATOR", "ATTRIBUTE_TITLE", "DOWNLOAD_FTP_URL")
    m = paste('-mo ', '\"', names(metadata), "=", as.vector(metadata), '\"', sep="", collapse = " ")
    command = paste0('gdal_edit.py ', m,' ', out.tif.lst[i])
    system(command, intern=TRUE)
  }
}

## ---- Final plots ----


ssa.lst = c("Cameroon","Ethiopia","Ghana","Kenya","Mali","Malawi","Nigeria","Rwanda", "Senegal","Tanzania","Uganda","Eritrea","Djibouti","South Sudan","Sudan","Central African Republic","Benin","Togo","Guinea","Guinea-Bissau","Liberia","Sierra Leone","Equatorial Guinea","Gabon","Congo","Republic of Congo","Burundi","Zimbabwe","Mozambique","Namibia","Swaziland","Botswana","Lesotho","Ivory Coast","Gambia","Niger","Chad","Cape Verde","Burkina Faso","Zambia","Madagascar","Comoros","Somalia","Djibouti","Sao Tome and Principe","Angola","South African Republic","Somaliland")
wg.url <- url("http://gsif.isric.org/lib/exe/fetch.php?media=admin.af.rda")
load(wg.url)
proj4string(admin.af) <- "+proj=longlat +datum=WGS84"
country.ssa2 <- admin.af[-c(14,26,28,49,15),]
writeOGR(country.ssa2, "SubSaharan_Africa_admin.shp", "SubSaharan_Africa_admin", "ESRI Shapefile")
country.af <- as(admin.af, "SpatialLines")
country.ssa <- country.af[-c(14,26,28,49,15),]

## Subsaharan Africa 
bb = c(-18, -35, 52, 17.5)
sel.sim = grep(pattern=glob2rx("^Simulated_*"), nut$PID)
nutMV.ssa = nut[-sel.sim,]
nutMV.ssa = nutMV.ssa[nutMV.ssa@coords[,2]<bb[4],]
png(file = "Fig_AfNutrients_spatial_coverage.png", res = 150, width = 1920, height = 1680)
par(mfrow=c(2,2), mar=c(0,0,4,0), oma=c(0,0,0,0))
plot(country.ssa, col="darkgrey", ylim=bb[c(2,4)], xlim=bb[c(1,3)], main="ext. P") 
points(nutMV.ssa[!is.na(nutMV.ssa$P),], pch=21, bg=alpha("red", 0.6), cex=.6, col="black")
plot(country.ssa, col="darkgrey", ylim=bb[c(2,4)], xlim=bb[c(1,3)], main="ext. K") 
points(nutMV.ssa[!is.na(nutMV.ssa$K),], pch=21, bg=alpha("red", 0.6), cex=.6, col="black")
plot(country.ssa, col="darkgrey", ylim=bb[c(2,4)], xlim=bb[c(1,3)], main="ext. Mg") 
points(nutMV.ssa[!is.na(nutMV.ssa$Mg),], pch=21, bg=alpha("red", 0.6), cex=.6, col="black")
plot(country.ssa, col="darkgrey", ylim=bb[c(2,4)], xlim=bb[c(1,3)], main="ext. Fe") 
points(nutMV.ssa[!is.na(nutMV.ssa$Fe),], pch=21, bg=alpha("red", 0.6), cex=.6, col="black")
dev.off()

## Fig_AfNutrients_final_maps.png
g1km <- stack(list.files(pattern=glob2rx("*agg30cm*.tif")))
names(g1km) = sapply(names(g1km), function(x){strsplit(x, "_")[[1]][1]})
g1km <- as(g1km, "SpatialGridDataFrame")
#g1km$P.a <- (0.8*g1km$P + 0.2*g1km$P.O)/2
names(g1km)[which(names(g1km)=="ORCDRC")] = "C"
## back-transform:
g1km$C = g1km$C*1000
g1km$C = ifelse(is.na(g1km$N), NA, g1km$C)
## P.a, S, B, Cu, Zn
g1km$P = g1km$P/100
#g1km$S = g1km$S/100
g1km$Cu = g1km$Cu/100
g1km$B = g1km$B/100
g1km$Zn = g1km$Zn/100

png(file = "Fig_AfNutrients_final_maps_macro.png", res = 160, width = 1440, height = 1920)
tMV.varsP0 <- c("K", "Mg", "Ca", "Na", "P.T", "N")
tMV.varsP1 <- c("ext. K", "ext. Mg", "ext. Ca", "ext. Na", "P tot.", "Org. N")
par(mfrow=c(3,2), mar=c(0,0,3,0), oma=c(0,0,0,0))
for(i in 1:length(tMV.varsP0)){
  xr = quantile(g1km@data[,tMV.varsP0[i]], c(0.05,0.95), na.rm=TRUE)
  g1km$fix <- ifelse(g1km@data[,tMV.varsP0[i]]<xr[1], xr[1], ifelse(g1km@data[,tMV.varsP0[i]]>xr[2], xr[2], g1km@data[,tMV.varsP0[i]]))
  image(raster(g1km["fix"]), col=SAGA_pal[[1]][-c(1:6)], xlim=c(-18,52), ylim=c(-35,26), main=tMV.varsP1[i], asp=1, axes=FALSE, xlab="", ylab="", zlim=xr)
  #lines(country.af, col="black")
  lines(country.ssa, col="black")
  rx1 <- c(round(xr[2]), rep("", 12), round(xr[1]))
  legend("bottomleft", rx1, fill=rev(SAGA_pal[[1]][-c(1:5)]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
}
dev.off()

png(file = "Fig_AfNutrients_final_maps_micro.png", res = 160, width = 1440, height = 1920)
tMV.varsP0.m <- c("Zn", "B", "Cu", "Fe",  "Mn", "Al")
tMV.varsP1.m <- c("ext. Zn", "ext. B", "ext. Cu", "ext. Fe", "ext. Mn", "ext. Al")
par(mfrow=c(3,2), mar=c(0,0,3,0), oma=c(0,0,0,0))
for(i in 1:length(tMV.varsP0.m)){
  xr = quantile(g1km@data[,tMV.varsP0.m[i]], c(0.05,0.95), na.rm=TRUE)
  g1km$fix <- ifelse(g1km@data[,tMV.varsP0.m[i]]<xr[1], xr[1], ifelse(g1km@data[,tMV.varsP0.m[i]]>xr[2], xr[2], g1km@data[,tMV.varsP0.m[i]]))
  image(raster(g1km["fix"]), col=SAGA_pal[[1]][-c(1:6)], xlim=c(-18,52), ylim=c(-35,26), main=tMV.varsP1.m[i], asp=1, axes=FALSE, xlab="", ylab="", zlim=xr)
  lines(country.ssa, col="black")
  rx1 <- c(round(xr[2]), rep("", 12), round(xr[1]))
  legend("bottomleft", rx1, fill=rev(SAGA_pal[[1]][-c(1:6)]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
}
dev.off()

## plot using Leaflet:
library(leaflet)
library(raster)
library(htmlwidgets)
## macro-nutrients
for(j in 1:length(tMV.varsP0)){
  if(!file.exists(paste0(tMV.varsP0[j], "_5km.html"))){
    xr = quantile(g1km@data[,tMV.varsP0[j]], c(0.05,0.95), na.rm=TRUE)
    g1km$fix <- ifelse(g1km@data[,tMV.varsP0[j]]<=xr[1], xr[1], ifelse(g1km@data[,tMV.varsP0[j]]>=xr[2], xr[2], g1km@data[,tMV.varsP0[j]]))
    #r = projectRaster(raster(g1km["fix"]), res=5000, crs="+init=EPSG:3857", method="ngb")
    r = raster(g1km["fix"])
    pal <- colorNumeric(SAGA_pal[[1]], values(r), na.color = "transparent")
    m <- leaflet() %>% addTiles() %>% addRasterImage(r, colors=pal, opacity=0.6, project=FALSE) %>% addLegend(pal=pal, values=values(r), title=paste0(tMV.varsP1[j], " (ppm)"))
    saveWidget(m, file=paste0(tMV.varsP0[j], "_5km.html"))
  }
}

## micro-nutrients
tMV.varsP0.m <- c("Zn", "B", "Cu", "Fe",  "Mn", "Al")
tMV.varsP1.m <- c("ext. Zn", "ext. B", "ext. Cu", "ext. Fe", "ext. Mn", "ext. Al")
for(j in 1:length(tMV.varsP0.m)){
  if(!file.exists(paste0(tMV.varsP0.m[j], "_5km.html"))){
    xr = quantile(g1km@data[,tMV.varsP0.m[j]], c(0.05,0.95), na.rm=TRUE)
    g1km$fix <- ifelse(g1km@data[,tMV.varsP0.m[j]]<=xr[1], xr[1], ifelse(g1km@data[,tMV.varsP0.m[j]]>=xr[2], xr[2], g1km@data[,tMV.varsP0.m[j]]))
    #r = projectRaster(raster(g1km["fix"]), res=5000, crs="+init=EPSG:3857", method="ngb")
    r = raster(g1km["fix"])
    pal <- colorNumeric(SAGA_pal[[1]], values(r), na.color = "transparent")
    m <- leaflet() %>% addTiles() %>% addRasterImage(r, colors=pal, opacity=0.6, project=FALSE) %>% addLegend(pal=pal, values=values(r), title=paste0(tMV.varsP1.m[j], " (ppm)"))
    saveWidget(m, file=paste0(tMV.varsP0.m[j], "_5km.html"))
  }
}

tMV.varsP0.m <- c("N", "C", "P.T")
tMV.varsP1.m <- c("org. N", "org. C", "tot. P")
for(j in 1:length(tMV.varsP0.m)){
  if(!file.exists(paste0(tMV.varsP0.m[j], "_5km.html"))){
    xr = quantile(g1km@data[,tMV.varsP0.m[j]], c(0.05,0.95), na.rm=TRUE)
    g1km$fix <- ifelse(g1km@data[,tMV.varsP0.m[j]]<=xr[1], xr[1], ifelse(g1km@data[,tMV.varsP0.m[j]]>=xr[2], xr[2], g1km@data[,tMV.varsP0.m[j]]))
    r = projectRaster(raster(g1km["fix"]), res=5000, crs="+init=EPSG:3857", method="ngb")
    pal <- colorNumeric(SAGA_pal[[1]], values(r), na.color = "transparent")
    m <- leaflet() %>% addTiles() %>% addRasterImage(r, colors=pal, opacity=0.6, project=FALSE) %>% addLegend(pal=pal, values=values(r), title=paste0(tMV.varsP1.m[j], " (ppm)"))
    saveWidget(m, file=paste0(tMV.varsP0.m[j], "_5km.html"))
  }
}

## Clusters analysis results
g1kmc = readGDAL("/data/AFMicroNutrients/NCluster_M_AF_5km.tif")
g1kmc$SSEI = readGDAL("/data/AFMicroNutrients/SSI_NCluster_AF_5km.tif")$band1
#spplot(g1kmc)
pal.cl = rainbow(20)[sample.int(20,20)]
xr = range(g1kmc@data[,2], na.rm=TRUE)
png(file = "Fig_nutrient_clusters_map.png", res = 160, width = 1050, height = 1780)
par(mfrow=c(2,1), mar=c(0,0,3,0), oma=c(0,0,0,0))
image(raster(g1kmc[1]), col=pal.cl, main="Clusters", asp=1, axes=FALSE, xlab="", ylab="", cex.main=1)
lines(country.ssa, col="black")
op <- par(family="Courier")
legend("bottomleft", rev(paste0("c", 1:20)), fill=rev(pal.cl), horiz=FALSE, bty="n", cex=.7, y.intersp=.8)
par(op)
image(raster(g1kmc[2]), col=bpy.colors(25), main="Scaled Shannon Entropy Index", asp=1, axes=FALSE, xlab="", ylab="", cex.main=1)
lines(country.ssa, col="black")
legend("bottomleft", c(round(xr[2]), rep("", 23), round(xr[1])), fill=rev(bpy.colors(25)), horiz=FALSE, bty="n", cex=.8, y.intersp=.4)
dev.off()

## Fig_AfNutrients_correlations_USDA.pdf
#tif250m.usda <- stack(list.files(path="X:/ftp.soilgrids.org/data/recent", pattern=glob2rx("TAXOUSDA_*_250m_ll.tif$"), full.names = TRUE))
#taxousda <- raster("X:/ftp.soilgrids.org/data/recent/TAXOUSDA_250m_ll.tif")
#sel.nut <- nut[complete.cases(nut@data[,c("Mg","K")]),c("Mg","K","LOC_ID")]
#y = spTransform(nut[!duplicated(nut$LOC_ID),], CRS(proj4string(taxousda)))
#dfc2 <- raster::extract(y, x=taxousda)
g1km$Number <- as.factor(readGDAL("TAXOUSDA_5km_ll.tif")$band1)
levels(g1km$Number)
col.legend <- read.csv("TAXOUSDA_legend.csv")
g1km$TAXOUSDA <- plyr::join(data.frame(Number=g1km@data[["Number"]]), col.legend[,c("Number","Group")], type="left")$Group
summary(g1km$TAXOUSDA)
#boxplot(K~TAXOUSDA, g1km@data)

cl_dfc <- data.frame(TAXOUSDA = c(paste(g1km$TAXOUSDA), paste(g1km$TAXOUSDA)), var=c(rep("ext_Fe", nrow(g1km)), rep("org_N", nrow(g1km))), content=c(g1km$Fe, g1km$N))
cl_dfc <- cl_dfc[complete.cases(cl_dfc),]
str(cl_dfc)

## plot differences in nutrients:
sel <- cl_dfc$TAXOUSDA %in% c("Aqualfs", "Ustalfs", "Udults", "Ustolls", "Usterts", "Fluvents", "Calcids", "Xeralfs", "Psamments", "Orthents", "Ustox","Humults")
summary(sel)
qplot(TAXOUSDA, content, fill=factor(var), data=cl_dfc[sel,], geom="boxplot", log="y")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Correlation between TAXOUSDA and nutrients:
g1km$Usterts <- readGDAL("TAXOUSDA_Usterts_5km_ll.tif")$band1
g1km$Humults <- readGDAL("TAXOUSDA_Humults_5km_ll.tif")$band1
g1km$Usterts <- ifelse(is.na(g1km$N), NA, g1km$Usterts)
g1km$Humults <- ifelse(is.na(g1km$N), NA, g1km$Humults)

png(file = "Fig_AfNutrients_correlations_Usterts.png", res = 160, width = 1820, height = 820)
par(mfrow=c(1,2), mar=c(0,0,3,0), oma=c(0,0,0,0))
xr = quantile(g1km@data[,"Mg"], c(0.05,0.95), na.rm=TRUE)
g1km$fix <- ifelse(g1km@data[,"Mg"]<xr[1], xr[1], ifelse(g1km@data[,"Mg"]>xr[2], xr[2], g1km@data[,"Mg"]))
image(raster(g1km["fix"]), col=SAGA_pal[[1]][-c(1:5)], xlim=c(-18,52), ylim=c(-35,26), main="Ext. Mg", asp=1, axes=FALSE, xlab="", ylab="", zlim=xr)
lines(country.ssa, col="black")
rx1 <- c(round(xr[2]), rep("", 13), round(xr[1]))
legend("bottomleft", rx1, fill=rev(SAGA_pal[[1]][-c(1:5)]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
g1km$fix2 <- ifelse(g1km$Usterts<0, 0, ifelse(g1km$Usterts>30, 30, g1km$Usterts))
image(raster(g1km["fix2"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], xlim=c(-18,52), ylim=c(-35,26), main="Usterts", asp=1, axes=FALSE, xlab="", ylab="", zlim=c(0,30))
lines(country.ssa, col="black")
legend("bottomleft", c(30, rep("", 18), 0), fill=rev(SAGA_pal[["SG_COLORS_YELLOW_RED"]]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
dev.off()

png(file = "Fig_AfNutrients_correlations_Humults.png", res = 160, width = 1820, height = 820)
par(mfrow=c(1,2), mar=c(0,0,3,0), oma=c(0,0,0,0))
xr = quantile(g1km@data[,"N"], c(0.05,0.95), na.rm=TRUE)
g1km$fix <- ifelse(g1km@data[,"N"]<xr[1], xr[1], ifelse(g1km@data[,"N"]>xr[2], xr[2], g1km@data[,"N"]))
image(raster(g1km["fix"]), col=SAGA_pal[[1]][-c(1:5)], xlim=c(-18,52), ylim=c(-35,26), main="Org. N", asp=1, axes=FALSE, xlab="", ylab="", zlim=xr)
lines(country.ssa, col="black")
rx1 <- c(round(xr[2]), rep("", 13), round(xr[1]))
legend("bottomleft", rx1, fill=rev(SAGA_pal[[1]][-c(1:5)]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
g1km$fix2 <- ifelse(g1km$Humults<0, 0, ifelse(g1km$Humults>30, 30, g1km$Humults))
image(raster(g1km["fix2"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], xlim=c(-18,52), ylim=c(-35,26), main="Humults", asp=1, axes=FALSE, xlab="", ylab="", zlim=c(0,30))
lines(country.ssa, col="black")
legend("bottomleft", c(30, rep("", 18), 0), fill=rev(SAGA_pal[["SG_COLORS_YELLOW_RED"]]), horiz=FALSE, bty="n", cex=1, y.intersp=.4)
dev.off()


