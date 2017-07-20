## OCP data preparation of soil properties for modeling (https://github.com/ISRICWorldSoil/AfricaSoilNutrients)
## Data sources: SoilGrids250km global (ftp://ftp.soilgrids.org/data/recent/) and Africa soil nutrients (http://gsif.isric.org/doku.php/wiki:africa_nutrient_maps)
## prepared by: Tom.Hengl@isric.org and Dennis.Waalvoort@wur.nl

library(rgdal)
library(raster)
## Spatial resolution 250m in latlon:
raster("/data/GEOG/TAXNWRB_250m_ll.tif")

## Mali and Senegal study areas:
system('ogrinfo ../shapes/ML_OCP_AOI123.shp')
system(paste0('saga_cmd -c=56 grid_gridding 0 -INPUT \"../shapes/ML_OCP_AOI123.shp\" -FIELD \"CULTURE\" -GRID \"ML_OCP_AOI123.sgrd\" -TARGET_DEFINITION 0 -TARGET_USER_SIZE 0.002083333'))
system(paste0('saga_cmd -c=56 grid_gridding 0 -INPUT \"../shapes/SN_OCP_AOI123_IRR.shp\" -FIELD \"CULTURE\" -GRID \"SN_OCP_AOI123_IRR.sgrd\" -TARGET_DEFINITION 0 -TARGET_USER_SIZE 0.002083333'))

## Aggregate values from SoilGrids vertically using the Trapezoidal rule (http://dx.doi.org/10.1371/journal.pone.0169748):
agg_layers <- function(tif, tif.out, d=c(0,5,15,30), ot="Int16", dstnodata=-32768){
  s <- stack(tif)
  s <- as(s, "SpatialGridDataFrame")
  ## Trapezoidal rule
  for(j in 1:(length(d)-1)){
    s@data[,paste0("sd",j)] <- rowMeans(s@data[,j:(j+1)])
  }
  s$sum_sd <- rowSums(as.matrix(s@data[,paste0("sd",1:(length(d)-1))]) %*% diag(diff(d)))/sum(diff(d))
  writeGDAL(s["sum_sd"], tif.out, type=ot, mvFlag=dstnodata, options="COMPRESS=DEFLATE")
}

## clip soilgrids map using local area:
clip_sg = function(tif, mask, tif.out){
  if(!file.exists(tif.out)){
    r = raster(mask)
    if(length(tif)>1){
      for(j in tif){
        system(paste0('gdalwarp ', j, ' ', basename(j), ' -t_srs \"', proj4string(r),'\" -co \"COMPRESS=DEFLATE\" -tr ', res(r)[1], ' ', res(r)[2],' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
      }
      agg_layers(basename(tif), tif.out, d=c(0,5,15,30))
      unlink(basename(tif))
    } else {
      system(paste0('gdalwarp ', tif, ' ', tif.out, ' -t_srs \"', proj4string(r),'\" -co \"COMPRESS=DEFLATE\" -tr ', res(r)[1], ' ', res(r)[2],' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
    }
    ## clip:
    r2 = readGDAL(tif.out)
    r2$mask = readGDAL(mask)$band1
    r2$band1 = ifelse(is.na(r2$mask), NA, r2$band1)
    writeGDAL(r2["band1"], tif.out, type = "Int16", options = c("COMPRESS=DEFLATE"), mvFlag = -32768)
  }
}

## test 1
#mask = "ML_OCP_AOI123.sdat"
#tif="/data/GEOG/Nutrients/K_M_agg30cm_AF_250m.tif"
#tif.out="K_M_agg30cm_ML_OCP_AOI123_250m.tif"
## test 2
#tif=paste0("/data/GEOG/PHIHOX_M_sl", 1:4, "_250m_ll.tif")
#tif.out="PHIHOX_M_agg30cm_ML_OCP_AOI123_250m.tif"

## list of layers to be processed:
tif.lst = c(lapply(c("ORCDRC","PHIHOX","CLYPPT","CECSOL"), function(i){paste0("/data/GEOG/", i, "_M_sl", 1:4,"_250m_ll.tif")}), as.list(paste0("/data/GEOG/Nutrients/", c("Al","Fe","Mg","Ca","P.B","K","Na","N", "P", "Cu", "Zn", "B"), "_M_agg30cm_AF_250m.tif")))
## 16 in total
## ORCDRC layer sl1 in error! (https://github.com/ISRICWorldSoil/SoilGrids250m/issues/35)
tif.lst[[1]][1] = tif.lst[[1]][2]

## run in loop (it could also work in parallel but then the file names need to be tweaked)
## TAKES ca 3 mins
for(i in 1:length(tif.lst)){
  for(k in c("ML_OCP_AOI123.sdat","SN_OCP_AOI123_IRR.sdat")){
    tif.out = paste0(strsplit(basename(tif.lst[[i]][1]), "_")[[1]][1], "_M_agg30cm_", strsplit(k, ".sdat")[[1]][1], "_250m.tif")
    clip_sg(tif=tif.lst[[i]], mask=k, tif.out)
  }
}
