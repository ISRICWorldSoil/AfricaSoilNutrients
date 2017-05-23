## Functions for Af soil nutrient mapping
## by: Tom.Hengl@isric.org

kern.density <- function(i, input.points, wowin, COLUMN, sigma=80000, out.dir, ext){
  out.tif = paste0(out.dir, "/dens_", i, "_", ext, ".tif")
  if(!file.exists(out.tif)){
    s = input.points[grep(i, paste(input.points@data[,COLUMN]), ignore.case = FALSE, fixed = TRUE),]
    if(length(s)>60){
      s$mag = 1
      s.ppp <- ppp(s@coords[,1], s@coords[,2], marks=s$mag, window=wowin)
      densMAG <- density.ppp(s.ppp, sigma=sigma, weights=s.ppp$marks)
      densMAG = as(densMAG, "SpatialGridDataFrame")
      densMAG$vf = densMAG$v/max(densMAG$v, na.rm=TRUE)*100
      #plot(raster(densMAG["vf"]), col=rev(bpy.colors(30)))
      proj4string(densMAG) = proj4string(input.points)
      writeGDAL(densMAG["vf"], out.tif, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      gc()
    } else {
      message('Less than 60 points... skipping')
    }
  }
}

extract.tif = function(x, y){
  r = raster(x)
  if(!is.na(proj4string(r))){
    y = spTransform(y, proj4string(r))
  }  
  out = raster::extract(y=y, x=r)
  return(out)
}

gdalwarp_250m <- function(x, out.dir="./stacked250m/"){
  out = paste0(out.dir, gsub(".sdat", ".tif", gsub(".vrt", ".tif", basename(x))))
  if(!file.exists(out)){ 
    if(length(grep(x, pattern = "SKY"))>0|length(grep(x, pattern = "dens_"))>0){
      resm = "cubicspline"
    } else {
      if(length(grep(x, pattern = "LC-L4"))>0|length(grep(x, pattern = "GES"))>0|length(grep(x, pattern = "LF_Desc"))>0){
        resm = "near"
      } else {
        resm = "bilinear"
      }
    }
    system(paste0(gdalwarp, ' ', x, ' ', out, ' -co \"COMPRESS=DEFLATE\" -r \"', resm, '\" -tr ', cellsize, ' ', cellsize, ' -co \"BIGTIFF=YES\" -multi -wm 2000 -t_srs \"', af.prj, '\" -te ', paste0(extent(r250m)[c(1,3,2,4)], collapse=" ")))
  }
}

summary_GAUL_tiles <- function(i, tile.tbl, out.path="/data/AFMicroNutrients/tiled", cnt="/data/AFMicroNutrients/stacked250m/GAUL_COUNTRIES.tif", countries){
  out.file = paste0(out.path, "/T", tile.tbl[i,"ID"], "/GAUL_T", tile.tbl[i,"ID"], ".tif")
  m = readGDAL(fname=cnt, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
  m$band1 = ifelse(m$band1<0, NA, m$band1)
  if(sum(is.na(m$band1))<nrow(m)){
    m$GAUL_COUNTRY = factor(m$band1, levels=as.character(1:nrow(countries)), labels=countries$NAMES)
    out = summary(m$GAUL_COUNTRY[!is.na(m$GAUL_COUNTRY)], maxsum = nrow(countries))
    out = data.frame(Value=attr(out, "names"), Count=as.numeric(out))
    out = out[out$Count>0,]
    out$ID = i
    if(!file.exists(out.file)){
      writeGDAL(m[1], out.file, type="Int16", options="COMPRESS=DEFLATE", mvFlag = -32768)
    }
    return(out)
  }
}

make_newdata <- function(i, tile.tbl, in.path="/data/AFMicroNutrients/stacked250m", out.path="/data/AFMicroNutrients/tiled", pr.lst, covs.tif, mask="/data/AFMicroNutrients/stacked250m/GAUL_COUNTRIES.tif", sel.country, ov.quant, GESUSG6.lev, LC.lev, LF.lev){
  out.rds <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/T", tile.tbl[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(fname=mask, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    names(m) = "GAUL_COUNTRIES"
    m = as(m, "SpatialPixelsDataFrame")
    sel.p = m$GAUL_COUNTRIES %in% sel.country
    if(sum(sel.p)>1){
      m = m[which(sel.p),] 
      x = spTransform(m, CRS("+proj=longlat +datum=WGS84"))
      m$LONWGS84 <- x@coords[,1]
      m$LATWGS84 <- x@coords[,2]
      for(j in 1:length(covs.tif)){
        m@data[,j+1] <- signif(readGDAL(covs.tif[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index], 4)
      }
      names(m)[2:(length(covs.tif)+1)] = file_path_sans_ext(gsub("-", "_", basename(covs.tif)))
      m$ID = m@grid.index
      m$GESUSG6 = factor(paste(m$GESUSG6), levels = GESUSG6.lev)
      m.GESUSG6 = data.frame(model.matrix(~GESUSG6-1, m@data))
      names(m.GESUSG6) = make.names(names(m.GESUSG6))
      m.GESUSG6 <- m.GESUSG6[match(rownames(m@data),rownames(m.GESUSG6)),]
      m.GESUSG6$ID = m$ID
      m$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1 = factor(paste(m$ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1), levels = LC.lev)
      m.LC = data.frame(model.matrix(~ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1-1, m@data))
      names(m.LC) = make.names(names(m.LC))
      m.LC <- m.LC[match(rownames(m@data),rownames(m.LC)),]
      m.LC$ID = m$ID
      m$EF_LF_Desc_250m = factor(paste(m$EF_LF_Desc_250m), levels = LF.lev)
      m.LF = data.frame(model.matrix(~EF_LF_Desc_250m-1, m@data))
      names(m.LF) = make.names(names(m.LF))
      m.LF <- m.LF[match(rownames(m@data),rownames(m.LF)),]
      m.LF$ID = m$ID
      m@data = plyr::join_all(list(m@data, m.GESUSG6, m.LC, m.LF), by="ID") #[,pr.lst]
      ## mask out water bodies (for this best use landsat NIR image):
      m = m[m$Landsat2014_NIR>17,]
      ## Fill-in missing values (if necessary):
      sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
      if(sum(sel.mis)>0){
        for(j in which(sel.mis)){
          if(!is.factor(m@data[,j])){
            if(length(grep(pattern="MNGMAP", names(m)[j]))>0|length(grep(pattern="COAUSG6", names(m)[j]))>0|length(grep(pattern="giems_d15_v10", names(m)[j]))>0|length(grep(pattern="Water_", names(m)[j]))>0|length(grep(pattern="GESUSG", names(m)[j]))>0){  
              repn = 0
            } else {
              repn = quantile(m@data[,j], probs=.5, na.rm=TRUE)
              if(is.na(repn)){
                repn = ov.quant[[names(m)[j]]]
              }
            }
            m@data[,j] = ifelse(is.na(m@data[,j]), repn, m@data[,j])
          }
        }
      }
      saveRDS(m, out.rds)
      gc(); gc()
    }
  }
}

agg_layers <- function(i, varn, d = c(0,5,15,30), in.path="/data/AFMicroNutrients/tiled", ot="Int16", dstnodata=-32768){
  in.lst <- paste0(path=in.path, "/", i, "/", varn, "_M_sl", 1:length(d), "_", i, ".tif")
  out.tif <- paste0(path=in.path, "/", i, "/", varn, "_M_agg", d[length(d)]-d[1], "cm_", i, ".tif")
  if(!file.exists(out.tif)&all(file.exists(in.lst))){
    s <- raster::stack(in.lst)
    s <- as(s, "SpatialGridDataFrame")
    ## Trapezoidal rule
    for(l in 1:(length(d)-1)){
      s@data[,paste0("sd",l)] <- rowMeans(s@data[,l:(l+1)])
    }
    s$sum_sd <- rowSums(as.matrix(s@data[,paste0("sd",1:(length(d)-1))]) %*% diag(diff(d)))/sum(diff(d))
    writeGDAL(s["sum_sd"], out.tif, type=paste(ot), mvFlag=paste(dstnodata), options="COMPRESS=DEFLATE")
    gc(); gc()
  }
}

tile_SoilGrids = function(i, tile.tbl, varn, in.tifs, out.path="/data/AFMicroNutrients/tiled", t_srs){
  out.tifs = paste0(out.path, "/T", i, "/", gsub("_250m_ll.tif", paste0("_T", i, ".tif"), basename(in.tifs)))
  if(any(!file.exists(out.tifs))){
    sapply(1:length(in.tifs), function(x){ system(paste0(gdalwarp, ' ', in.tifs[x], ' ', out.tifs[x], ' -co \"COMPRESS=DEFLATE\" -t_srs \"', t_srs, '\" -te ', tile.tbl[i,"xl"], ' ', tile.tbl[i,"yl"], ' ', tile.tbl[i,"xu"], ' ', tile.tbl[i,"yu"], ' -tr 250 250'), intern = TRUE, show.output.on.console = FALSE)})
  }
}

make_mosaick <- function(i, varn, tile.name, in.path="/data/AFMicroNutrients/tiled", ot="Int16", dstnodata=-32768, compress=TRUE){
  out.tif <- paste0("/data/GEOG/Nutrients/", varn, '_', i, '_', tile.name, '_250m.tif')
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
    if(compress==TRUE){
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -ot \"', paste(ot), '\" -dstnodata \"',  paste(dstnodata), '\" -r \"near\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\"')) 
    } else {
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -ot \"', paste(ot), '\" -dstnodata \"',  paste(dstnodata), '\" -r \"near\" -co \"BIGTIFF=YES\" -wm 2000'))
    }
    system(paste0(gdaladdo, ' ', out.tif, ' 2 4 8 16 32 64 128'))
    system(paste0(gdal_translate, ' -of GTiff -r \"average\" -tr 1000 1000 ', vrt.tmp, ' ', gsub("_250m.tif", "_1km.tif", out.tif), ' -ot \"', paste(ot), '\" -a_nodata \"', paste(dstnodata), '\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
    unlink(vrt.tmp)
    unlink(out.tmp)
  }
}

predict_mcluster_tiles <- function(i, gm, tile.tbl, in.path="/data/AFMicroNutrients/tiled", out.path="/data/AFMicroNutrients/tiled", nlevs){
  out.tifs <- paste0(out.path, "/", i, "/NCluster_", 1:nlevs, "_", i, ".tif")
  if(any(!file.exists(out.tifs))){
    rds.file = paste0(in.path, "/", i, "/", i,".rds")
    if(file.exists(rds.file)){
      m = readRDS(rds.file)
      ## add missing column:
      m$PHIHOX_M = rowMeans(m@data[,grep("PHIHOX", names(m))], na.rm=TRUE)
      m@data = data.frame(predict(gm, m@data)$predictions)
      for(j in 1:ncol(m)){
        out.tif <- paste0(out.path, "/", i, "/NCluster_", j, "_", i, ".tif")
        m@data[,j] <- round(m@data[,j]*100)
        writeGDAL(m[j], out.tif, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
      }
      ## match most probable class
      m$cl <- apply(m@data,1,which.max)
      writeGDAL(m["cl"], paste0(out.path, "/", i, "/NCluster_M_", i, ".tif"), type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
    }
    gc()
  }
}

entropy_tile <- function(i, in.path, varn, levs){
  out.p <- paste0(in.path, "/", i, "/SSI_", varn, "_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- paste0(in.path, "/", i, "/", varn, "_", levs, "_", i, ".tif")
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    gc()
    v <- unlist(alply(s@data, 1, .fun=function(x){entropy.empirical(unlist(x))})) 
    s$SSI <- round(v/entropy.empirical(rep(1/length(levs),length(levs)))*100)
    writeGDAL(s["SSI"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

gdalwarp_1km <- function(x, out.dir="./stacked1km/", cellsize=1000, r1km){
  out = paste0(out.dir, gsub(".sdat", ".tif", gsub(".vrt", ".tif", basename(x))))
  if(!file.exists(out)){ 
    system(paste0('gdalwarp ', x, ' ', out, ' -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -co \"BIGTIFF=YES\" -multi -wm 2000 -t_srs \"', proj4string(r1km), '\" -te ', paste0(extent(r1km)[c(1,3,2,4)], collapse=" ")))
  }
}
