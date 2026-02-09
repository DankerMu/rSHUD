#' Write ESRI shapefile out
#' \code{writeshape}
#' @param shp Spatial file
#' @param crs projection
#' @param file file path, without '.shp'.
#' @export
#' @examples
#' library(sp)
#' library(rgeos)
#' library(rgdal)
#' sp1 = readWKT("POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))")
#' raster::crs(sp1) =sp::CRS("+init=epsg:4326")
#' writeshape(sp1, file=file.path(tempdir(), 'sp1'))
#' sp2=readOGR(file.path(tempdir(), 'sp1.shp'))
#' plot(sp2)
writeshape <- function(shp, file = NULL, crs = NULL) {
  msg = "writeshape::"

  if (is.null(file) || !nzchar(file)) {
    return(invisible(NULL))
  }
  if (tolower(tools::file_ext(file)) == "shp") {
    file <- sub("\\.[Ss][Hh][Pp]$", "", file)
  }

  path <- dirname(file)
  fn <- basename(file)
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }

  v <- .rshud_as_spatvector(shp, crs = crs)
  val <- tryCatch(terra::values(v, dataframe = TRUE), error = function(e) NULL)
  if (is.null(val) || ncol(val) == 0) {
    terra::values(v) <- data.frame(ID = seq_len(nrow(v)))
  }

  out <- file.path(path, paste0(fn, ".shp"))
  terra::writeVector(v, out, overwrite = TRUE)
  message(msg, out, " is saved")
  invisible(out)
}
#' Re-project coordinates betwen GCS and PCS
#' \code{ProjectCoordinate} 
#' @param  x 2-column matrix of coordinates.
#' @param  proj4string proj4string
#' @param  P2G if TRUE, a cartographic projection into lat/long, otherwise projects from lat/long into a cartographic projection.
#' @return Basic model infomation, figures and tables
#' @export
ProjectCoordinate <- function(x, proj4string, P2G=TRUE){
  x <- as.matrix(x)
  if (ncol(x) != 2) stop("ProjectCoordinate: x must be a 2-column matrix")

  crs_out <- .rshud_as_crs(proj4string)
  if (P2G) {
    pts <- sf::st_as_sf(data.frame(X = x[, 1], Y = x[, 2]), coords = c("X", "Y"), crs = crs_out)
    pts <- sf::st_transform(pts, 4326)
    y <- sf::st_coordinates(pts)
    colnames(y) <- c("Lon", "Lat")
  } else {
    pts <- sf::st_as_sf(data.frame(Lon = x[, 1], Lat = x[, 2]), coords = c("Lon", "Lat"), crs = 4326)
    pts <- sf::st_transform(pts, crs_out)
    y <- sf::st_coordinates(pts)
    colnames(y) <- c("X", "Y")
  }
  y
}
#' SpatialData to Raster
#' \code{sp2raster}
#' @param sp SpatialPolygon
#' @param mask Raster mask of mesh domain.
#' @param ngrids Number of grid along x direction.
#' @param resolution Resolution, defaul = NULL, resolution = extent / ngrids
#' @param field Index of field
#' @return Raster map
#' @export
sp2raster <- function (sp, mask = get('MASK', envir = .shud),
                       ngrids=200, 
                       resolution=NULL, field=1) {
  if( is.null(mask) ){
    ext <-  raster::extent (sp)
    xlim=ext[1:2]
    ylim=ext[3:4]
    if ( resolution<=0 || is.null(resolution)){
      dx=diff(xlim) / ngrids;
    }else{
      dx=resolution
    }
    r <- raster::raster(ext, res=dx)
  }else{
    r = mask
  }
  ## Rasterize the shapefile
  rr <-raster::rasterize(sp, r, field=field)
  return(rr)
}

#' Generate the raster mask of Mesh domain
#' \code{shud.mask}
#' @param pm \code{shud.mesh}
#' @param n Number of grid
#' @param rr Default mask in .shud environment
#' @param cellsize Resolution, defaul = NULL, resolution = extent / ngrids
#' @param proj Projection parameter
#' @return Raster map
#' @export
shud.mask  <- function (pm = readmesh(), proj=NULL,
                        rr = get('MASK', envir=.shud),
                        n=10000, cellsize=NULL){
  # mesh=readmesh(shp=TRUE); ngrids=100; resolution=0
  if(is.null(rr)){
    spm =sp.mesh2Shape(pm)
    sp0=rgeos::gUnaryUnion(spm)
    if(is.null(cellsize)){
      # grd <- as.data.frame(sp::spsample(spm, "regular", n=n))
      # grd <- as.data.frame(sp::spsample(spm, "regular", nsig=2, n=n))
      grd <- sp::makegrid(sp0, n = n)
    }else{
      # grd <- as.data.frame(sp::spsample(spm, "regular", cellsize = cellsize))
      grd <- sp::makegrid(sp0, cellsize = cellsize, pretty = FALSE)
    }
    names(grd)       <- c("X", "Y")
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object
    rr=raster::raster(grd); rr[]=1
    rr=raster::mask(rr, sp0)
    if(!is.null(proj)){
      raster::crs(rr) = proj
    }
    assign('MASK', rr, envir=.shud)
  }else{
    rr = rr
  }
  rr
}

#' SpatialData to Raster
#' \code{MeshData2Raster}
#' @param x vector or matrix, length/nrow is number of cells.
#' @param rmask mask of the mesh file
#' @param stack Whether export the stack, only when the x is a matrix, i.e. (Ntime x Ncell).
#' @param proj Projejction parameter
#' @param pm shud mesh
#' @param method method for interpolation, default = 'idw'
#' @param plot Whether plot the result.
#' @return Raster map
#' @export
MeshData2Raster <- function(x=getElevation(),
                            rmask=shud.mask(proj=proj), 
                            pm=readmesh(), proj=NULL,
                            stack=FALSE, method='ide',
                            plot =FALSE){
  
  if(stack){
    ret <- raster::stack(apply(x, 1, FUN = MeshData2Raster) )
  }else{
    if( is.matrix(x) | is.data.frame(x)){
      x = as.numeric(x[nrow(x),])
    }
    if(any(is.na(x)) ){
      x[is.na(x)] = 0
    }
    if (any(is.infinite(x))){
      x[is.infinite(x)] = 0
    }
    xy=getCentroid(pm=pm)[,2:3]
    
    if(grepl('idw', tolower(method))){
      val= data.frame(xy, x)
      colnames(val) = c('X', 'Y', 'Z')
      sp::coordinates(val) = c('X', 'Y')
      grd=methods::as(rmask, 'SpatialGrid')
      # if(grepl(method, 'idw')){
      # Interpolate the grid cells using a power value of 2 (idp=2.0)
      dat <- gstat::idw(Z ~ 1, val, newdata=grd, idp=2.0)
      r = raster::raster(dat)
    }
    if(grepl('linear', tolower(method))){
      xr = raster::rasterToPoints(rmask)
      ext=raster::extent(rmask);res=raster::res(rmask); hr = res/2
      r0=rmask; r0[]=1
      xyo=raster::rasterToPoints(r0)
      xx=interp::interp(x=xy[,1], y=xy[,2], z=x, xo=xyo[,1], yo=xyo[,2] )
      r = raster::setValues(rmask, as.numeric(xx$z))
    }
    if(grepl('ide', tolower(method))){
      tps <- fields::Tps(xy, x)
      # use model to predict values at all locations
      r <- raster::interpolate(rmask, tps)
    }
    ret <- raster::mask(r,rmask)
  }
  if(plot){
    raster::plot(ret)
  }
  if(!is.null(proj)){ raster::crs(ret) <- proj }
  return(ret)
}


#' Remove the holes in polygons
#' \code{removeholes}
#' @param sp SpatialPolygons or SpatialPolygonDataFrame
#' @return Raster map
#' @export
#' @examples
#' library(sp)
#' p.out = Polygon(cbind(c(4,4,6,7,4),c(5,3,2,5,5))  )
#' p.hole = Polygon(cbind(c(5,6,6,5,5),c(4,4,3,3,4) ), hole = TRUE)
#' sp <- SpatialPolygons(list(Polygons(list(p.out, p.hole), "1")))
#' s = removeholes(sp)
#' par(mfrow=c(1,2))
#' plot(sp)
#' plot(s)
removeholes <- function(sp){
  if (inherits(sp, "SpatVector")) {
    x_sf <- sf::st_as_sf(sp)
    y_sf <- removeholes(x_sf)
    return(terra::vect(y_sf))
  }

  if (inherits(sp, "sf") || inherits(sp, "sfc")) {
    x_sf <- if (inherits(sp, "sf")) sp else sf::st_sf(geometry = sp)
    geom <- sf::st_geometry(x_sf)
    new_geom <- sf::st_sfc(lapply(geom, .rshud_geom_drop_holes), crs = sf::st_crs(x_sf))
    if (inherits(sp, "sfc")) return(new_geom)
    x_sf$geometry <- new_geom
    return(x_sf)
  }

  if (inherits(sp, "SpatialPolygons") || inherits(sp, "SpatialPolygonsDataFrame")) {
    x <- sp
    nx <- length(x)
    rl <- vector("list", nx)
    for (i in seq_len(nx)) {
      spg <- x@polygons[[i]]@Polygons
      npg <- length(spg)
      ypg <- list()
      k <- 1
      for (j in seq_len(npg)) {
        if (!spg[[j]]@hole) {
          ypg[[k]] <- spg[[j]]
          k <- k + 1
        }
      }
      rl[[i]] <- sp::Polygons(ypg, ID = as.character(i))
    }
    ret <- sp::SpatialPolygons(rl)
    crs0 <- tryCatch(sp::proj4string(x), error = function(e) NA_character_)
    if (!is.na(crs0) && nzchar(crs0)) sp::proj4string(ret) <- crs0
    return(ret)
  }

  stop("removeholes: unsupported input type: ", paste(class(sp), collapse = "/"))
}
#' Generatue fishnet
#' \code{fishnet} Generate fishnet by the coordinates
#' @param xx  x coordinates
#' @param yy  y coordinates
#' @param crs projections parameters, defaul = epsg:4326
#' @param type option = 'polygon', 'points', 'line'
#' @return spatildata (.shp)
#' @export
#' @examples
#' library(raster)
#' ext=c(0, 8, 2, 10)
#' dx = 2; dy = 4
#' xx=seq(ext[1], ext[2], by=dx)
#' yy=seq(ext[3], ext[4], by=dy)
#' sp1 =fishnet(xx=ext[1:2], yy=ext[3:4])
#' sp2 =fishnet(xx=xx + .5 * dx, yy=yy + 0.5 * dy)
#' sp3 =fishnet(xx=xx, yy=yy, type = 'point')
#' plot(sp1, axes=TRUE, xlim=c(-1, 1)*dx +ext[1:2], ylim=c(-1, 1)*dy + ext[3:4])
#' plot(sp2, axes=TRUE, add=TRUE, border=2)
#' plot(sp3, axes=TRUE, add=TRUE, col=3, pch=20)
#' grid()
fishnet <- function(xx, yy,
                    crs = 4326,
                    type = "polygon") {
  type0 <- tolower(type)
  crs_sf <- .rshud_as_crs(crs)
  xx <- as.numeric(xx)
  yy <- as.numeric(yy)
  if (length(xx) == 0 || length(yy) == 0) stop("fishnet: xx/yy cannot be empty")

  if (grepl("line", type0)) {
    ymin <- min(yy)
    ymax <- max(yy)
    xmin <- min(xx)
    xmax <- max(xx)

    n <- length(xx) + length(yy)
    geoms <- vector("list", n)
    x1 <- numeric(n)
    y1 <- numeric(n)
    x2 <- numeric(n)
    y2 <- numeric(n)

    k <- 1
    for (i in seq_along(xx)) {
      x1[[k]] <- xx[[i]]
      y1[[k]] <- ymin
      x2[[k]] <- xx[[i]]
      y2[[k]] <- ymax
      geoms[[k]] <- sf::st_linestring(matrix(c(x1[[k]], y1[[k]], x2[[k]], y2[[k]]), ncol = 2, byrow = TRUE))
      k <- k + 1
    }
    for (j in seq_along(yy)) {
      x1[[k]] <- xmin
      y1[[k]] <- yy[[j]]
      x2[[k]] <- xmax
      y2[[k]] <- yy[[j]]
      geoms[[k]] <- sf::st_linestring(matrix(c(x1[[k]], y1[[k]], x2[[k]], y2[[k]]), ncol = 2, byrow = TRUE))
      k <- k + 1
    }

    df <- data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
    return(sf::st_sf(df, geometry = sf::st_sfc(geoms, crs = crs_sf)))
  }

  if (grepl("point", type0)) {
    xm <- expand.grid(X = xx, Y = yy)
    ret <- sf::st_as_sf(xm, coords = c("X", "Y"), crs = crs_sf, remove = FALSE)
    return(ret)
  }

  if (!grepl("polygon", type0)) stop("fishnet: type must be polygon|point|line")

  xx <- sort(unique(xx))
  yy <- sort(unique(yy))
  nx <- length(xx)
  ny <- length(yy)
  if (nx < 2 || ny < 2) stop("fishnet: polygon requires at least 2 x and 2 y values")

  ncell <- (nx - 1) * (ny - 1)
  geoms <- vector("list", ncell)
  xmin <- numeric(ncell)
  xmax <- numeric(ncell)
  ymin <- numeric(ncell)
  ymax <- numeric(ncell)
  xcenter <- numeric(ncell)
  ycenter <- numeric(ncell)

  k <- 1
  for (i in seq_len(nx - 1)) {
    for (j in seq_len(ny - 1)) {
      xi1 <- xx[[i]]
      xi2 <- xx[[i + 1]]
      yj1 <- yy[[j]]
      yj2 <- yy[[j + 1]]

      xmin[[k]] <- min(xi1, xi2)
      xmax[[k]] <- max(xi1, xi2)
      ymin[[k]] <- min(yj1, yj2)
      ymax[[k]] <- max(yj1, yj2)
      xcenter[[k]] <- (xmin[[k]] + xmax[[k]]) / 2
      ycenter[[k]] <- (ymin[[k]] + ymax[[k]]) / 2

      coords <- matrix(
        c(xmin[[k]], ymin[[k]],
          xmin[[k]], ymax[[k]],
          xmax[[k]], ymax[[k]],
          xmax[[k]], ymin[[k]],
          xmin[[k]], ymin[[k]]),
        ncol = 2,
        byrow = TRUE
      )
      geoms[[k]] <- sf::st_polygon(list(coords))
      k <- k + 1
    }
  }

  df <- data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    xcenter = xcenter,
    ycenter = ycenter
  )
  ret <- sf::st_sf(df, geometry = sf::st_sfc(geoms, crs = crs_sf))
  return(ret)
}

#' Add holes into Polygons
#' \code{AddHoleToPolygon}
#' @param poly SpatialPolygons
#' @param hole Hole Polygon
#' @export
AddHoleToPolygon <-function(poly,hole){
  # https://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
  # invert the coordinates for Polygons to flag it as a hole
  coordsHole <-  hole@polygons[[1]]@Polygons[[1]]@coords
  newHole <- sp::Polygon(coordsHole,hole=TRUE)
  
  # punch the hole in the main poly
  listPol <- poly@polygons[[1]]@Polygons
  listPol[[length(listPol)+1]] <- newHole
  punch <-sp::Polygons(listPol,poly@polygons[[1]]@ID)
  
  # make the polygon a SpatialPolygonsDataFrame as the entry
  new <- sp::SpatialPolygons(list(punch),proj4string=poly@proj4string)
  new <- sp::SpatialPolygonsDataFrame(new,data=as(poly,"data.frame"))
  return(new)
}
#' Cut sptialLines with threshold.
#' \code{sp.CutSptialLines}
#' @param sl SpatialLines or SpatialLineDataFrame
#' @param tol Tolerence. If the length of segment is larger than tolerance, cut the segment until the maximum segment is shorter than tolerance.
#' @export
#' @examples
#' library(rSHUD)
#' library(sp)
#' x=1:1000/100
#' l1 = Lines(Line(cbind(x, sin(x)) ), ID='a' )
#' sl = SpatialLines( list(l1) )
#' tol1=5;
#' tol2 =2
#' sl1 = sp.CutSptialLines(sl, tol1)
#' sl2 = sp.CutSptialLines(sl, tol2)
#' par(mfrow=c(1,2))
#' plot(sl1, col=1:length(sl1));title(paste0('Tol=', tol1))
#' plot(sl2, col=1:length(sl2));title(paste0('Tol=', tol2))
#'
#' data(sh)
#' riv=sh$riv
#' x = sp.CutSptialLines(riv, tol=5)
#' par(mfrow=c(2,1))
#' plot(riv, col=1:length(riv), lwd=3);
#'
#' plot(riv, col='gray', lwd=3);
#' plot(add=TRUE, x, col=1:length(x))
sp.CutSptialLines <- function(sl, tol){
  msg='sp.CutSptialLines::'
  ll = rgeos::gLength(sl, byid = TRUE)
  if(all(ll < tol) ){
    ret = sl
  }else{
    nsp = length(sl)
    xsl = list(); ik=1
    for(i in 1:nsp){
      sx = sl[i, ]
      pxy = extractCoords(sx,unique = TRUE)
      np = nrow(pxy)
      dacc = cumsum( sp::LineLength(pxy, sum = FALSE))
      # dacc =getDist(pxy)
      tol = max(c(tol, min(dacc) ) )
      len= rgeos::gLength(sx)
      if(len > tol){
        nsplit = ceiling(len / tol)
      }else{
        nsplit = 1
      }
      dd = len / nsplit
      v0 = 1  # Vetex 0, Vetex 1
      message(msg, i, '/', nsp, '\t', nsplit, '\t', round(dd, 2) )
      for(k in 1:nsplit){
        if(v0 >=np){
          break
        }
        # message(msg, '\t', k, '/', nsplit)
        dk = dd * k
        v1 = order(abs(dacc - dk), decreasing = FALSE)[1] + 1
        if(v1 + 1>np){
          v1 = np
        }
        message(msg, v0,'\t', v1)
        if(v0 == v1){
          next;
        }
        # plot(sl[i, ]);points(pxy); points(pxy[c(v0, v1), ], pch=2, col=2)
        xsl[[ik]]= sp::Lines(sp::Line( pxy[c(v0:v1), ]), ID=ik)
        ik=ik+1
        # points(pxy[v0:v1,], col=k)
        v0=v1
      }
    }
    nsl = length(xsl)
    tmp = sp::SpatialLines(xsl, proj4string = raster::crs(sl))
    ilen = rgeos::gLength(tmp, byid=TRUE)
    att=data.frame('INDEX'=1:length(tmp), 'Length'=ilen)
    ret = sp::SpatialLinesDataFrame(tmp, data = att)
  }
  return(ret)
}

#' Extract values on Raster map. The line is a straight line between (0,1). 
#' \code{extractRaster}
#' @param r Raster
#' @param xy coordinates of the line, dim=(Npoints, 2); x and y must be in [0, 1]
#' @param ext extension of value xy.
#' @param plot Whether plot result.
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors
#' @importFrom graphics grid hist lines par plot points
#' @importFrom methods as
#' @importFrom stats dist rnorm time
#' @importFrom utils read.table
#' @export
#' @examples
#' library(raster)
# r <- raster(ncol=36, nrow=18)
# r[] <- 1:ncell(r)
# extractRaster(r)
extractRaster<-function(r, xy=NULL, ext = raster::extent(r), plot=T){
  if(is.null(xy)){
    ndim = dim(r)
    x=0:ndim[2] / ndim[2]
    y = rep(0.5, length(x))
    xy = cbind(x,y)
  }
  x = ext[1] + xy[,1] * (ext[2]- ext[1] )
  y = ext[3] + xy[,2] * (ext[4]- ext[3] )
  if(plot){
    raster::plot(r);
    points(x, y, col=2)
    nx=length(x)
    points(x,y)
    graphics::arrows(x[1], y[1], x[nx], y[nx], lty=3, lwd=1.5, col=2)
    # lines(x,y, lwd=1.5, col=2, lty=2)
  }
  v = raster::extract(r, cbind(x,y))
  ret = cbind('x'=x,'y'=y,'z'=v)
  return(ret)
}

#' Simplify SpatialData.
#' @param x SpatialData
#' @return Simplified SpatialData
#' @export
SimpleSpatial <-function(x){
  # n1=length(x@polygons)
  # nj=unlist(lapply(1:n1, function(i){ length(x@polygons[[i]]@Polygons) } ))
  # x@polygons[[1]]@Polygons[[1]]@coords
  msg='SimpleSpatial'
  ni = length(x@polygons)
  k=1
  sl=list()
  for(i in 1:ni){
    nj = length(x@polygons[[1]]@Polygons)
    for(j in 1:nj){
      cd = x@polygons[[i]]@Polygons[[j]]@coords
      np=nrow(cd)
      message(msg, i,'-',j ,'\t', np)
      sl[[k]] = paste('POLYGON((', paste( paste(cd[,1], cd[,2]), collapse = ',' ), '))')
      if(k==1){
        str = sl[[k]]
      }else{
        str = paste(str, ',', sl[[k]] )
      }
      k=k+1
    }
  }
  r=rgeos::readWKT(paste('GEOMETRYCOLLECTION(', str, ')'))
}

#' Find the points in distance less than tol.
#' \code{PointInDistance}
#' @param pt 2-column coordinates (x,y).
#' @param tol Tolerance
#' @return Index of points that within tol
#' @export
PointInDistance <- function(pt, tol){
  msg='PointInDistance'
  dm = as.matrix(stats::dist(pt, diag  = TRUE))
  dm[dm==0]=NA
  # View(dm)
  dmin=apply(dm, 1, min, na.rm=T)
  id=which(dmin < 100)
  id1=id2=NULL
  tmp1=tmp=dmin[id]
  i=2
  for(i in 1:length(id)){
    if(id[i] %in% id2){next }
    id1 = c(id1, i); id1
    tmp1[id1] = NA; tmp1
    id[which(tmp1 %in% tmp[id1])]
    id2=unique(c(id2,  id[which(tmp1 %in% tmp[id1])]))
    id2
    tmp1=tmp
  }
  cbind('P1'=id[id1], P2=id2)
}

#' Conver the MULTIPOLYGONS to SINGLEPOLYGONS.
#' @param x the spatialpolygon*
#' @param id Index of the sorted (decreasing) polygons to return. default = 0;
#' @export
SinglePolygon <- function(x, id=0){
  n1 = length(x)
  y1 = list()
  y2 = list()
  k=1
  for(i in 1:n1){
    message('level 1:', i, '/', n1)
    x1 = x@polygons[[i]]
    
    n2=length(x1@Polygons)
    for(j in 1:n2){
      message('level 2:', j, '/', n2)
      x2 = x1@Polygons[[j]]
      y1[[k]] = Polygons( list(x2), ID=k)
      k=k+1
    }
  }
  
  y=sp::SpatialPolygonsDataFrame(SpatialPolygons(y1), data=data.frame('ID'=2:k-1))
  
  if(id < 1){
    return(y)
  }else{
    ia=rgeos::gArea(y, byid = TRUE)
    id=order(ia, decreasing = TRUE)[1]
    return(y[id,])
    
  }
}

#' Remove the duplicated lines, which means the FROM and TO points are identical.
#' \code{rmDuplicatedLines}
#' @param x ShapeLine*
#' @param ... More options in duplicated()
#' @return ShapeLine* without duplicated lines.
#' @export
rmDuplicatedLines <- function(x, ...){
  # x=spi.riv
  cd = extractCoords(x)
  # dim(cd)
  ft = FromToNode(x, cd)
  id=which(duplicated(ft[, -1], MARGIN = 1, ...))
  if(length(id)>0){
    r =x[-id,]
  }else{
    r = x
  }
  return(r)
}

#' Generate Thiesson Polygons from a point data
#' \code{voronoipolygons}
#' @param x ShapePoints* in PCS
#' @param pts Coordinates (x,y) of the points
#' @param rw extent
#' @param crs projection parameters
#' @return ShapePolygon*
#' @export
#' @examples 
#' library(rgeos)
#' library(rSHUD)
#' n=10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x=readWKT(str)
#' vx = voronoipolygons(x)
#' raster::plot(vx, axes=TRUE)
#' raster::plot(add=TRUE, x, col=2)
#' #' ====END====
#' 
#' x=1:5
#' y=1:5
#' xy=expand.grid(x,y)
#' vx=voronoipolygons(pts=xy)
#' plot(vx, axes=TRUE)
#' points(xy)
#' #' ====END====
#' 
#' library(rgdal)
#' library(raster)
#' library(rgeos)
#' n=10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x=readWKT(str)
#' y=readWKT(paste('MULTIPOINT(', paste(paste('(', xx+2, yy+2, ')'), collapse = ','), ')'))
#' e1 = extent(y)
#' e2 = extent(x)
#' rw = c(min(e1[1], e2[1]),
#'        max(e1[2], e2[2]),
#'        min(e1[3], e2[3]),
#'        max(e1[4], e2[4]) ) + c(-1, 1, -1, 1)
#' vx=voronoipolygons(x=x, rw=rw)
#' plot(vx); plot(add=TRUE, x, col=2); plot(add=TRUE, y, col=3)
voronoipolygons = function(x, pts = x@coords, rw=NULL, crs=NULL) {
  z = deldir::deldir(pts[,1], pts[,2], rw=rw)
  w = deldir::tile.list(z)
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = sp::Polygons(list(sp::Polygon(pcrds)), ID=as.character(i))
  }
  SP = sp::SpatialPolygons(polys)
  voronoi = sp::SpatialPolygonsDataFrame(SP, 
                                         data=data.frame(x=pts[,1], 
                                                         y=pts[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                     function(x) slot(x, 'ID'))))
  if(!is.null(crs)){
    raster::crs(voronoi) = crs;
  }
  return(voronoi)
}



#' Generate the coverage map for forcing sites.
#' \code{ForcingCoverage}
#' @param sp.meteoSite ShapePoints* in PCS
#' @param pcs Projected Coordinate System
#' @param gcs Geographic Coordinate System
#' @param dem DEM raster
#' @param wbd watershed boundary
#' @param enlarge enlarge factor for the boundary.
#' @return ShapePolygon*
#' @export
ForcingCoverage <- function(
  sp.meteoSite = NULL,
  filenames = NULL,
  pcs,
  gcs = 4326,
  dem,
  wbd,
  enlarge = 10000
) {
  msg <- "ForcingCoverage::"

  if (missing(pcs) || is.null(pcs)) stop(msg, "pcs is required")
  if (missing(wbd) || is.null(wbd)) stop(msg, "wbd is required")

  pcs_crs <- .rshud_as_crs(pcs)
  gcs_crs <- .rshud_as_crs(gcs)

  wbd_sf <- .rshud_as_sf(wbd)
  wbd_pcs <- sf::st_transform(wbd_sf, pcs_crs)

  sites_pcs <- NULL
  if (is.null(sp.meteoSite)) {
    sites_pcs <- sf::st_sf(geometry = sf::st_centroid(sf::st_union(wbd_pcs)))
  } else {
    sites_sf <- .rshud_as_sf(sp.meteoSite)
    sites_pcs <- sf::st_transform(sites_sf, pcs_crs)
  }
  sites_pcs <- sf::st_cast(sites_pcs, "POINT", warn = FALSE)

  sites_gcs <- sf::st_transform(sites_pcs, gcs_crs)
  ll <- sf::st_coordinates(sites_gcs)
  xy <- sf::st_coordinates(sites_pcs)

  n <- nrow(xy)
  if (n == 0) stop(msg, "No forcing sites")

  z <- rep(NA_real_, n)
  if (!missing(dem) && !is.null(dem)) {
    r <- terra::rast(dem)
    pts_v <- terra::vect(sites_pcs)
    ex <- terra::extract(r, pts_v)
    if (!is.null(ex) && ncol(ex) >= 2) z <- as.numeric(ex[[2]])
  }

  if (is.null(filenames)) {
    id_for_files <- NULL
    if (!is.null(sp.meteoSite)) {
      if (inherits(sp.meteoSite, "sf") && "ID" %in% names(sp.meteoSite)) {
        id_for_files <- sp.meteoSite$ID
      } else if (inherits(sp.meteoSite, "SpatVector")) {
        vals <- tryCatch(terra::values(sp.meteoSite, dataframe = TRUE), error = function(e) NULL)
        if (!is.null(vals) && "ID" %in% names(vals)) id_for_files <- vals$ID
      } else if (any(grepl("^Spatial", class(sp.meteoSite)))) {
        if (requireNamespace("sp", quietly = TRUE)) {
          id_for_files <- tryCatch(sp.meteoSite@data$ID, error = function(e) NULL)
        }
      }
    }
    if (is.null(id_for_files) || length(id_for_files) != n) id_for_files <- seq_len(n)
    filenames <- paste0(id_for_files, ".csv")
  }
  if (length(filenames) != n) stop(msg, "filenames must have length ", n)

  att <- data.frame(
    ID = seq_len(n),
    Lon = ll[, 1],
    Lat = ll[, 2],
    X = xy[, 1],
    Y = xy[, 2],
    Z = z,
    Filename = filenames
  )
  att[is.na(att)] <- -9999

  bbox_w <- sf::st_bbox(wbd_pcs)
  bbox_s <- sf::st_bbox(sites_pcs)
  rw <- c(
    min(bbox_w[["xmin"]], bbox_s[["xmin"]]),
    max(bbox_w[["xmax"]], bbox_s[["xmax"]]),
    min(bbox_w[["ymin"]], bbox_s[["ymin"]]),
    max(bbox_w[["ymax"]], bbox_s[["ymax"]])
  )

  if (is.null(enlarge)) {
    enlarge <- min(diff(rw[1:2]), diff(rw[3:4])) * 0.02
  }
  rw <- rw + c(-1, 1, -1, 1) * enlarge

  cover <- NULL
  if (n < 2) {
    cover <- fishnet(xx = rw[1:2], yy = rw[3:4], crs = pcs_crs, type = "polygon")
  } else {
    pts_v <- terra::vect(sites_pcs)
    terra::values(pts_v) <- data.frame(site_index = seq_len(n))
    vor <- terra::voronoi(pts_v, bnd = terra::ext(rw[1], rw[2], rw[3], rw[4]))
    cover <- sf::st_as_sf(vor)
    cover <- cover[order(cover$site_index), , drop = FALSE]
  }

  if (nrow(cover) != n) stop(msg, "Coverage polygons count mismatch: ", nrow(cover), " vs ", n)
  out <- sf::st_sf(att, geometry = sf::st_geometry(cover))
  sf::st_crs(out) <- pcs_crs
  return(out)
}

#' Find the subset of a extent in a grid.
#' \code{grid.subset}
#' @param ext txtent of the grid
#' @param res resolution of the grid
#' @param ext.sub the extent of subset
#' @param x x coordinates of the grids
#' @param y y coordinates of the grids
#' @return list, list(xid, yid, x, y)
#' @export
#' @examples 
grid.subset <- function(ext, res,
                        ext.sub,
                        dx =matrix(res, 2,1)[1],
                        dy =matrix(res, 2,1)[2],
                        x = seq(ext[1]+dx/2, ext[2]-dx/2, by=dx),
                        y = seq(ext[3]+dy/2, ext[4]-dy/2, by=dy)
                        ){
  xmin = min(x - dx/2); xmax = max(x + dx/2)
  ymin = min(y - dy/2); ymax = max(y + dy/2)
  # ext= c(min(x), max(x), min(y), max(y))

  if(ext.sub[1] < xmin | ext.sub[2] > xmax | ext.sub[3] < ymin | ext.sub[4] > ymax){
    warning(paste('Extent required is larger than the boundbox of dataset'))
    message(paste(ext.sub, collaps=','))
    message(paste(c(xmin,xmax,ymin, ymax), collaps=','))
  }
  xid = min(which(abs(x  - ext.sub[1]) <= dx/2)):max(which(abs(x  - ext.sub[2]) <= dx/2))
  yid = min(which(abs(y  - ext.sub[3]) <= dy/2)):max(which(abs(y  - ext.sub[4]) <= dy/2))
  nx = length(xid); ny = length(yid)
  x.cord = x[xid]; y.cord = y[yid]
  rt = list('xid' = xid,
            'yid' = yid,
            'x'=x.cord,
            'y'=y.cord)
  return(rt)
}


#' Generate a shapefile from coordinates.
#' \code{xy2shp}
#' @param xy matrix
#' @param df attribute table
#' @param crs projection parameters
#' @param shape Shape of the result in points, lines or polygons
#' @return SpatialPointsDataFrame, SpatialLinesDataFrame, or SpatialPolygonsDataFrame
#' @export
#' @examples 
#' library(raster)
#' xy=list(cbind(c(0, 2, 1), c(0, 0, 2)),  cbind(c(0, 2, 1), c(0, 0, 2))+2)
#' sp1 = xy2shp(xy=xy, shape = 'polygon')
#' raster::plot(sp1, axes=TRUE, col='gray')
#' 
#' sp2 = xy2shp(xy=xy, shape = 'lines')
#' raster::plot(sp2, add=TRUE, lty=2, lwd=3,col='red')
#' sp3 = xy2shp(xy=xy, shape = 'POINTS')
#' raster::plot(sp3, add=TRUE, pch=1, cex=2)
#' 
xy2shp <- function(xy, df=NULL, crs=NULL, shape='points'){
  shape0 <- tolower(shape)
  crs_sf <- if (is.null(crs)) sf::NA_crs_ else .rshud_as_crs(crs)

  xy_list <- if (is.list(xy)) xy else list(xy)

  if (grepl("point", shape0)) {
    coords <- do.call(rbind, xy_list)
    coords <- as.matrix(coords)
    if (ncol(coords) != 2) stop("xy2shp: points require 2 columns (x,y)")

    if (is.null(df)) df <- data.frame(ID = seq_len(nrow(coords)))
    if (nrow(df) != nrow(coords)) stop("xy2shp: df row count must match points")

    geom <- sf::st_sfc(
      lapply(seq_len(nrow(coords)), function(i) sf::st_point(as.numeric(coords[i, 1:2]))),
      crs = crs_sf
    )
    return(sf::st_sf(df, geometry = geom))
  }

  if (grepl("line", shape0)) {
    if (is.null(df)) df <- data.frame(ID = seq_len(length(xy_list)))
    if (nrow(df) != length(xy_list)) stop("xy2shp: df row count must match lines")

    geom <- sf::st_sfc(
      lapply(xy_list, function(m) {
        m <- as.matrix(m)
        if (ncol(m) != 2) stop("xy2shp: lines require 2 columns (x,y)")
        sf::st_linestring(m[, 1:2, drop = FALSE])
      }),
      crs = crs_sf
    )
    return(sf::st_sf(df, geometry = geom))
  }

  if (grepl("polygon", shape0)) {
    if (is.null(df)) df <- data.frame(ID = seq_len(length(xy_list)))
    if (nrow(df) != length(xy_list)) stop("xy2shp: df row count must match polygons")

    geom <- sf::st_sfc(
      lapply(xy_list, function(m) {
        m <- as.matrix(m)
        if (ncol(m) != 2) stop("xy2shp: polygons require 2 columns (x,y)")
        ring <- m[, 1:2, drop = FALSE]
        if (nrow(ring) < 3) stop("xy2shp: polygon ring must have >= 3 points")
        if (!all(ring[1, ] == ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
        sf::st_polygon(list(ring))
      }),
      crs = crs_sf
    )
    return(sf::st_sf(df, geometry = geom))
  }

  stop("xy2shp: shape must be points|lines|polygon")
}
