test_that("fishnet polygon returns sf", {
  xx <- seq(0, 2, by = 1)
  yy <- seq(0, 2, by = 1)
  fn <- fishnet(xx = xx, yy = yy, crs = 4326, type = "polygon")
  expect_s3_class(fn, "sf")
  expect_equal(nrow(fn), 4)
  expect_true(all(c("xmin", "xmax", "ymin", "ymax", "xcenter", "ycenter") %in% names(fn)))
  expect_equal(sf::st_crs(fn)$epsg, 4326)
  expect_true(all(sf::st_geometry_type(fn) == "POLYGON"))
})

test_that("fishnet points returns sf", {
  xx <- 1:3
  yy <- 1:2
  fn <- fishnet(xx = xx, yy = yy, crs = 4326, type = "point")
  expect_s3_class(fn, "sf")
  expect_equal(nrow(fn), length(xx) * length(yy))
  expect_true(all(sf::st_geometry_type(fn) == "POINT"))
})

test_that("xy2shp creates sf geometries", {
  pts <- matrix(c(0, 0, 1, 2), ncol = 2, byrow = TRUE)
  s <- xy2shp(pts, df = data.frame(ID = 1:2), crs = 4326, shape = "points")
  expect_s3_class(s, "sf")
  expect_equal(nrow(s), 2)
  expect_true(all(sf::st_geometry_type(s) == "POINT"))

  ln <- list(matrix(c(0, 0, 1, 1, 2, 0), ncol = 2, byrow = TRUE))
  s2 <- xy2shp(ln, df = data.frame(ID = 1), crs = 4326, shape = "lines")
  expect_s3_class(s2, "sf")
  expect_true(all(sf::st_geometry_type(s2) == "LINESTRING"))
})

test_that("removeholes drops interior rings for sf polygons", {
  outer <- matrix(c(0, 0, 0, 10, 10, 10, 10, 0, 0, 0), ncol = 2, byrow = TRUE)
  hole <- matrix(c(3, 3, 7, 3, 7, 7, 3, 7, 3, 3), ncol = 2, byrow = TRUE)
  p <- sf::st_polygon(list(outer, hole))
  x <- sf::st_sf(geometry = sf::st_sfc(p, crs = 4326))
  y <- removeholes(x)
  expect_s3_class(y, "sf")
  expect_equal(length(sf::st_geometry(y)[[1]]), 1)
})

test_that("ForcingCoverage returns one polygon per site", {
  wbd_ring <- matrix(c(0, 0, 0, 10, 10, 10, 10, 0, 0, 0), ncol = 2, byrow = TRUE)
  wbd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(wbd_ring)), crs = 4326))

  sites <- sf::st_as_sf(
    data.frame(ID = c(1, 2), x = c(2, 8), y = c(2, 8)),
    coords = c("x", "y"),
    crs = 4326
  )

  fc <- ForcingCoverage(sp.meteoSite = sites, pcs = 4326, dem = NULL, wbd = wbd, enlarge = 0)
  expect_s3_class(fc, "sf")
  expect_equal(nrow(fc), 2)
  expect_true(all(c("ID", "Lon", "Lat", "X", "Y", "Z", "Filename") %in% names(fc)))
  expect_true(all(sf::st_geometry_type(fc) == "POLYGON"))
})

test_that("xyz2Raster returns terra SpatRaster", {
  x <- c(0, 1, 2)
  y <- c(10, 11)
  arr <- matrix(1:6, nrow = 3, ncol = 2)
  r <- xyz2Raster(x = x, y = y, arr = arr, flip = TRUE, plot = FALSE)
  expect_true(inherits(r, "SpatRaster"))
  expect_equal(terra::ncol(r), length(x))
  expect_equal(terra::nrow(r), length(y))

  arr3 <- array(1:(3 * 2 * 2), dim = c(3, 2, 2))
  r2 <- xyz2Raster(x = x, y = y, arr = arr3, flip = TRUE, plot = FALSE)
  expect_equal(terra::nlyr(r2), 2)
})

test_that("sp.Tri2Shape does not require rgeos", {
  tri <- list(
    T = matrix(c(1, 2, 3), nrow = 1),
    P = matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  )
  shp <- sp.Tri2Shape(tri, dbf = data.frame(ID = 1), crs = NULL)
  expect_true(inherits(shp, "SpatialPolygonsDataFrame"))
  expect_equal(nrow(shp), 1)
  expect_true("Area" %in% names(shp@data))
  expect_equal(as.numeric(shp@data$Area), 0.5, tolerance = 1e-8)
})
