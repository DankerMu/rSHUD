.rshud_as_crs <- function(crs) {
  if (is.null(crs)) return(sf::NA_crs_)
  if (inherits(crs, "crs")) return(crs)
  if (is.numeric(crs) && length(crs) == 1 && is.finite(crs)) return(sf::st_crs(crs))
  if (is.character(crs) && length(crs) == 1 && nzchar(crs)) return(sf::st_crs(crs))
  if (inherits(crs, "CRS")) {
    if (!requireNamespace("sp", quietly = TRUE)) {
      stop("sp is required to convert sp::CRS -> sf::crs")
    }
    return(sf::st_crs(sp::proj4string(crs)))
  }
  stop("Unsupported CRS format: ", paste(class(crs), collapse = "/"))
}

.rshud_crs_wkt <- function(crs) {
  crs_sf <- .rshud_as_crs(crs)
  if (is.na(crs_sf)) return(NA_character_)
  if (!is.null(crs_sf$wkt) && nzchar(crs_sf$wkt)) return(crs_sf$wkt)
  if (!is.null(crs_sf$proj4string) && nzchar(crs_sf$proj4string)) return(crs_sf$proj4string)
  NA_character_
}

.rshud_as_sf <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "sf")) return(x)
  if (inherits(x, "sfc")) return(sf::st_sf(geometry = x))
  if (inherits(x, "SpatVector")) return(sf::st_as_sf(x))
  if (any(grepl("^Spatial", class(x)))) return(sf::st_as_sf(x))
  stop("Unsupported spatial object: ", paste(class(x), collapse = "/"))
}

.rshud_as_spatvector <- function(x, crs = NULL) {
  if (inherits(x, "SpatVector")) {
    v <- x
  } else if (inherits(x, "sfc")) {
    v <- terra::vect(sf::st_sf(geometry = x))
  } else {
    v <- terra::vect(x)
  }
  if (!is.null(crs)) {
    wkt <- .rshud_crs_wkt(crs)
    if (!is.na(wkt)) terra::crs(v) <- wkt
  }
  v
}

.rshud_geom_drop_holes <- function(g) {
  gt <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
  if (length(gt) == 0) return(g)
  gt <- gt[[1]]

  if (gt == "POLYGON") {
    if (length(g) == 0) return(g)
    return(sf::st_polygon(list(g[[1]])))
  }
  if (gt == "MULTIPOLYGON") {
    if (length(g) == 0) return(g)
    polys <- lapply(g, function(poly) {
      if (length(poly) == 0) return(poly)
      list(poly[[1]])
    })
    return(sf::st_multipolygon(polys))
  }
  g
}

