test_that("snap_coords_to_ref restores exact coord_ref values", {
  skip_if_not_installed("sp")

  pts <- matrix(
    c(
      1000000.123, 2000000.456,
      1000001.234, 2000001.567,
      1000002.345, 2000002.678,
      1000003.456, 2000003.789
    ),
    ncol = 2,
    byrow = TRUE
  )

  sl <- sp::SpatialLines(list(sp::Lines(list(sp::Line(pts)), ID = "1")))
  coord_ref <- extractCoords(sl, unique = TRUE)

  sl_drift <- sl
  eps <- 1e-9
  sl_drift@lines[[1]]@Lines[[1]]@coords[, 1:2] <-
    sl_drift@lines[[1]]@Lines[[1]]@coords[, 1:2] + eps

  tol <- 1e-6
  snapped <- rSHUD:::snap_coords_to_ref(sl_drift, coord_ref = coord_ref, tol = tol)
  coords_snapped <- snapped@lines[[1]]@Lines[[1]]@coords[, 1:2, drop = FALSE]

  expect_equal(coords_snapped, pts)
})

test_that("FromToNode simplify returns non-zero Fr/To nodes", {
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")
  skip_if_not_installed("rgeos")

  make_line <- function(from, to, n, noise_sd = 0.01) {
    xy <- cbind(
      seq(from[1], to[1], length.out = n),
      seq(from[2], to[2], length.out = n)
    )
    xy <- xy + matrix(stats::rnorm(n * 2, sd = noise_sd), ncol = 2)
    xy[1, ] <- from
    xy[n, ] <- to
    xy
  }

  set.seed(1)
  p0 <- c(1000000.123, 2000000.456)
  p1 <- c(1000100.789, 1999900.321)
  p2 <- c(1000200.321, 1999800.987)
  p3 <- c(1000050.555, 1999850.111)

  seg1 <- make_line(p0, p1, n = 50)
  seg2 <- make_line(p1, p2, n = 60)
  seg3 <- make_line(p1, p3, n = 40)

  sl <- sp::SpatialLines(list(
    sp::Lines(list(sp::Line(seg1)), ID = "1"),
    sp::Lines(list(sp::Line(seg2)), ID = "2"),
    sp::Lines(list(sp::Line(seg3)), ID = "3")
  ))

  coord_ref <- extractCoords(sl, unique = TRUE)
  ft <- FromToNode(sl, coord = coord_ref, simplify = TRUE)

  expect_true(all(ft[, "FrNode"] > 0))
  expect_true(all(ft[, "ToNode"] > 0))
  expect_true(all(ft[, "FrNode"] <= nrow(coord_ref)))
  expect_true(all(ft[, "ToNode"] <= nrow(coord_ref)))

  pt.list <- unlist(sp::coordinates(sl), recursive = FALSE)
  expect_equal(length(pt.list), length(sl))

  fr_true <- do.call(rbind, lapply(pt.list, function(x) x[1, 1:2, drop = FALSE]))
  to_true <- do.call(rbind, lapply(pt.list, function(x) x[nrow(x), 1:2, drop = FALSE]))
  fr <- coord_ref[ft[, "FrNode"], , drop = FALSE]
  to <- coord_ref[ft[, "ToNode"], , drop = FALSE]

  dimnames(fr) <- NULL
  dimnames(fr_true) <- NULL
  dimnames(to) <- NULL
  dimnames(to_true) <- NULL

  expect_equal(fr, fr_true)
  expect_equal(to, to_true)
})
