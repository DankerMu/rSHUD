#' determine index of downstream
#' \code{sp.RiverOrder}
#' @param sp SpatialLines*
#' @param coord coordinates of \code{sp}.
#' @return Index of downstream for each segements
#' @export
sp.RiverDown <- function(sp, coord = extractCoords(sp)){
  msg = 'sp.RiverDown::'
  ext = raster::extent(sp)
  sp = rgeos::gSimplify(sp, tol = (ext[2] - ext[1] )*0.01)
  ft = rbind(FromToNode(sp, coord = coord)[, 2:3])
  nsp = length(sp)
  idown = rep(-3, nsp)
  for(i in 1:nsp){
    pto = ft[i,2]
    id = which(ft[,1] == pto)
    nid = length(id)
    if(nid == 1){
      idown[i] = id
    }else if(length(id) > 1){
      message(msg, i, '/', nsp, '\t', nid, ' Downstream found')
      print(id)
      idown[i] = id[1]
      # xid = c(id, i)
      # plot(sp[xid, ]); points(coord[ft[i, ], ], col=1:2);
      # points(coord[ft[id[1],], ], col=1:2, pch=2)
      # points(coord[ft[id[2],], ], col=1:2, pch=3)
    }else{
      # VOID
    }
  }
  idown
}
# sp.slt = spr
# rivdown = sp.RiverDown(sp.slt)
#' determine River path, disolve the river network based on downstream-relationship.
#' \code{sp.RiverPath}
#' @param sp SpatialLines*
#' @param idown downstream index of \code{sp}.
#' @param coord coordinates of \code{sp}.
#' @param tol.simplify tolerance for simplification of river lines.
#' @return List of River Path(SegIDs, PointIDs, sp)
#' @export
#' @examples 
# data(sh)
# riv=sh$riv
# x=sp.CutSptialLines(riv, tol=200)
# px(x)
# xx=sp.RiverPath(x)
# px(xx$sp)
sp.RiverPath <- function(sp, 
                         coord = extractCoords(sp, unique = TRUE), 
                         idown = sp.RiverDown(sp, coord = coord),
                         tol.simplify = 0){
  msg = 'sp.RiverPath::'
  ext = raster::extent(sp)
  if(tol.simplify > 0){
    sp = rgeos::gSimplify(sp, tol = tol.simplify)
  }
  pt.list <- unlist(coordinates(sp), recursive = FALSE)
  id.list<-lapply(pt.list, function(x) { xy2ID(xy = x, coord = coord)})

  ft = FromToNode(sp, coord = coord)[, 2:3]
  nsp = length(sp)
  p.frq = data.frame(table(unlist(id.list)) )
  p0 = ft[!(ft[, 1] %in% ft[,2]), 1]
  tmp =  p.frq[p.frq[,2]>2, 1] # more than twice.
  p1 = as.numeric(levels(tmp))[tmp]
  p2 = ft[!(ft[, 2] %in% ft[,1]), 2] # TO that not in FROM
  
  # jid = which(table(as.matrix(ft)) > 2) #id of joint points
  riv.keyid = c(which(idown<0), which(ft[,2] %in% p1) )
  updown = cbind(1:nsp, idown)
  goUp <- function(updown, id0){
    uid = which(idown == id0) # id of my upstream
    if(length(uid) == 1){ #up stream exist
      ret = c(goUp(updown, uid), id0)
    }else{
      ret = id0
    }
    ret
  }
  message(msg, 'to Upstream')
  StreamPath = lapply(riv.keyid, function(x) goUp(updown, x))
  #========Point ID==========
  nstr = length(StreamPath)
  pl = list()
  message(msg, 'searching from/to nodes')
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = unique(unlist(lapply(1:length(sid), function(x) id.list[[sid[x]]])))
    # print(sid)
    # print(pid)
    pl[[i]] = pid
  }
  #=========Spatial Lines===============
  ll = list()
  message(msg, 'build new spatialLines')
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = unique(unlist(lapply(1:length(sid), function(x) id.list[[sid[x]]])))
    ll[[i]] = sp::Lines( sp::Line(coord[pid, ] ), ID = i)
  }
  spx = sp::SpatialLines(ll, proj4string = raster::crs(sp))
  ret <- list(SegIDs = StreamPath,
              PointIDs= pl,
              sp = spx)
}


#' calculate river order
#' \code{sp.RiverOrder}
#' @param sp SpatialLines*
#' @param coord coordinates of \code{sp}.
#' @return Stream Order of SpatialLines*
#' @export
#' @examples 
#' library(rgdal)
#' data(sac)
#' riv=sac$riv
#' ord = sp.RiverOrder(riv)
#' print(ord)
#' idx = sort(unique(ord))
#' raster::plot(riv, col=ord)
#' legend('topleft', legend=idx, col=idx, lwd=1, lty=1)
sp.RiverOrder <- function(sp, coord = extractCoords(sp)){
  msg='sp.RiverOrder::'
  ext = raster::extent(sp)
  sp = rgeos::gSimplify(sp, tol = (ext[2] - ext[1] )*0.01)
  get1st <- function(x){
    fr = x[,2]
    to = x[,3]
    tb= table(x[,2:3])
    pid = as.numeric(names(tb))
    ncount = as.numeric(tb)
    p.out = to[which(! to %in% fr)]
    p.jnt = pid[which(ncount> 2)]
    p.key = c(p.out, p.jnt)
    # print(p.key)
    tid=fr[ !fr %in% to] #target from-id
    # plot(riv.simp)
    # points(coord[ ])
    # points(coord[p.key,], col=2)
    nd = length(tid)
    ret = NULL
    # message(msg, 'Number of Streams = ', nd)
    for(i in 1:nd){
      cr = tid[i]
      for(j in 1:100000){
        # message(msg, 'Stream ',i,'\tSeg ', j, '\tFrom Node ', cr)
        rid = which(fr %in% cr)
        ret=c(ret, rid)
        sid= to[rid] #ID of to-point;
        # if(length(sid) < 1 | length(sid)>1){
        #   print(sid)
        # }
        if( any(sid  %in% p.key) ){ #if sid IS IN key-points(outlets or joint point).
          break;
        }else{
          cr = sid;
        }
      }
    }
    ret;
  }

  # ext=raster::extent(sp)
  # sp.tmp =  rgeos::gSimplify(sp, tol=max(diff(ext[1:2] ), diff(ext[3:4])) )
  ft0 = FromToNode(sp, coord)
  ft = rbind(unique(ft0[, 2:3]))
  if(length(sp) != nrow(ft)){
    message(msg, 'ERROR: duplicated river reches extis.')
    stop('STOP WITH ERROR')
  }
  x = cbind(1:length(sp), ft)
  y =x
  x.ord =x[,1]*0
  for(i in 1:10000){
    ids = y[,1]
    message(msg, 'Order = ', i)
    id = get1st(y)
    x.ord[ ids[id] ] = i
    y=rbind(y[-1 * id ,])
    ny = length(y)
    # message(msg, 'Left Segments =', ny)
    if(ny<=0){ break }
  }
  x.ord
}
#' return the Index of nodes in SpatialData
#' \code{NodeIDList}
#' @param sp SpatialLines*
#' @param coord Coordinate of vertex in \code{sp}
#' @return list of Index of each SpatialData, line or polygon
#' @export
NodeIDList <- function(sp, coord = extractCoords(sp, unique = TRUE) ){
  pt.list <- unlist(sp::coordinates(sp), recursive = FALSE)
  id.list<-lapply(pt.list, function(x) { xy2ID(x, coord = coord)})
}

snap_coords_to_ref <- function(sp_simplified, coord_ref, tol) {
  if (!inherits(sp_simplified, "SpatialLines")) {
    stop("snap_coords_to_ref: `sp_simplified` must be SpatialLines*.", call. = FALSE)
  }
  if (is.data.frame(coord_ref)) {
    coord_ref <- as.matrix(coord_ref)
  }
  if (!is.matrix(coord_ref) || ncol(coord_ref) < 2) {
    stop("snap_coords_to_ref: `coord_ref` must be a matrix with 2 columns.", call. = FALSE)
  }
  coord_ref <- coord_ref[, 1:2, drop = FALSE]
  tol <- as.numeric(tol)
  if (length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    return(sp_simplified)
  }

  lines <- sp_simplified@lines
  n_parts <- sum(vapply(lines, function(x) length(x@Lines), integer(1)))
  if (n_parts < 1) {
    return(sp_simplified)
  }

  coords_list <- vector("list", n_parts)
  line_i <- integer(n_parts)
  line_j <- integer(n_parts)
  k <- 0L
  for (i in seq_along(lines)) {
    parts <- lines[[i]]@Lines
    for (j in seq_along(parts)) {
      k <- k + 1L
      coords_list[[k]] <- parts[[j]]@coords
      line_i[[k]] <- i
      line_j[[k]] <- j
    }
  }

  row_counts <- vapply(coords_list, nrow, integer(1))
  total_rows <- sum(row_counts)
  if (total_rows < 1) {
    return(sp_simplified)
  }

  coord_simplified <- do.call(rbind, lapply(coords_list, function(x) x[, 1:2, drop = FALSE]))

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = coord_ref, query = coord_simplified, k = 1)
    idx <- nn$nn.idx[, 1]
    d <- nn$nn.dists[, 1]
    keep <- is.finite(d) & (d <= tol)
    if (any(keep)) {
      coord_simplified[keep, ] <- coord_ref[idx[keep], , drop = FALSE]
    }
  } else {
    tol2 <- tol * tol

    cx_ref <- floor(coord_ref[, 1] / tol)
    cy_ref <- floor(coord_ref[, 2] / tol)
    key_ref <- paste(cx_ref, cy_ref, sep = ":")
    buckets <- split(seq_len(nrow(coord_ref)), key_ref)

    cx_q <- floor(coord_simplified[, 1] / tol)
    cy_q <- floor(coord_simplified[, 2] / tol)
    key_q <- paste(cx_q, cy_q, sep = ":")
    q_groups <- split(seq_len(nrow(coord_simplified)), key_q)

    offsets <- c(-1L, 0L, 1L)
    for (key in names(q_groups)) {
      idx_q <- q_groups[[key]]
      if (length(idx_q) < 1) {
        next
      }
      cx <- cx_q[idx_q[1]]
      cy <- cy_q[idx_q[1]]
      neighbor_keys <- paste(
        rep(cx + offsets, each = 3),
        rep(cy + offsets, times = 3),
        sep = ":"
      )
      cand_idx <- unlist(buckets[neighbor_keys], use.names = FALSE)
      if (length(cand_idx) < 1) {
        next
      }

      ref_cand <- coord_ref[cand_idx, , drop = FALSE]
      q_mat <- coord_simplified[idx_q, , drop = FALSE]

      a2 <- rowSums(q_mat * q_mat)
      b2 <- rowSums(ref_cand * ref_cand)
      d2 <- outer(a2, b2, "+") - 2 * (q_mat %*% t(ref_cand))
      d2 <- pmax(d2, 0)

      best_pos <- max.col(-d2, ties.method = "first")
      best_d2 <- d2[cbind(seq_len(nrow(q_mat)), best_pos)]
      keep <- is.finite(best_d2) & (best_d2 <= tol2)
      if (any(keep)) {
        coord_simplified[idx_q[keep], ] <- ref_cand[best_pos[keep], , drop = FALSE]
      }
    }
  }

  ends <- cumsum(row_counts)
  starts <- ends - row_counts + 1L
  for (k in seq_len(n_parts)) {
    i <- line_i[[k]]
    j <- line_j[[k]]
    idx <- starts[k]:ends[k]
    coords <- lines[[i]]@Lines[[j]]@coords
    coords[, 1:2] <- coord_simplified[idx, , drop = FALSE]
    lines[[i]]@Lines[[j]]@coords <- coords
  }
  sp_simplified@lines <- lines

  bbox <- matrix(
    c(
      min(coord_simplified[, 1]),
      max(coord_simplified[, 1]),
      min(coord_simplified[, 2]),
      max(coord_simplified[, 2])
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("x", "y"), c("min", "max"))
  )
  sp_simplified@bbox <- bbox
  sp_simplified
}
#' return the FROM and TO nodes index of the SpatialLines
#' \code{FromToNode}
#' @param sp SpatialLines*
#' @param coord Coordinate of vertex in \code{sp}
#' @param simplify whether simplify the spatialLines
#' @return FROM and TO nodes index of the SpatialLines
#' @export
FromToNode <- function(sp, coord = extractCoords(sp, unique = TRUE), simplify=TRUE){
  if(simplify){
    ext = raster::extent(sp)
    tol = (ext[2] - ext[1]) * 0.01
    coord <- force(coord)  # force eval before simplify reassigns sp
    sp = rgeos::gSimplify(sp, tol = tol)
    sp = snap_coords_to_ref(sp, coord, tol = tol)
  }
  id.list = NodeIDList(sp, coord=coord)
  frto = cbind(
    unlist(lapply(id.list, function(x) x[1])),
    unlist(lapply(id.list, function(x) x[length(x)] ) )
  )
  frto = cbind(1:length(sp), frto)
  colnames(frto)=c('ID', 'FrNode', 'ToNode')
  ret = rbind(frto)
  return(ret)
}

#' parameters for river types
#' \code{RiverType}
#' @param n number of types
#' @param width Width of rivers
#' @return data.frame of parameters
#' @export
RiverType <- function(n, width = 2 * (1:n) ){
  cn = c('Index', 'Depth', 'BankSlope',
         'Width', 'Sinuosity', 'Manning',
         'Cwr', 'KsatH', 'BedThick')
  nc = length(cn)
  rtype = cbind(1:n,
                5.0 + 0.5 * 1:n, #Depth
                0, #bankslope
                width, #Width
                1.0, #Sinuosity
                0.04, # manning's n, s/m^(1/3)
                # 4.63e-07, #manning's n, day/m^(1/3)
                0.6, #CWR
                0.1, #KsatH
                .1 # Bed Sediment Thickness.
  )
  colnames(rtype) = cn
  rtype
}
#' Cut the river by shud.mesh
#' \code{sp.RiverSeg}
#' @param sp.mesh shud mesh
#' @param sp.riv river
#' @return SpatialLinesDataFrame
#' @export
sp.RiverSeg <- function(sp.mesh, sp.riv){
  sp = sp::SpatialPolygonsDataFrame(sp.mesh, 
                                 data = data.frame('iEle'=1:length(sp.mesh) ),
                                 match.ID = FALSE)
  sr =sp::SpatialLinesDataFrame(sp.riv, 
                                data = data.frame('iRiv' = 1:length(sp.riv)),
                                match.ID = FALSE )
  seg = raster::intersect(sr, sp)
  seg
}

#' Find Shared points in SpatialLine*
#' \code{SharedPoints}
#' @param sp SpatialLines
#' @return Topologic relationship between components of SpatialLines.
#' @export
SharedPoints <- function(sp){
  msg='SharedPoints::'
#  sp=riv
  pl = extractCoords(sp, aslist = TRUE)
  nsp = length(sp)
  ret = matrix(0, nrow = nsp, ncol=nsp)
  for(i in 1:(nsp-1) ){
    for(j in (i+1):nsp ){
      plot(rbind(pl[[i]], pl[[j]]), asp=1); 
      points(pl[[i]], col=2)
      points(pl[[j]], col=3)
      cp = CommonPoints(pl[[i]], pl[[j]]) 
      if( is.null(cp) ){
        
      }else{
        message(msg, i, '/', j)
        # print(ret)
        # print(cp)
        # readline("go")
        ret[i,j] = cp[,1]
        ret[j,i] = cp[,2]
      }
    }
  }
  ret
}
