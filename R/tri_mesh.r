#' @importFrom utils head
path2seg <- function(x) {
  head(suppressWarnings(matrix(x, nrow = length(x) + 1, ncol = 2, byrow = FALSE)), -2L)
}


#' Generate triangle mesh
#'
#' Create triangle mesh structures from various inputs.
#'
#' Methods exist for SpatialPolygons, ...
#' @param x input data
#' @param ... arguments passed to methods
#'
#' @return a list of tibble data frames, using the gris-map_table model
#' @export
tri_mesh <- function(x, ...) {
  UseMethod("tri_mesh")
}

#' @rdname tri_mesh
#' @export
#' @importFrom sp over SpatialPoints proj4string CRS
#' @importFrom dplyr inner_join
#' @importFrom RTriangle pslg triangulate
#' @importFrom sp geometry
#' @importFrom spbabel map_table
#' @importFrom tibble tibble
#' @examples
#' if (require(rworldxtra)) {
#'
#' data(countriesHigh)
#' sv <- "Australia"
#' a <- subset(countriesHigh, SOVEREIGNT == sv)
#' b <- tri_mesh(a)
#' }
tri_mesh.SpatialPolygons <- function(x, ...) {
  pr4 <- proj4string(x)
  x0 <- x
  tabs <- spbabel::map_table(x)
  ## remove the branches (FIX ME - later we can keep them as long as triang doesn't change)
  tabs$b <- tabs$bXv <- NULL
  tabs$v <- tabs$t <- tabs$tXv <- NULL
  
  ## loop over $o (and for now, map_table each one separately as it easier to develop)
  for (i_obj in seq(nrow(x))) {
    tab0 <- spbabel::map_table(x[i_obj, ])
    #spbabel:::semi_cascade
    tab0$v$countingIndex <- seq(nrow(tab0$v))
    nonuq <- dplyr::inner_join(tab0$bXv, tab0$v, "vertex_")
    ps <- RTriangle::pslg(P = as.matrix(tab0$v[, c("x_", "y_")]),
                          S = do.call(rbind, lapply(split(nonuq, nonuq$branch_),
                                                    function(x) path2seg(x$countingIndex))))
    ## FIXME: need to pick sensible behaviour for a
    ## "tr" is the raw triangulation
    ## "tri" is the tibble of triangles (for now)
    tr <- RTriangle::triangulate(ps)
    ## process the holes if needed
    ## may be quicker than testing entire object
    if (any(!tab0$b$island_)) {
      holes <- spbabel::sp(dplyr::inner_join(dplyr::inner_join(dplyr::filter_(tab0$b, quote(!island_)), tab0$bXv, "branch_"), 
                                             tab0$v, "vertex_"))
      centroids <- matrix(unlist(lapply(split(tr$P[t(tr$T), ], rep(seq(nrow(tr$T)), each = 3)), .colMeans, 3, 2)), 
                          ncol = 2, byrow = TRUE)
      
      badtris <- !is.na(over(SpatialPoints(centroids), sp::geometry(holes)))
      if (any(badtris)) tr$T <- tr$T[!badtris, ]
    }
    
    v <- tibble::tibble(x_ = tr$P[,1], y_ = tr$P[,2], vertex_ = spbabel:::id_n(nrow(tr$P)))#  seq(nrow(tr$P)))
    ## tri (don't call something "t")                            ## not that this is the global tabs (FIXME)
    tri <- tibble::tibble(triangle_ = seq(nrow(tr$T)), object_ = tabs$o$object_[i_obj])
    tXv <- tibble::tibble(triangle_ = rep(tri$triangle_, each = 3), vertex_ = v$vertex_[as.vector(t(tr$T))])
    
    tabs$v <- dplyr::bind_rows(tabs$v, v)
    tabs$tXv <- dplyr::bind_rows(tabs$tXv, tXv)
    tabs$t <- dplyr::bind_rows(tabs$tri, tri)
    
  } ## end of loop i_obj that appends to v, t, tXv
  
  ## early prototype needs to have only one object
  ## but we can just keep it all now
 # tabs$o <- tabs$o[1,]  ## FIX ME
  

  ## finally add longitude and latitude
  tabs$meta <- tibble(proj = pr4, x = "x_", y = "y_")
  class(tabs) <- "trimesh"
  tabs
}

th3d <- function() {
  structure(list(vb = NULL, it = NULL, primitivetype = "triangle",
                 material = list(), normals = NULL, texcoords = NULL), .Names = c("vb",
                                                                                  "it", "primitivetype", "material", "normals", "texcoords"), class = c("mesh3d",
                                                                                                                                                        "shape3d"))
}

#' plot the triangles in the tables
#'
#' plot
#' @param x object from tri_mesh
#' @param ... args for underlying plotting
#'
#' @return the rgl mesh3d object, invisibly
#' @export
#' @importFrom rgl shade3d
#' @examples
#' example(tri_mesh)
#' if(exists("b")) { 
#'  plot(b)
#'  }
plot.trimesh <- function(x, ...) {
  tt <- th3d()
  x0 <- spbabel:::semi_cascade(x, tables = c("o", "t", "tXv", "v"))
  v <- x0$v %>% mutate(v_count = row_number())
  tXv <- x0$tXv %>% inner_join(v, "vertex_")
  tt$vb <- t(cbind(x0$v$x_, x0$v$y_, 0, 1))
  tt$it <- t(matrix(tXv$v_count, ncol = 3, byrow = TRUE))
  if (!requireNamespace("rgl", quietly = TRUE))
    stop("rgl required")
  rgl::shade3d(tt, ...)
  invisible(tt)
}

#' Title
#'
#' @param x object from tri_mesh
#' @param halo show the radius
#' @param ... passed to plot
#' @param rad radius
#'
#' @return the mesh object, invisibly
#' @export
#' @examples
#' example(tri_mesh)
#' if(exists("b")) { 
#'  ##globe(b, halo = TRUE)
#'  }
globe <- function(x, ...) {
  UseMethod("globe")
}

#' @rdname globe
#' @export
globe.trimesh <- function(x, halo = FALSE, ..., rad = 1) {
  gproj <- sprintf("+proj=geocent +a=%f +b=%f", rad, rad)
  p4 <- x$meta$proj[1]
  ll <- cbind(as.matrix(x$v[, c("x_", "y_")]), 0)
  if (grepl("longlat", p4)) ll <- ll * pi / 180
  xyz <- proj4::ptransform(ll, src.proj = p4, dst.proj = gproj)
  tt <- th3d()
  tt$vb <- t(cbind(xyz, 1))
  tt$it <- t(matrix(x$tXv$vertex_, ncol = 3, byrow = TRUE))
  if (!requireNamespace("rgl", quietly = TRUE))
    stop("rgl required")
  rgl::shade3d(tt, ...)
  if (halo) rgl::spheres3d(0, 0, 0, radius = rad * 0.99, fog = FALSE, specular = "black", col = "dodgerblue", alpha = 0.4)
  invisible(tt)
}