# loading dataset for test ------------------------------------------------

shp <- rgdal::readOGR(here::here("data", "shapes_tiranideos.shp"))
shape.america <- rgdal::readOGR(here::here("data", "shape_america2.shp"))
grid.resol <- 2

# modified version of grid filter function (Hidasi-Neto) ------------------

GridFilter <- function(shape,
                       resol = 1,
                       prop = 0)
  {
  grid <- raster::raster(extent(shape))
  terra::res(grid) <- resol
  sp::proj4string(grid) <- sp::proj4string(shape)
  gridpolygon <- rasterToPolygons(grid)
  drylandproj <- sp::spTransform(shape, CRS("+proj=laea"))
  gridpolproj <- sp::spTransform(gridpolygon, CRS("+proj=laea"))
  gridpolproj$layer <- c(1:length(gridpolproj$layer))
  areagrid <- rgeos::gArea(gridpolproj, byid=T)
  dry.grid <- raster::intersect(drylandproj, gridpolproj)
  areadrygrid <- rgeos::gArea(dry.grid, byid=T)
  info <- cbind(dry.grid$layer, areagrid[dry.grid$layer], areadrygrid)
  dry.grid$layer <- info[,3]/info[,2]
  dry.grid <- sp::spTransform(dry.grid, CRS(proj4string(shape)))
  dry.grid.filtered <- dry.grid[dry.grid$layer >= prop,]
}


# build grid and composition function -------------------------------------

grid_comp <- function(shp, 
                      grid.resol, 
                      prop, 
                      plot.rich = TRUE){
  filter_grid <- GridFilter(shape = shape.america, resol = grid.resol, prop = 1)
  slot(filter_grid, "data") <- cbind("ID" = 1:length(filter_grid),
                                     slot(filter_grid, "data"))
  pa.shp <- letsR::lets.presab.grid(shapes = shp, 
                                    grid = filter_grid, 
                                    sample.unit = "ID")
  ext <- extent(pa.shp$grid)
  xmn <- ext[1]
  xmx <- ext[2]
  ymn <- ext[3]
  ymx <- ext[4]
  r <- raster::raster(resolution = grid.resol, xmn = xmn, xmx = xmx,
                      ymn = ymn, ymx = ymx)
  rich <- terra::rasterize(pa.shp$grid, r, field = "ID")
  values(rich)[which(!is.na(values(rich)) == TRUE)] <- rowSums(pa.shp$PAM)
  xy <- terra::xyFromCell(rich, 1:ncell(rich))
  colnames(xy) <- c("Longitude(x)", "Latitude(y)")
  pam <- cbind(ID = which(!is.na(values(rich)) == TRUE), xy[which(!is.na(values(rich)) == TRUE), ], pa.shp$PAM)
  if(plot.rich == TRUE){
    plot(rich)
    map(add = TRUE)
    return(pam)
  } else{
    pam
  }
}


# testing with real data --------------------------------------------------

grid_test <- grid_comp(shp = shp, 
                       grid.resol = 2,
                       prop = 0, 
                       plot.rich = TRUE)




# obtaining long data format for coordinates ------------------------------


coords_grid <- grid_test[, c(2, 3)]
rich_spp <- rowSums(grid_test[, c(4:ncol(grid_test))])
grid_comp_pa <- grid_test[, c(4:ncol(grid_test))][which(rich_spp > 0), ]
rich_spp_pa <- rowSums(grid_comp_pa)
coords_grid_pa <- coords_grid[which(rich_spp > 0), ]

coords_rich_pa <- cbind(coords_grid_pa, rich_spp_pa)
coords_long <- do.call(rbind, apply(coords_rich_pa, MARGIN = 1, function(x){
  matrix(rep(x[1:2], times = x[3]),
         nrow = x[3],
         ncol = 2, 
         dimnames = list(1:x[3], c("x", "y")),
         byrow = T)
}, simplify = FALSE))
names_long <- unlist(apply(grid_comp_pa, MARGIN = 1, function(x) names(which(x == 1))))
units_long <- unlist(lapply(rownames(grid_comp_pa), MARGIN = 1, function(x) names(x)))
data_long <- data.frame(coords_long, names_long, ID = rep(rownames(grid_comp_pa), rich_spp_pa))

