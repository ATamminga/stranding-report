DetrendDEM <- function(src.file, map.bound, water.poly, cell.size, dem.out) {
    
  #read bounding box 
  bound.poly <- readOGR(map.bound)
  coord.ref <- bound.poly@proj4string
  bound.box <- extent(bound.poly)
    
  #read xyz file of bed elevations
  bed.points <- read.csv(src.file, header=TRUE)
  bed.points <- na.omit(bed.points)
  coordinates(bed.points) <- ~x+y
  proj4string(bed.points) <- coord.ref
  bed.points <- crop(bed.points, bound.box)
  
  #read high water surface shapefile
  wse.poly <- readOGR(water.poly)
  wse.poly <- crop(wse.poly, bound.box)
  
  #create variables describing data range
  x.out <- seq(bound.box@xmin, bound.box@xmax, cell.size)
  y.out <- seq(bound.box@ymin, bound.box@ymax, cell.size)
  
  #interpolate dem
  dem <- interp(bed.points$x, bed.points$y, bed.points$z, x.out, y.out)
  dem <- raster(dem, crs = coord.ref)
  
  #detrend dem 
  wse.points <- spsample(wse.poly, 10000, "regular")
  wse.points <- data.frame(coordinates(wse.points), extract(dem, wse.points))
  names(wse.points) <- c("x","y","z")
  wse.points <- na.omit(wse.points)
  coordinates(wse.points) <- ~x+z
  plane.fit <- lm(wse.points$z ~ wse.points$x + wse.points$y)
  plane.coefs <- tidy(plane.fit)
  wse.plane <- expand.grid(x = x.out, y = y.out) %>% 
    mutate(z = plane.coefs$estimate[2] * x + plane.coefs$estimate[3]*y + plane.coefs$estimate[1])
  coordinates(wse.plane) <- ~x + y
  gridded(wse.plane) <- TRUE
  wse.plane <- raster(wse.plane)
  dem.detrend <- dem - wse.plane
  dem.detrend <- dem.detrend - minValue(dem.detrend) + .001
  dem.mask <- rasterize(wse.poly, dem.detrend)
  dem.detrend <- mask(dem.detrend, dem.mask)
  
  writeRaster(dem.detrend, dem.out, overwrite = TRUE)
  
  return(dem.detrend)
  
}

RankElevs <- function(dem.in, map.bound, n.sample){
  
  #read bounding box 
  bound.poly <- readOGR(map.bound)
  coord.ref <- bound.poly@proj4string
  bound.box <- extent(bound.poly)
  
  #read dem
  dem <- raster(dem.in)
  names(dem) <- c('z')
  
  #create sample polygons
  sample.length <- (bound.box@xmax-bound.box@xmin)/(n.sample)
  breaks <- as.list(seq(from = bound.box@xmin, to = bound.box@xmax-sample.length, by = sample.length ))
  
  defineExtent <- function(Y) {
    sample.polygon <- extent(breaks[[Y]], breaks[[Y]]+sample.length, bound.box@ymin, bound.box@ymax)
  }
  
  polys<- lapply(1:length(breaks), defineExtent)
  
  SubsetDEM <- function(X) {
    dem.subset <- crop(dem, polys[[X]]) %>%
      as.data.frame(xy=TRUE) %>% 
      filter(!is.na(z)) %>% 
      arrange(z) %>% 
      mutate(z = z - min(z), 
             ecdf = dense_rank(z)/length(z),
             wetted.area = ecdf*length(ecdf)*(prod(res(dem))),
             wetted.width = wetted.area/sample.length,
             id = as.factor(X)) %>% 
      filter(ecdf < .99)
    model <- dem.subset %>%
      nls(formula = "wetted.width ~ k * z^a", start=list(k=25,a=1))
    model.coef <- coef(model)
    dem.subset <- dem.subset %>%
      mutate(model.fit = predict(model))
    outputs <- list("dem.subset" = dem.subset, "model" = model, "box" = polys[[X]], "model.coef" = model.coef)
    # return(outputs)
    
  }
  
    dems.out <-mclapply(1:length(polys), SubsetDEM, mc.cores=7)

    

}

PlotElevs <- function(elevs.in, water.poly, map.bound){
  
  #read bounding box 
  bound.poly <- readOGR(map.bound)
  coord.ref <- bound.poly@proj4string
  bound.box <- extent(bound.poly)
  
  wse.poly <- readOGR(water.poly)
  
  poly.colors <- brewer.pal(9, "Set1")
  
  
  dems <- do.call(rbind, simplify2array(elevs.in)[1,])
  models <- simplify2array(elevs.in)[2,]
  extents <- simplify2array(elevs.in)[3,]
  model.coef <- simplify2array(elevs.in)[4,]
  
  model.coef<- data.frame(do.call(rbind, unname(model.coef))) %>% 
    mutate("id" = as.factor(seq(1, length(extents), 1)))
  
  
  polys <- do.call(bind, lapply(1:length(extents), FUN = function(x) as(extents[[x]], "SpatialPolygons")))
  
  theme_map <- function(...) {
    theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.background = element_rect(fill = "#f5f5f2", color = NA), 
        # panel.background = element_rect(fill = "#f5f5f2", color = NA), 
        # legend.background = element_rect(fill = "#f5f5f2", color = NA)
        ...
      )
  }
  
  p <- dems %>% 
    ggplot()+
    facet_wrap(~id)+
    geom_line(aes(x=wetted.width, y = z, color = id), size = .5)+
    geom_line(aes(x=model.fit, y =z), lty="dotted")+
    geom_text(data = model.coef, aes(x=400, y = 8.5, label = paste("b ==",round(k, digits = 2),"* D ^",round(a, digits=2))), size = 3, parse = TRUE)+
    ylim(0,9)+
    theme_bw()+
    xlab("b (width, m)")+
    ylab("D (elev. above datum, m)")+
    scale_color_manual("ID", values = poly.colors)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    
  
  g <- dems %>% 
    ggplot()+
    geom_line(aes(x=model.fit, y = z, color = id))+
    theme_bw()+
    scale_color_manual("ID", values = poly.colors)+
    xlim(0, 1000)+
    ylim(0, 10)+
    xlab("b (width, m)")+
    ylab("D (elev. above datum, m)")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  q <- dems %>% 
    ggplot()+
    geom_tile(aes(x=x, y=y, fill = z))+
    geom_path(data = tidy(polys), aes(x=long, y =lat, group=group, color = id))+
    geom_path(data = tidy(wse.poly), aes(x=long,y=lat, group=group), size = .15)+
   # theme_map()+
    scale_fill_viridis(limits = c(0, 10), name = "elev (m)")+
    coord_equal()+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_color_manual("ID", values = poly.colors, guide = FALSE)+
    xlim(bound.box@xmin, bound.box@xmax)+
    ylim(bound.box@ymin, bound.box@ymax)
    


  return(plots <- list(p,g, q))
  
}
