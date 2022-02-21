

# libraries ---------------------------------------------------------------

library(shiny)
library(leaflet)
library(maps)
library(raster)
library(sp)
library(sf)
library(terra)
library(ggplot2)

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


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Biogeographylo - Application to analyze spatial patterns in Historical Biogeography"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "shp_input", 
                label = "Input shapefile of occurrences", 
                multiple = TRUE,
                accept = c('.shp','.dbf','.sbn','.sbx','.shx','.prj')
      ),
      fileInput(inputId = "phylo_input",
                label = "Input phylogenetic tree",
                accept = c('.txt', '.new')
      ),
      actionButton("shp_tirani", "Example shapefile"),
      textInput(inputId = "grid_size", label = "Grid size cells in decimal degrees:"),
      textInput(inputId = "prop_cell", label = "Proportion of the cells used in the grid"),
      checkboxGroupInput(inputId = "metrics_calc", 
                         label = "Which metrics should be calculated?", 
                         choices = c("PD" = "PD_faith", 
                                     "PE" = "PE_rosauer", 
                                     "PD model-based" = "PD_model",
                                     "PE model-based" = "PE_model"
                         ),
                         selected = c("PD_faith",
                                      "PE_rosauer", 
                                      "PD_model", 
                                      "PE_model"
                         ),
                         width = "100%"
      )
      
    ),
    
    # plot phylogeny and map
    mainPanel(
      fluidRow(
        plotOutput(outputId = "mapPlot", height = 500)
      ),
      fluidRow(
        plotOutput(outputId = "phyloPlot")
      )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # input shapefile of occurrence records from user
  map <- reactive({
    req(input$shp_input)
    
    # shpdf is a data.frame with the name, size, type and
    # datapath of the uploaded files
    shpdf <- input$shp_input
    
    # The files are uploaded with names
    # 0.dbf, 1.prj, 2.shp, 3.xml, 4.shx
    # (path/names are in column datapath)
    # We need to rename the files with the actual names:
    # fe_2007_39_county.dbf, etc.
    # (these are in column name)
    
    # Name of the temporary directory where files are uploaded
    tempdirname <- dirname(shpdf$datapath[1])
    
    # Rename files
    for (i in 1:nrow(shpdf)) {
      file.rename(
        shpdf$datapath[i],
        paste0(tempdirname, "/", shpdf$name[i])
      )
    }
    
    # Now we read the shapefile with readOGR() of rgdal package
    # passing the name of the file with .shp extension.
    
    # We use the function grep() to search the pattern "*.shp$"
    # within each element of the character vector shpdf$name.
    # grep(pattern="*.shp$", shpdf$name)
    # ($ at the end denote files that finish with .shp,
    # not only that contain .shp)
    map <- readOGR(paste(tempdirname,
                         shpdf$name[grep(pattern = "*.shp$", shpdf$name)],
                         sep = "/"
    ))
    map <- sf::st_as_sf(x = map)
    map
  })
  
  # input example shapefile of occurrence records from Tiranidae
  map <- reactive({
    req(input$shp_tirani)
    map <- rgdal::readOGR(here::here("data", "shape_america2.shp")
    )
    map <- sf::st_as_sf(x = map)
    map
  })
  
  # user option to input phylogenetic tree
  phylo_out <- reactive({
    req(input$phylo_input)
    phylo_out <- ape::read.tree(file = input$phylo_input$datapath)
    phylo_out
  })
  
  # example with tiranidae phylogenet tree
  phylo_out <- reactive({
    req(input$shp_tirani)
    phylo_out <- ape::read.tree(file = here::here("data", "Tree_TF400Howard_Pruned.tre"))
    phylo_out
  })
  
  # map output with shapefile
  output$mapPlot <- renderPlot({
    
    ggplot() +
      geom_sf(data = map(), aes(geometry = geometry))
    
  })
  
  # output$mapPlot <- renderPrint({
  #   map()
  # })
  
  # phylo output with phylogenetic tree
  output$phyloPlot <- renderPlot({
    ape::plot.phylo(phylo_out(), type = "fan", show.tip.label = FALSE, show.node.label = T)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
