
library(shiny)
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

shp <- rgdal::readOGR(here::here("data", "shapes_tiranideos.shp"))
shape.america <- rgdal::readOGR(here::here("data", "shape_america2.shp"))
grid.resol <- 4

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
    titlePanel("Biogeographylo - an application to analyse spatial patterns in evolutionary history"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"), # here goes the phylogeny
           plotOutput("gridmap") # here goes the grid map 
                                    # both of them must be reactive objects
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 


shinyApp(ui = ui, server = server)
