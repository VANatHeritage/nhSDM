library(raster)
library(NHsdm)
library(rgdal)
library(leaflet)
library(shiny)

ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("map", width = "100%", height = "100%"),
  absolutePanel(top = 10, right = 10,
                sliderInput("sep", "Seperation distance", min = 0, max = 10000, value = 1000, step = 100)
))

server <- function(input, output, session) {
  # a<-raster("D:/SDM/Tobacco/env_vars/Tobacco/AnnMnTemp.tif")
  spP <- readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe_expl")
  
  pg <- reactive({
    pg <- poly_group(spP, sep.dist = input$sep, union = FALSE)
    pg <- spTransform(pg, CRSobj = "+init=epsg:4326")
    pg$group <- as.factor(pg$group)
    pg
  })
  
  output$map <- renderLeaflet({
    pg <- pg()
    factpal <- colorFactor(rainbow(length(levels(pg$group))), pg$group)
    
    leaflet(pg) %>%
      addTiles() %>%
      addPolygons(color = ~factpal(group), stroke = TRUE, weight = 20, fillColor = ~factpal(group), fillOpacity = 0.5, popup = ~group)
  })
}

shinyApp(ui, server)
