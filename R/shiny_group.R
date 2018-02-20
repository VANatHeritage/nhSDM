# shiny_group

#' Interactively view (in a leaflet map) and select grouping distance
#'
#' @param spf input spatial features (sp or sf spatial object)
#'
#' @import leaflet
#' @import shiny
#' @importFrom sp spTransform
#' @importFrom grDevices rainbow
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe_expl")
#' shiny_group(spf)
#' }


shiny_group <- function(spf) {
  shinyApp(ui = bootstrapPage(
    tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
    leafletOutput("map", width = "100%", height = "100%"),
    absolutePanel(top = 10, right = 10,
                  sliderInput("sep", "Separation distance", min = 0, max = 50000, value = 1000, step = 100)
  )), server = function(input, output, session) {

    pg <- reactive({
      pg <- nh_group(spf, sep.dist = input$sep, union = FALSE)
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
  })
}
