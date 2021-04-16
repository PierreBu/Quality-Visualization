library(shiny)
library(dplyr)
library(ggplot2)
library(patchwork)


# Global variables --------------------------------------------------------

RED_LIMIT <- 25
Y_SPACING <- 0.2
ANNOTATION_SIZE <- 2

# Functions in global environment -----------------------------------------

adaptReferenceToEstimatedMatrix <- function(referenceToBeAdaptedToEstimated, referenceSampleNames, referenceElements, estimatedSampleNames, estimatedElements) {
  
  adaptedMatrix <- referenceToBeAdaptedToEstimated[match(estimatedSampleNames, referenceSampleNames),
                                                   match(estimatedElements, referenceElements)]
  
  return(adaptedMatrix)
}

removeSamplesThatAreNotStandards <- function(estimatedMatrixWithNonStandards, referenceSamplesNames, estimatedSamplesNames) {
  
  if (length(setdiff(estimatedSamplesNames, referenceSamplesNames)) > 0) {
    estimatedMatrixWithNonStandardsRemoved <- estimatedMatrixWithNonStandards[-which(estimatedSamplesNames %in% setdiff(estimatedSamplesNames, referenceSamplesNames)), ]
  }
  else {
    estimatedMatrixWithNonStandardsRemoved <- estimatedMatrixWithNonStandards
  }
  
  return(estimatedMatrixWithNonStandardsRemoved)
}

CalculateRelativeDeviationValue <- function(sampleMatrix, referenceMatrix) {
  
  relativeDeviationMatrix <- (sampleMatrix - referenceMatrix) / referenceMatrix * 100
  
  return(relativeDeviationMatrix)
}

CalculateRelativeDeviationStandardDeviation <- function(sampleValueMatrix, referenceValueMatrix,
                                                        sampleStandardDeviationMatrix, referenceStandardDeviationMatrix) {
  
  relativeDeviationStandardDeviationMatrix <- sqrt((sampleStandardDeviationMatrix / sampleValueMatrix)^2 + (referenceStandardDeviationMatrix / referenceValueMatrix)^2) * sampleValueMatrix / referenceValueMatrix * 100
  
  return(relativeDeviationStandardDeviationMatrix)
}


GetAsMatrix <- function(dataFilePath){
  concentrationDf <- read.table(dataFilePath, header = TRUE, sep = ";")
  sampleNames <- as.character(concentrationDf[, 1])
  concentrationMatrix <- as.matrix(concentrationDf[,-1])
  concentrationMatrix <- apply(concentrationMatrix, 2, as.numeric)
  row.names(concentrationMatrix) <- sampleNames
  
  return(concentrationMatrix)
}

createYSpacing <- function(elementSampleNames) {

  numberOfSamples <- length(elementSampleNames)
  samplesUniqueNames <- unique(elementSampleNames)

  ySpacing <- rep(NA, numberOfSamples)
  lowerYLimit <- 0

  for (name in samplesUniqueNames) {
    numberOfOccurenceOfCurrentName <- sum(elementSampleNames == name)
    ySpacing[elementSampleNames == name] <- seq(from = lowerYLimit,
                                                to = Y_SPACING * (numberOfOccurenceOfCurrentName - 1) + lowerYLimit,
                                                by = Y_SPACING)
    lowerYLimit <- max(ySpacing, na.rm = TRUE) + Y_SPACING
  }

  return(- ySpacing)
}


# UI definition -----------------------------------------------------------

ui <- fluidPage(
  fileInput("file", "File upload"),
  actionButton("setAsEstimatedConcentration", "Set as estimated concentration file"),
  actionButton("setAsEstimatedStandardDeviation", "Set as estimated standard deviation file"),
  actionButton("setAsReferenceConcentration", "Set as reference concentration file"),
  actionButton("setAsReferenceStandardDeviation", "Set as reference standard deviation file"),
  numericInput(inputId = "ncol", label = "Column number", value = 8),
  textInput("graphName", "Graph name", ""),
  actionButton("createGraphs", "Create graphs")
)

# Server definition -------------------------------------------------------

server <- function(input, output, session) {

# Base matrixes definition ------------------------------------------------

  reactiveBaseMatrixes <- reactiveValues(referenceConcentration = matrix(),
                                         referenceStandardDeviation = matrix(), 
                                         estimatedConcentration = matrix(),
                                         estimatedStandardDeviation = matrix())

# Adapted matrixes definition and implementation --------------------------
  
  reactiveAdaptedMatrixes <- reactiveValues(referenceConcentration = matrix(),
                                            referenceStandardDeviation = matrix(),
                                            estimatedConcentration = matrix(),
                                            estimatedStandardDeviation = matrix())
  
  reactiveAdaptedMatrixes$estimatedConcentration <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    removeSamplesThatAreNotStandards(reactiveBaseMatrixes$estimatedConcentration,
                                     reactiveBaseReferenceParameters$sampleNames(),
                                     reactiveBaseEstimatedParameters$sampleNames())
    
  })
  
  reactiveAdaptedMatrixes$referenceConcentration <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    adaptReferenceToEstimatedMatrix(reactiveBaseMatrixes$referenceConcentration,
                                    reactiveBaseReferenceParameters$sampleNames(), reactiveBaseReferenceParameters$sampleElements(),
                                    reactiveAdaptedEstimatedParameters$sampleNames(), reactiveAdaptedEstimatedParameters$sampleElements())
      
  })
  
  reactiveAdaptedMatrixes$referenceStandardDeviation <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    if (all(is.na(reactiveBaseMatrixes$referenceStandardDeviation))) 
    {
      matrix(rep(0, reactiveAdaptedEstimatedParameters$sampleNumber() * reactiveAdaptedEstimatedParameters$elementNumber()),
             nrow = reactiveAdaptedEstimatedParameters$sampleNumber(), ncol = reactiveAdaptedEstimatedParameters$elementNumber())
    }
    else {
      adaptReferenceToEstimatedMatrix(reactiveBaseMatrixes$referenceStandardDeviation,
                                      reactiveBaseReferenceParameters$sampleNames(), reactiveBaseReferenceParameters$sampleElements(),
                                      reactiveAdaptedEstimatedParameters$sampleNames(), reactiveAdaptedEstimatedParameters$sampleElements())
    }
  })
  
  reactiveAdaptedMatrixes$estimatedStandardDeviation <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    if (all(is.na(reactiveBaseMatrixes$estimatedStandardDeviation))) 
    {
      matrix(rep(0, reactiveAdaptedEstimatedParameters$sampleNumber() * reactiveAdaptedEstimatedParameters$elementNumber()),
             nrow = reactiveAdaptedEstimatedParameters$sampleNumber(), ncol = reactiveAdaptedEstimatedParameters$elementNumber())
    }
    else {
      removeSamplesThatAreNotStandards(reactiveBaseMatrixes$estimatedStandardDeviation,
                                       reactiveBaseReferenceParameters$sampleNames(),
                                       reactiveBaseEstimatedParameters$sampleNames())
    }
  })

# Calculated matrixes definition and implementation -----------------------
  
  reactiveCalculatedMatrixes <- reactiveValues(relativeDeviationValue = matrix(),
                                               relativeDeviationStandardDeviation = matrix())
  
  reactiveCalculatedMatrixes$relativeDeviationValue <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    CalculateRelativeDeviationValue(reactiveAdaptedMatrixes$estimatedConcentration(),
                                    reactiveAdaptedMatrixes$referenceConcentration())
    
  })
  
  reactiveCalculatedMatrixes$relativeDeviationStandardDeviation <- reactive({
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    CalculateRelativeDeviationStandardDeviation(reactiveAdaptedMatrixes$estimatedConcentration(),
                                                reactiveAdaptedMatrixes$referenceConcentration(),
                                                reactiveAdaptedMatrixes$estimatedStandardDeviation(),
                                                reactiveAdaptedMatrixes$referenceStandardDeviation())
  })

# Estimated parameters definition and implementation ----------------------
  
  reactiveBaseEstimatedParameters <- reactiveValues(sampleNames = NULL,
                                                    sampleElements = NULL,
                                                    sampleNumber = NULL,
                                                    elementNumber = NULL)
  
  reactiveBaseEstimatedParameters$sampleNames <- reactive({row.names(reactiveBaseMatrixes$estimatedConcentration)})
  reactiveBaseEstimatedParameters$sampleElements <- reactive({colnames(reactiveBaseMatrixes$estimatedConcentration)})
  reactiveBaseEstimatedParameters$sampleNumber <- reactive({length(reactiveBaseEstimatedParameters$sampleNames())})
  reactiveBaseEstimatedParameters$elementNumber <- reactive({length(reactiveBaseEstimatedParameters$sampleElements())})
  
# Reference parameters definition and implementation ----------------------

  reactiveBaseReferenceParameters <- reactiveValues(sampleNames = NULL,
                                                    sampleElements = NULL,
                                                    sampleNumber = NULL,
                                                    elementNumber = NULL)
  
  reactiveBaseReferenceParameters$sampleNames <- reactive({row.names(reactiveBaseMatrixes$referenceConcentration)})
  reactiveBaseReferenceParameters$sampleElements <- reactive({colnames(reactiveBaseMatrixes$referenceConcentration)})
  reactiveBaseReferenceParameters$sampleNumber <- reactive({length(reactiveBaseReferenceParameters$sampleNames())})
  reactiveBaseReferenceParameters$elementNumber <- reactive({length(reactiveBaseReferenceParameters$sampleElements())})
  
# Adapted Estimated parameters definition and implementation --------------
  
  reactiveAdaptedEstimatedParameters <- reactiveValues(sampleNames = NULL,
                                                       sampleElements = NULL,
                                                       sampleNumber = NULL,
                                                       elementNumber = NULL)
  
  reactiveAdaptedEstimatedParameters$sampleNames <- reactive({row.names(reactiveAdaptedMatrixes$estimatedConcentration())})
  reactiveAdaptedEstimatedParameters$sampleElements <- reactive({colnames(reactiveAdaptedMatrixes$estimatedConcentration())})
  reactiveAdaptedEstimatedParameters$sampleNumber <- reactive({length(reactiveBaseEstimatedParameters$sampleNames())})
  reactiveAdaptedEstimatedParameters$elementNumber <- reactive({length(reactiveBaseEstimatedParameters$sampleElements())})
  
# Adapted Reference parameters definition and implementation --------------
  
  reactiveAdaptedReferenceParameters <- reactiveValues(sampleNames = NULL,
                                                       sampleElements = NULL,
                                                       sampleNumber = NULL,
                                                       elementNumber = NULL)
  
  reactiveAdaptedReferenceParameters$sampleNames <- reactive({row.names(reactiveAdaptedMatrixes$referenceConcentration())})
  reactiveAdaptedReferenceParameters$sampleElements <- reactive({colnames(reactiveAdaptedMatrixes$referenceConcentration())})
  reactiveAdaptedReferenceParameters$sampleNumber <- reactive({length(reactiveBaseReferenceParameters$sampleNames())})
  reactiveAdaptedReferenceParameters$elementNumber <- reactive({length(reactiveBaseReferenceParameters$sampleElements())})
  
# Reactive required conditions --------------------------------------------

  areEstimatedAndReferenceMatrixesDefined <- reactive({
    
    (all(is.na(reactiveBaseMatrixes$estimatedConcentration)) == FALSE && all(is.na(reactiveBaseMatrixes$referenceConcentration)) == FALSE)
    
  })

# Input observers ---------------------------------------------------------
  
  observeEvent(input$setAsReferenceConcentration, {
    req(input$file)
    
    reactiveBaseMatrixes$referenceConcentration <- GetAsMatrix(input$file$datapath)
  })
  
  observeEvent(input$setAsReferenceStandardDeviation, {
    req(input$file)
    
    reactiveBaseMatrixes$referenceStandardDeviation <- GetAsMatrix(input$file$datapath)
  })
  
  observeEvent(input$setAsEstimatedConcentration, {
    req(input$file)

    updateTextInput(session = session, inputId = "graphName", label = "Graph name", value = substr(input$file$name, 1, nchar(input$file$name)-4))
    
    reactiveBaseMatrixes$estimatedConcentration <- GetAsMatrix(input$file$datapath)
  })
  
  observeEvent(input$setAsEstimatedStandardDeviation, {
    req(input$file)
    
    reactiveBaseMatrixes$estimatedStandardDeviation <- GetAsMatrix(input$file$datapath)
  })
  
  observeEvent(input$createGraphs, {
    
    req(areEstimatedAndReferenceMatrixesDefined())
    
    myPlotList <- list()
    
    ySpacing <- createYSpacing(row.names(reactiveCalculatedMatrixes$relativeDeviationValue()))
    
    for (elementName in reactiveAdaptedEstimatedParameters$sampleElements()) {
      
      myPlot <- ggplot()
      
      minX <- min(as.numeric(reactiveCalculatedMatrixes$relativeDeviationValue()[ , elementName] - reactiveCalculatedMatrixes$relativeDeviationStandardDeviation()[ , elementName]), na.rm = TRUE)
      maxX <- max(as.numeric(reactiveCalculatedMatrixes$relativeDeviationValue()[ , elementName] + reactiveCalculatedMatrixes$relativeDeviationStandardDeviation()[ , elementName]), na.rm = TRUE)
      xRange <- maxX - minX
      firstOccurenceOfEachStandard <- match(unique(reactiveAdaptedEstimatedParameters$sampleNames()), reactiveAdaptedEstimatedParameters$sampleNames())
      
      plotDf <- data.frame(element = reactiveCalculatedMatrixes$relativeDeviationValue()[ , elementName],
                           standardDeviation = reactiveCalculatedMatrixes$relativeDeviationStandardDeviation()[ , elementName],
                           height = ySpacing,
                           sampleName = reactiveAdaptedEstimatedParameters$sampleNames(),
                           referenceConcentration = reactiveAdaptedMatrixes$referenceConcentration()[ , elementName])
      
      myPlot <- myPlot +
        geom_errorbarh(plotDf, mapping = aes(xmin = element - standardDeviation, xmax = element + standardDeviation, y = height), height = 0.2, size = 0.1) +
        geom_point(plotDf, mapping = aes(x = element, y = height, color = abs(element), shape = sampleName)) +
        annotate(geom="text", x = minX - 0.2 * xRange,
                 y = ySpacing[firstOccurenceOfEachStandard],
                 label = reactiveAdaptedMatrixes$referenceConcentration()[firstOccurenceOfEachStandard, elementName],
                 color="black", hjust = 0, size = ANNOTATION_SIZE) +
        xlab(paste(elementName, " deviation (%)", sep = "")) +
        xlim(minX - 0.2 * xRange, maxX) +
        labs(color='Deviation (%)') +
        labs(shape='Standard Name') +
        scale_color_gradient(low = "blue", high = "red", limits=c(0, RED_LIMIT), na.value = "red") + 
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      myPlotList[[elementName]] <- myPlot
    }
    myWrapPlots <- wrap_plots(myPlotList, guides = "collect", ncol = input$ncol) & theme(legend.position = 'bottom')
    myWrapPlots + plot_annotation(
      title = input$graphName,
      subtitle = 'Quality Vizualisation plots show the deviation of estimated from certified standard concentrations',
      caption = 'Contact: burckel@ipgp.fr'
    )
    ggsave(filename = paste(input$graphName, ".pdf", sep = ""),
            width = 2 * input$ncol, height = 1.5 * ceiling(length(reactiveAdaptedEstimatedParameters$sampleElements()) / input$ncol),
            dpi = 150,  units = "in", device = "pdf")
    print("Finished")
  })
}

shinyApp(ui, server)