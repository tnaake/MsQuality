#' @name plotMetric
#'
#' @title Visualize a quality metric
#'
#' @description
#' The function \code{plotMetric} visualizes the metric values per sample. The 
#' function accepts the output of \code{calculateMetrics} or,
#' \code{calculateMetricsFromSpectra}, or 
#' \code{calculateMetricsFromMsExperiment} and a vector specifying the metric 
#' to display.
#' 
#' @details
#' \code{plotMetric} will select all columns that start with
#' \code{metric}. The different levels in the \code{name} column in the 
#' returned tibble correspond to the columns that were selected and do not 
#' contain the \code{metric} prefix. In case there is no additional specification 
#' (e.g. for the metric \code{rtDuration} only the column \code{rtDuration} will 
#' be selected), the \code{name} column will include the \code{metric} 
#' (\code{rtDuration}).
#' 
#' @param qc \code{matrix}/\code{data.frame}
#' @param metric \code{character}
#' @param plotly \code{logical(1)}
#' 
#' @return \code{gg} \code{plotly}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom ggplot2 ggplot geom_point aes_string scale_colour_brewer theme_bw
#' @importFrom ggplot2 xlab ggtitle guides guide_legend theme element_text
#' @importFrom ggplot2 element_blank
#' @importFrom plotly ggplotly partial_bundle
#' 
#' @export
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' msexp <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(msexp) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' 
#' library(Spectra)
#' ## import the data and add it to the msexp object
#' spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#'     
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange)
#' qc <- calculateMetricsFromMsExperiment(msexp = msexp, metrics = metrics, 
#'     msLevel = 1, relativeTo = "Q1", change = "jump")
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' 
#' ## do the actual plotting
#' plotMetric(qc, metric = "areaUnderTic", plotly = TRUE)
plotMetric <- function(qc, metric = "areaUnderTic", plotly = TRUE) {
    
    qc_tbl_l <- plotMetricTibble(qc = qc, metric = metric)
    
    g <- ggplot(qc_tbl_l) +
        geom_point(aes_string(x = "rowname", y = "value", col = "name")) +
        scale_colour_brewer(palette = "Set1") + theme_bw() +
        xlab("sample") + ggtitle(metric) +
        guides(shape = guide_legend(override.aes = list(size = 5))) +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        theme(
            axis.text.x = element_text(angle = 90, size = 10), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    if (plotly)
        g |> 
            ggplotly(tooltip = c("x", "y")) |>
            partial_bundle()
    else 
        g
}

#' @name plotMetricTibble
#'
#' @title Helper function for plotMetric
#'
#' @description
#' The function \code{plotMetricTibble} is a helper function for the function
#' \code{plotMetric}. It returns a tibble in long format that is interpretable
#' by \code{ggplot2}.
#' 
#' @details
#' \code{plotMetricRibble} will select all columns that start with
#' \code{metric}. The different levels in the \code{name} column in the returned 
#' tibble correspond to the columns that were selected and do not contain the
#' \code{metric} prefix. In case there is no additional specification 
#' (e.g. for the metric \code{rtDuration} only the column \code{rtDuration} will 
#' be selected), the \code{name} column will include the \code{metric} 
#' (\code{rtDuration}). 
#' 
#' @param qc \code{data.frame}
#' @param metric \code{character}
#' 
#' @return \code{tibble}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom stringr str_remove
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ProtGenerics spectra
#' 
#' @export
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' msexp <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(msexp) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(msexp) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
#'     sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#'
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange)
#' qc <- calculateMetricsFromMsExperiment(msexp = msexp, metrics = metrics, 
#'     msLevel = 1, relativeTo = "Q1", change = "jump")
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' plotMetricTibble(qc, metric = "areaUnderTic")
plotMetricTibble <- function(qc, metric) {
    
    cols <- grep(colnames(qc), pattern = metric)
    
    if (length(cols) == 0) stop("'metric' not in qc")
    
    qc_df <- as.data.frame(qc)
    qc_df <- qc_df[, cols, drop = FALSE]
    
    ## remove from the colnames the "metric" part and the first _ separator
    ## assign then the new colnames to qc_df (in case the metric doesn't have
    ## any suffixes, e.g. "rtDuration", leave it as it is)
    colnames_df <- str_remove(colnames(qc_df), pattern = metric)
    colnames_df <- str_remove(colnames_df, pattern = "^_")
    
    colnames_df[colnames_df == ""] <- colnames(qc_df)[colnames_df == ""]
    colnames(qc_df) <- colnames_df
    
    ## add the rownames to the new column rowname
    qc_df <- rownames_to_column(qc_df)
    
    ## convert the table into the long format  
    qc_df_l <- pivot_longer(qc_df, cols = seq_len(ncol(qc_df))[-1])
    qc_df_l$rowname <- factor(qc_df_l$rowname, levels = qc_df$rowname)
    
    return(qc_df_l)
}

#' @name shinyMsQuality
#'
#' @title Shiny application to visualize quality metrics
#'
#' @description
#' The function \code{shinyMsQuality} function starts a shiny application to 
#' visualize the quality metrics interactively. It allows to display all metrics
#' contained in \code{qc}. 
#' 
#' The function accepts the output of \code{calculateMetrics},
#' \code{calculateMetricsFromSpectra}, or \code{calculateMetricsFromMsExperiment}
#' 
#' @details
#' The plots within the shiny application can be saved by clicking on the 
#' download button.
#' 
#' @param qc \code{matrix}, contains the calculated quality metrics, the columns
#' contain the metrics and the rows the samples
#' 
#' @return \code{shiny}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom shiny shinyUI selectInput fluidRow req runApp reactive
#' @importFrom shiny downloadButton downloadHandler
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' @importFrom shinydashboard dashboardBody
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom htmlwidgets saveWidget
#' @importFrom stringr str_split
#' 
#' @examples
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' msexp <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(msexp) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(msexp) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
#'     sampleIndex = c(1, 2), withIndex = c(1, 1))
#' 
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#'     
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange)
#' qc <- calculateMetricsFromMsExperiment(msexp = msexp, metrics = metrics, 
#'     msLevel = 1, relativeTo = "Q1", change = "jump")
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' 
#' if (interactive())
#'     shinyMsQuality(qc = qc)
#'
shinyMsQuality <- function(qc) {
    
    if (!is.matrix(qc)) stop("'qc' is not a matrix")
    if (!is.numeric(qc)) stop("'qc' has to be numeric")

    metrics <- str_split(colnames(qc), pattern = "_", simplify = TRUE)[, 1]
    metrics <- unique(metrics)

    ## define the user interface of the application: create a sidebar that 
    ## allows to select the metrics and a body that displays the plot and 
    ## allows for downloading the plot
    ui <- shinyUI(dashboardPage(skin = "black",
        dashboardHeader(title = "MsQuality"),
        dashboardSidebar(
            selectInput(inputId = "metric",
                               label = "Select metric",
                               choices = metrics, multiple = FALSE)
        ),
        dashboardBody(fluidRow(
            plotlyOutput(outputId = "metricPlot"),
            downloadButton(outputId = "downloadPlot", "Download plot")
        ))
    ))

    server <-  function(input, output, session) {

        ## reactive expression that stores the plotly object
        metricPlot_r <- reactive({
            plotMetric(qc = qc, metric = input$metric)
        })

        ## expression that renders the plotly object
        output$metricPlot <- renderPlotly({
            req(metricPlot_r())
            metricPlot_r()
        })

        ## add functionality to download the plot
        output$downloadPlot <- downloadHandler(
            filename = function() {
                paste("metrics_", input$metric, ".html", sep = "")
            },
            content = function(file) {
                saveWidget(metricPlot_r(), file)
            }
        )
    }

    app <- list(ui = ui, server = server)
    runApp(app)
}
