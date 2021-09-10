
#' @name plotMetric
#'
#' @title Visualize a quality metric
#'
#' @description
#' The function `plotMetric` visualizes the metric values per sample. The 
#' function accepts the output of `calculateMetrics`,
#' `calculateMetricsFromSpectra`, or `calculateMetricsFromMsExperiment` and
#' a vector specifying the metric to display.
#' 
#' @details
#' 
#' 
#' @param qc `matrix`/`data.frame`
#' @param metric `character`
#' 
#' @return `plotly`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_brewer theme_bw
#' @importFrom ggplot2 xlab ggtitle guides guide_legend theme element_text
#' @importFrom ggplot2 element_blank
#' @importFrom plotly ggplotly 
#' 
#' @export
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' mse <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(mse) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(mse) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' mse <- linkSampleData(mse, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' mse <- linkSampleData(mse, with = "experimentFiles.annotations",
#'                       sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(mse) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' qc <- calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     params = params_l)
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' plotMetric(qc, metric = "areaUnderTIC") 
plotMetric <- function(qc, metric = "areaUnderTIC") {
    
    qc_tbl_l <- plotMetric_tibble(qc = qc, metric = metric)
    
    g <- ggplot2::ggplot(qc_tbl_l) +
        ggplot2::geom_point(ggplot2::aes(x = rowname, y = value, col = name)) +
        ggplot2::scale_colour_brewer(palette= "Set1") + ggplot2::theme_bw() +
        ggplot2::xlab("sample") + ggplot2::ggtitle(metric) +
        ggplot2::guides(shape = ggplot2::guide_legend(
            override.aes = list(size = 5))) +
        ggplot2::guides(colour = ggplot2::guide_legend(
            override.aes = list(size= 5))) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, size = 10), 
            panel.grid.major = ggplot2::element_blank(), 
            panel.grid.minor = ggplot2::element_blank())
    
    g |> plotly::ggplotly(tooltip = c("x", "y"))
}

#' @name plotMetric_tibble
#'
#' @title Helper function for plotMetric
#'
#' @description
#' The function `plotMetric_tibble` is a helper function for the function
#' `plotMetric`. It returns a tibble in long format that is interpretable
#' by `ggplot2`.
#' 
#' @details
#' `plotMetric_tibble` will select all columns that start with
#' `metric`. The different levels in the `name` column in the returned tibble 
#' correspond to the columns that were selected and do not contain the
#' `metric` prefix. In case there is no additional specification 
#' (e.g. for the metric `rtDuration` only the column `rtDuration` will 
#' be selected), the `name` column will include the `metric` (`rtDuration`). 
#' 
#' @param qc `data.frame`
#' @param metric `character`
#' 
#' @return `tibble`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom stringr str_remove
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' mse <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(mse) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(mse) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' mse <- linkSampleData(mse, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' mse <- linkSampleData(mse, with = "experimentFiles.annotations",
#'                       sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(mse) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' qc <- calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     params = params_l)
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' plotMetric_df(qc, metric = "areaUnderTIC")
plotMetric_tibble <- function(qc, metric) {
    
    cols <- grep(colnames(qc), pattern = metric)
    
    if (length(cols) == 0) stop("metric not in qc")
    
    qc_df <- as.data.frame(qc)
    qc_df <- qc_df[, cols, drop = FALSE]
    
    ## remove from the colnames the "metric" part and the first _ separator
    ## assign then the new colnames to qc_df (in case the metric doesn't have
    ## any suffixes, e.g. "rtDuration", leave it as it is)
    colnames_df <- stringr::str_remove(colnames(qc_df), pattern = metric)
    colnames_df <- stringr::str_remove(colnames_df, pattern = "^_")
    
    colnames_df[colnames_df == ""] <- colnames(qc_df)[colnames_df == ""]
    colnames(qc_df) <- colnames_df
    
    ## add the rownames to the new column rowname
    qc_df <- tibble::rownames_to_column(qc_df)
    
    ## convert the table into the long format  
    qc_df_l <- tidyr::pivot_longer(qc_df, cols = seq_len(ncol(qc_df))[-1])
    qc_df_l$rowname <- factor(qc_df_l$rowname, levels = qc_df$rowname)
    
    return(qc_df_l)
}

#' @name MsQuality
#'
#' @title Shiny application to visualize quality metrics
#'
#' @description
#' The function `MsQuality` function starts a shiny application to visualize
#' the quality metrics interactively. It allows to display all metrics
#' contained in `qc`. 
#' 
#' The function accepts the output of `calculateMetrics`,
#' `calculateMetricsFromSpectra`, or `calculateMetricsFromMsExperiment`
#' 
#' @details
#' The plots within the shiny application can be saved by clicking on the 
#' download button.
#' 
#' @param qc `data.frame`
#' 
#' @return `shiny`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom shiny shinyUI selectInput fluidRow req runApp
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' @importFrom shinydashboard dashboardBody
#' @importFrom plotly plotlyOutput renderPlotly
#' 
#' @examples
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' mse <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(mse) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(mse) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' mse <- linkSampleData(mse, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' mse <- linkSampleData(mse, with = "experimentFiles.annotations",
#'                       sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(mse) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' qc <- calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     params = params_l)
#' rownames(qc) <- c("Sample 1", "Sample 2")
#' 
#' \dontrun{
#' MsQuality(qc)
#' }
MsQuality <- function(qc) {

    metrics <- stringr::str_split(colnames(qc), 
                                        pattern = "_", simplify = TRUE)[, 1]
    metrics <- unique(metrics)

    ui <- shiny::shinyUI(shinydashboard::dashboardPage(skin = "black",
        shinydashboard::dashboardHeader(title = "MsQuality"),
        shinydashboard::dashboardSidebar(
            shiny::selectInput(inputId = "metric",
                               label = "Select metric",
                               choices = metrics, multiple = FALSE)
        ),
        shinydashboard::dashboardBody(shiny::fluidRow(
            plotly::plotlyOutput(outputId = "metricPlot"),
            shiny::downloadButton(outputId = "downloadPlot", "Download plot")
        ))
    ))

    server <-  function(input, output, session) {

        metricPlot_r <- reactive({plotMetric(qc = qc, metric = input$metric)})

        output$metricPlot <- plotly::renderPlotly({
            shiny::req(metricPlot_r())
            metricPlot_r()
        })

        output$downloadPlot <- shiny::downloadHandler(
            filename = function() {
                paste("metrics_", input$metric, ".html", sep = "")
            },
            content = function(file) {
                htmlwidgets::saveWidget(metricPlot_r(), file)
            }
        )
    }

    app <- list(ui = ui, server = server)
    shiny::runApp(app)
}
