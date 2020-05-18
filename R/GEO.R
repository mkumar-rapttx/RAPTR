#' Download files associated with a GEO accession number
#'
#' @param accessionNumber GEO accession number of study
#' @param dataDirectory path to directory to place files
#'
#' @return
#' @export
#'
#' @examples
downloadGEO = function(accessionNumber, dataDirectory = NULL)
{
    dataDirectory = if(is.null(dataDirectory)) paste0(accessionNumber,"/") else dataDirectory
    # Import CEL files
    dir.create(dataDirectory,showWarnings = FALSE)
    filePaths = GEOquery::getGEOSuppFiles(accessionNumber,
                                          makeDirectory = FALSE,
                                          baseDir = dataDirectory)
    untar(rownames(filePaths),exdir = dataDirectory)
    #untar(list.files(data_directory,pattern = ".gz$"),exdir = data_directory)
    system(paste0("gunzip ",dataDirectory,"*.gz"))

}

#' Download meta data associated with a GEO accession number
#'
#' @param accessionNumber GEO accession number of study
#'
#' @return
#' @export
#'
#' @examples
getMetaData = function(accessionNumber)
{
    gse = GEOquery::getGEO(accessionNumber)
    series_matrix = gse[[paste0(accessionNumber,"_series_matrix.txt.gz")]]

    meta_data =
        series_matrix@phenoData@data %>% # sample phenotype/meta data
        as_tibble()
}

#' Download probe data associated with a GEO accession number
#'
#' @param accessionNumber GEO accession number of study
#'
#' @return
#' @export
#'
#' @examples
getProbeData = function(accessionNumber)
{
    gse = GEOquery::getGEO(accessionNumber)
    series_matrix = gse[[paste0(accessionNumber,"_series_matrix.txt.gz")]]

    probe_data = series_matrix@assayData$exprs %>% # probe expression values
        as_tibble() %>%
        add_column(ID = rownames(series_matrix@assayData$exprs),.before = 1)
}

#' Download feature data associated with a GEO accession number
#'
#' @param accessionNumber GEO accession number of study
#'
#' @return
#' @export
#'
#' @examples
getFeatureData = function(accessionNumber)
{
    gse = GEOquery::getGEO(accessionNumber)
    series_matrix = gse[[paste0(accessionNumber,"_series_matrix.txt.gz")]]
    feature_data =
        series_matrix@featureData@data %>%
        as_tibble()
}
