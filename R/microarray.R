#' import a CEL file
#'
#' @param directory path to directory with CEL files
#' @param probe_file provide path to the probe file
#' @param skip_lines number of lines to skip in probe file
#'
#' @return
#' @export
#'
#' @examples
import_CEL = function(directory, probe_file = "../data/GPL570-55999.txt", skip_lines = 16)
{
    # Import CEL files
    CEL_files = list.files(directory,pattern = "\\.CEL$") # files that end with .CEL
    CEL_filepath = paste0(directory,CEL_files)
    dat = affy::ReadAffy(filenames = CEL_filepath)
    eset = rma(dat)
    microarray_expression = eset@assayData$exprs
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
convert_probe_ids = function()
{
    # Convert from probe id to gene symbol
    probe_info = read_delim(probe_file,delim = "\t",skip = skip_lines)

    ID = rownames(microarray_expression)
    microarray_expression = add_column(as_tibble(microarray_expression), ID, .before = 1)

    tbl = as_tibble(probe_info %>% select(ID,`Gene Symbol`) %>% right_join(microarray_expression, by = "ID") ) %>% filter(!is.na(`Gene Symbol`))
    GeneSymbols = tbl$`Gene Symbol`
    sample_id = colnames(tbl)[3:ncol(tbl)]

    microarray_expression = as_tibble(t(tbl %>% select(-ID,-`Gene Symbol`)))
    colnames(microarray_expression) = GeneSymbols
    microarray_expression = add_column(microarray_expression, sample_id, .before = 1)
}


#' Collapse a expression of probes to gene symbols
#'
#' @param accesion_number
#'
#' @return
#' @export
#'
#' @examples
collapse_probe_ids = function(accesion_number)
{
    gse = getGEO(accesion_number)

    probe_data = gse$GSE130588_series_matrix.txt.gz@assayData$exprs %>% # probe expression values
        as_tibble() %>%
        add_column(ID = rownames(gse$GSE130588_series_matrix.txt.gz@assayData$exprs),.before = 1)


    feature_data =
        gse$GSE130588_series_matrix.txt.gz@featureData@data %>%
        as_tibble()

    expression_data =
        feature_data %>%
        filter(!grepl("_[a-z^l]_at",ID)) %>% # filter out probes that match multiple segments
        select(ID, `Gene Symbol` ) %>%
        left_join(probe_data, by = "ID") %>%
        group_by(`Gene Symbol`) %>%
        select(-ID) %>%
        summarise_all(mean)

    return(expression_data)

}



