#' transpose_tibble
#'
#' @param tbl the tibble to transpose
#' @param row.names name of the variable that should be used as new column names
#' @param col_label string that will be used as variable name for the new variable containing old column names
#'
#' @return
#' @export
#'
transpose_tibble = function(tbl,row.names = NA, col_label = "colnames")
{

    if (!is.na(row.names))
    {
        mat = as.matrix(tbl %>% dplyr::select(-row.names))
        row.names = tbl %>% dplyr::select(row.names) %>% pull()
        rownames(mat) = row.names
    }
    else
    {
        mat= as.matrix(tbl)
        #col.names = colnames(tbl)
    }
    col.names = colnames(mat)
    t_tbl = as_tibble(t(mat))
    t_tbl = t_tbl %>% add_column(!!col_label := col.names,.before = 1)
    return(t_tbl)
}
