#' Split a string at every position of a delimiter, and returns the subset at a specified index.
#'
#' @param string_vector the string to split
#' @param return_position the index of the split string to return
#' @param sep_char the delimiter to use to split the string
#' @param fixed logical passed to strsplit. If TRUE match split exactly, otherwise use regular expressions.
#' @export
split_string = function(string_vector, return_position = 1, sep_char = "_", fixed = TRUE)
{
    split_sample = strsplit(as.character(string_vector),split = sep_char, fixed = fixed)
    split_val = sapply(split_sample,`[`,return_position)
    return(split_val)
}
