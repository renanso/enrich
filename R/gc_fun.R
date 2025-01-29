#' GC % calculation
#'
#' This function is a simple function to calculate GC content in the probe
#' @export
gc_fun <- function(x){ 
  num_g <- str_count(x, "G")
  num_c <- str_count(x, "C")
  return ((num_g + num_c) / str_length(x) * 100 )} 