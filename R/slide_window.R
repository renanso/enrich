#' Slide window
#'
#' This function will create kmers by sliding a window in the fragment. The window size
#' is 120 bp and the slide_size can be changes as you need. slide_size=1 indicates that a new 120-kmer 
#' will be generated at each nucleotide.
#' @export
slide_window <- function(sequence, window_size, slide_size) 
{
  chars <- strsplit(sequence, "")[[1]]
  windows <- list()
  for (i in seq(from=1, to=(length(chars) - window_size + 1), by = slide_size)) {
    windows[[i]] <- paste(chars[i:(i + window_size - 1)], collapse = "")
  }
  return(windows)
}

