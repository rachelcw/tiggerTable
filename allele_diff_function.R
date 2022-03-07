allele_diff <- function(germs) {
  germs <- lapply(germs, function(x) strsplit(x, '')[[1]])
  View(germs)
  germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
  View(germs_m)
  setdiff_mat <- function(x) {
    sum(!unique(x) %in% c('.', NA, "N"))
  }
  idx = which(apply(germs_m, 2, setdiff_mat) > 1)
  return(idx)
}