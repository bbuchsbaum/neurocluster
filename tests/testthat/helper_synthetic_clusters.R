cluster_purity <- function(pred, truth) {
  stopifnot(length(pred) == length(truth))
  valid <- !(is.na(pred))
  if (!any(valid)) return(0)
  tab <- table(truth[valid], pred[valid])
  if (length(tab) == 0) return(0)
  sum(apply(tab, 1, max)) / length(truth[valid])
}
