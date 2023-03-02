#' -----------------------------------------------------------------------------
data1 <- as.matrix(read.table("inst/extdata/data1.txt", sep = "", header = F))
colnames(data1) <- NULL
usethis::use_data(data1, overwrite = TRUE)
#' -----------------------------------------------------------------------------
data2 <- as.matrix(read.table("inst/extdata/data2.txt", sep = "", header = F))
colnames(data2) <- NULL
usethis::use_data(data2, overwrite = TRUE)
#' -----------------------------------------------------------------------------

