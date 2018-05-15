

#' @rdname rlnr
#' @export
n1PDFfixedt0_lnr <- function (dt, meanlog, sdlog)
{
  if (is.matrix(meanlog) & is.matrix(sdlog)) {
    out <- n1PDFfixedt0_lnr2(dt, meanlog, sdlog)
  } else if (is.vector(meanlog) & is.vector(sdlog)) {
    out <- n1PDFfixedt0_lnr1(dt, meanlog, sdlog)
  } else {
    stop("meanlog/sdlog type not supported.")
  }
  return(out[,1]);
}

