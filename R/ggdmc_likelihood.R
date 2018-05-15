##' likelihood.dmc for lba_B
##'
##' The likelihood function for LBA type norm, used in DMC.
##'
##' @param p.vector parameter vector
##' @param data data model instance
##' @param min.like minimal likelihood. 1e-10 is the default
##' @return a vector
##' @export
##' @examples
##' model <- ggdmc:::BuildModel(
##' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
##' match.map = list(M = list(s1 = 1, s2 = 2)),
##' factors   = list(S = c("s1", "s2")),
##' constants = c(st0 = 0, sd_v = 1),
##' responses = c("r1", "r2"),
##' type      = "norm")
##'
##' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1, mean_v.false = .25)
##'
##' sim1 <- ggdmc:::simulate.model(model, p.vector, 1024)
##' data <- ggdmc:::BindDataModel(sim1, model)
##'
##' model    <- attr(data, "model")
##' ise      <- attr(data, "cell.empty")
##' allpar   <- attr(model, "all.par")
##' parnames <- attr(model, "par.names")
##' type     <- attr(model, "type")
##' n1idx    <- attr(model, "n1.order")
##' mc       <- attr(model, "match.cell")
##' isr1     <- ggdmc:::check_rd(type, model)
##' cellidx  <- ggdmc:::cellidxmat(data)
##' pnames   <- names(p.vector)
##' matchcell<- attr(model, "match.cell")
##'
##' \dontrun{
##' setwd("/media/yslin/MERLIN/Documents/DMC-PDA/dmc-amsterdam17/")
##' source ("dmc/dmc.R")
##' load_model ("LBA","lba_B.R")
##'
##' den1 <- likelihood_norm(p.vector, data)
##' den2 <- likelihood.dmc(p.vector, data)
##' den3 <- density_norm(p.vector, pnames, allpar, parnames, model, type,
##' dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
##' n1idx, ise, cellidx, data$RT, matchcell, isr1)
##'
##' all.equal(den1, den2)
##' all.equal(den2, den3[,1])
##' ## all TRUEs
##' }
##' @export
likelihood_norm <- function(p.vector, data, min.like = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(p.vector)

  out <- density_norm(p.vector, pnames, allpar, parnames, model, type,
                      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
                      n1idx, ise, cellidx, data$RT, mc, isr1)
  pmax(out[,1], min.like)
}

##' @rdname likelihood_norm
##' @export
likelihood_norm_pda <- function(p.vector, data, min.like = 1e-10)
{
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(p.vector)

  nsim     <- attr(data, "n.pda")
  bw       <- attr(data, "bw")

  out <- density_norm_pda(p.vector, pnames, allpar, parnames, model, type,
                      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
                      n1idx, ise, cellidx, data$RT, mc, isr1, nsim, bw)
  pmax(out[,1], min.like)
}

##' @rdname likelihood_norm
##' @export
likelihood_rd <- function(p.vector, data, min.like = 1e-10) {
  model    <- attr(data,  "model")
  ise      <- attr(data,  "cell.empty")
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  mc       <- attr(model, "match.cell")
  isr1     <- check_rd(type, model)
  cellidx  <- cellIdx2Mat(data)
  pnames   <- names(p.vector)

  out <- density_rd(p.vector, pnames, allpar, parnames, model, type,
                    dimnames(model)[[1]], dimnames(model)[[2]],
                    dimnames(model)[[3]],
                    n1idx, ise, cellidx, data$RT, mc, isr1)
  pmax(out[,1], min.like)
}

#' Remove t0
#'
#' This function minus t0 from RT.
#'
#' @param RT,rt a numeric vector
#' @param t0 a numeric scalar
#' @return a vector
#'
#' @examples
#' rt <- rlnorm(10) + .2
#' dt <- ggdmc::removet0(rt, .2)
#' all.equal(dt, rt - .2)
#' @export
removet0 <- function(RT, t0) {
  if(!is.vector(RT)) stop("RT must be a vector")
  if(!is.numeric(t0)) stop("t0 must be a number")
  out <- remove_t0(RT, t0)
  return(out[,1])
}
