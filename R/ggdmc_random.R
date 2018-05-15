


##' Piecewise LBA Model Type 1 and Type 2
##'
##' Density and random generation for the PLBA Model Type 1 and Type 2.
##'
##' @param n number of observations.
##' @param x vector of quantiles.
##' @param A upper bound of start point. Can be an integer a 2-element vector.
##' @param b threshold. Can be an integer a 2-element vector.
##' @param mean_v piece 1 mean drift rate. Should be a 2-element vector
##' @param mean_w piece 2 mean drift rate. Should be a 2-element vector
##' @param sd_v common standard deviation of the piece 1 drift rates. If
##' sd_w is not present, this will be used as the piece 2 drift rate sd,
##' too. Negative value is not allowed.
##' @param sd_w standard deviation of the piece 2 drift rates for
##' @param rD rate delay
##' @param swt switch time
##' @param t0 nondecision time
##' @param h bandwidth
##' @param ncore number of CPU cores for running Open MP.
##' @param debug internal debug switch
##' @return a [RT R] matrix (C) or a data frame (R)
##' @examples
##' #############20
##' ## rplba1    ##
##' #############20
##' \dontrun{
##' n <- 2^20; n
##' A <- 1.5
##' b <- 2.7
##' mean_v <- c(3.3, 2.2)
##' mean_w <- c(1.5, 3.7)
##' sd_v <- c(1, 1)
##' rD  <- .3
##' swt <- .5
##' t0  <- .08
##' ncore <- 12
##' dat1 <- ggdmc:::rplba1R(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
##' dat2 <- ggdmc:::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt, ncore)
##' dat3 <- gpda:::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
##' dat4 <- ggdmc:::rplba1_test(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
##'
##' dat1r1 <- dat1[dat1[, 2] == 1, 1]
##' dat1r2 <- dat1[dat1[, 2] == 2, 1]
##' dat2r1 <- dat2[dat2[, 2] == 1, 1]
##' dat2r2 <- dat2[dat2[, 2] == 2, 1]
##' dat3r1 <- dat3[dat3[, 2] == 1, 1]
##' dat3r2 <- dat3[dat3[, 2] == 2, 1]
##' dat4r1 <- dat4[dat4[, 2] == 1, 1]
##' dat4r2 <- dat4[dat4[, 2] == 2, 1]
##'
##' xlim <- c(0, 3)
##' ## Check if two methods produce SPDF overlaping with each other
##' par(mfrow = c(4, 2), mar = c(4, 5.3, 0.82, 1))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R", xlim = xlim)
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R", xlim = xlim)
##' hist(dat2r1, breaks = "fd", freq = FALSE, main = "Choice1 C++", xlim = xlim)
##' hist(dat2r2, breaks = "fd", freq = FALSE, main = "Choice2 C++", xlim = xlim)
##' hist(dat3r1, breaks = "fd", freq = FALSE, main = "Choice1 GPU", xlim = xlim)
##' hist(dat3r2, breaks = "fd", freq = FALSE, main = "Choice2 GPU", xlim = xlim)
##' hist(dat4r1, breaks = "fd", freq = FALSE, main = "Choice1 test", xlim = xlim)
##' hist(dat4r2, breaks = "fd", freq = FALSE, main = "Choice2 test", xlim = xlim)
##'
##' par(mfrow = c(1, 2))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R, C++, & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r1, breaks = "fd", freq = FALSE, add = TRUE, col = "lightblue")
##' hist(dat3r1, breaks = "fd", freq = FALSE, add = TRUE, col = "lightgreen")
##'
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R, C++ & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r2, breaks = "fd", freq = FALSE, add = TRUE, col = "lightblue")
##' hist(dat3r2, breaks = "fd", freq = FALSE, add = TRUE, col = "lightgreen")
##'
##'
##' library(rbenchmark)
##' res <- benchmark(r1 = ggdmc:::rplba1R(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt),
##'   r2 = ggdmc:::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt),
##'   r3 = gpda:::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt),
##'   replications = 10)
##'
##' print(res[,1:4])
##'
##' # test replications elapsed relative
##' #   r1           10   2.830  134.762
##' #   r2           10   0.435   20.714
##' #   r3           10   0.021    1.000
##' #  --------------------------------#
##' #   r1           10  19.356  186.115
##' #   r2           10   2.946   28.327
##' #   r3           10   0.104    1.000
##'
##' res <- benchmark(r1 = ggdmc:::rplba1_test(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt),
##'                  r2 = ggdmc:::rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt, core),
##'                  replications = 10)
##'
##' print(res[,1:4])
##' ## test replications elapsed relative
##' ##   r1           10   3.484    1.177
##' ##   r2           10   2.959    1.000
##' }
##'
##' #############20
##' ## rplba2    ##
##' #############20
##' \dontrun{
##' n <- 2^15
##' ncore <- 4
##' A <- c(1.5, 1.5)
##' b <- c(2.7, 2.7)
##' mean_v <- c(3.3, 2.2)
##' mean_w <- c(1.5, 3.7)
##' sd_v <- c(1, 1)
##' sd_w <- c(1, 1)
##' rD <- .3
##' swt <- .5
##' t0 <- .08
##' dat1 <- ggdmc:::rplba2R(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##' dat2 <- ggdmc:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt, ncore)
##' dat3 <- gpda:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##' dat4 <- ggdmc:::rplba2_test(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
##'
##' dat1r1 <- dat1[dat1[, 2] == 1, 1]
##' dat1r2 <- dat1[dat1[, 2] == 2, 1]
##' dat2r1 <- dat2[dat2[, 2] == 1, 1]
##' dat2r2 <- dat2[dat2[, 2] == 2, 1]
##' dat3r1 <- dat3[dat3[, 2] == 1, 1]
##' dat3r2 <- dat3[dat3[, 2] == 2, 1]
##' dat4r1 <- dat4[dat4[, 2] == 1, 1]
##' dat4r2 <- dat4[dat4[, 2] == 2, 1]
##'
##' wesanderson::wes_palette("Royal1")
##' palettes  <- wesanderson::wes_palettes$GrandBudapest
##' palettes2 <- wesanderson::wes_palettes$GrandBudapest2
##' xlim <- c(0, 3)
##' ## Check if two methods produce SPDF overlaping with each other
##' par(mfrow = c(4, 2), mar = c(4, 5.3, 0.82, 1))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R", xlim = xlim)
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R", xlim = xlim)
##' hist(dat2r1, breaks = "fd", freq = FALSE, main = "Choice1 C++", xlim = xlim)
##' hist(dat2r2, breaks = "fd", freq = FALSE, main = "Choice2 C++", xlim = xlim)
##' hist(dat3r1, breaks = "fd", freq = FALSE, main = "Choice1 GPU", xlim = xlim)
##' hist(dat3r2, breaks = "fd", freq = FALSE, main = "Choice2 GPU", xlim = xlim)
##' hist(dat4r1, breaks = "fd", freq = FALSE, main = "Choice1 test", xlim = xlim)
##' hist(dat4r2, breaks = "fd", freq = FALSE, main = "Choice2 test", xlim = xlim)
##'
##' par(mfrow = c(1, 2))
##' hist(dat1r1, breaks = "fd", freq = FALSE, main = "Choice1 R, C++, & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[1])
##' hist(dat3r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[2])
##' hist(dat4r1, breaks = "fd", freq = FALSE, add = TRUE, col = palettes[4])
##'
##' hist(dat1r2, breaks = "fd", freq = FALSE, main = "Choice2 R, C++ & GPU",
##'   xlim = xlim, ylim = c(0, 3))
##' hist(dat2r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[1])
##' hist(dat3r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[2])
##' hist(dat4r2, breaks = "fd", freq = FALSE, add = TRUE, col = palettes2[3])
##'
##' library(rbenchmark)
##' res <- benchmark(r1 = ggdmc:::rplba2R(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt),
##'   r2 = ggdmc:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt, ncore),
##'   r3 = gpda:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt),
##'   replications = 10)
##'
##' print(res[,1:4])
##' ## test replications elapsed relative
##' ##   r1           10   0.480   30.000
##' ##   r2           10   0.077    4.813
##' ##   r3           10   0.016    1.000
##'
##' library(rbenchmark)
##' res <- benchmark(r1 = ggdmc:::rplba2R(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt),
##'                  r2 = ggdmc:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt),
##'                  replications = 10)
##'
##' print(res[,1:4])
##' ## test replications elapsed relative
##' ##   r1           10  25.454    8.415
##' ##   r2           10   3.025    1.000
##'
##' library(rbenchmark)
##' n <- 2^23
##' ncore <- 4
##' res <- benchmark(r1 = ggdmc:::rplba2_test(n, A, b, t0, mean_v, mean_w, sd_v,
##'                                           sd_w, rD, swt),
##'                  r2 = ggdmc:::rplba2(n, A, b, t0, mean_v, mean_w, sd_v,
##'                  sd_w, rD, swt, ncore, FALSE),
##'                  replications = 10)
##'
##' print(res[,1:4])
##' ## test replications elapsed relative
##' ##   r1           10  26.999    1.147
##' ##   r2           10  23.548    1.000
##' }
##' @export
rplba1R <- function(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt)
{
  if (length(mean_v) != 2) stop("Current version only provides 2-accumulator model.")
  if (length(sd_v) != 2) stop("The length of sd_v and mean_v must match.")
  eswt <- swt + rD
  v1 <- rtnorm(n, mean_v[1], sd_v[1], 0, Inf)[,1] ## Stage 1 LBA
  v2 <- rtnorm(n, mean_v[2], sd_v[2], 0, Inf)[,1]
  sp <- matrix(runif(2*n, 0, A), 2)
  dt1 <- rbind((b - sp[1,])/v1, (b - sp[2,])/v2) ## Race

  ## dt[dt<0] <- Inf
  choice    <- apply(dt1, 2, which.min)
  chosen_dt <- dt1[cbind(choice, 1:n)]  ## convert to vector choose (row, col)

  done <- (chosen_dt <= eswt)   ## Which are finished?
  n2 <- sum(!done)

  ## Distance left to travel for those not finished
  B1 <- b - (sp[1, !done] + eswt*v1[!done])
  B2 <- b - (sp[2, !done] + eswt*v2[!done])

  w1 <- rtnorm(n2, mean_w[1], sd_v[1], 0, Inf)[,1]   ## Stage 2 LBA
  w2 <- rtnorm(n2, mean_w[2], sd_v[2], 0, Inf)[,1]
  dt2 <- rbind(B1/w1, B2/w2)   ## Race

  choice[!done] <- apply(dt2, 2, which.min)
  chosen_dt[!done] <- eswt + dt2[cbind(choice[!done], 1:n2)]

  ## The last object automatically return
  data.frame(cbind(RT=t0 + chosen_dt, R = choice))
}

#' @rdname rplba1R
#' @export
rplba2R <- function(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt)
{
  # Calculate effective switch time
  eswt <- swt + rD

  # Stage 1 LBA Race
  v <- t(cbind(rtnorm(n, mean_v[1], sd_v[1], 0, Inf), rtnorm(n, mean_v[2], sd_v[2], 0, Inf)))
  sp <- matrix(runif(2*n, 0, A), 2)
  dt <- (b- sp) / v
  # dt[dt<0] <- Inf
  choice <- apply(dt, 2, which.min)
  chosen_dt <- dt[cbind(choice, 1:n)]

  # Which are finished?
  done <- (chosen_dt <= eswt)
  n2   <- sum(!done)

  # Distance left to travel for those not finished
  B <- b - (sp[, !done] + eswt * v[, !done])

  # Stage 2 LBA Race
  w <- t(cbind(rtnorm(n2, mean_w[1], sd_w[1], 0, Inf), rtnorm(n2, mean_w[2], sd_w[2], 0, Inf)))
  dt <- B / w
  choice[!done] <- apply(dt, 2, which.min)
  chosen_dt[!done] <- eswt + dt[ cbind(choice[!done], 1:n2) ]

  # save results
  ## The last object automatically return
  data.frame(cbind(RT = t0 + chosen_dt, R = choice))
}

#' @rdname rplba1R
#' @export
rplbaR3 <- function(n=10, pVec=c(A1=1.5, A2=1.5, B1=1.2, B2=1.2, C1=.3, C2=.3,
  v1=3.32, v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1,
  sw1=1, sw2=1, rD=0.3, tD=.3, swt=0.5, t0=0.08)) {

  # n <- 10
  # pVec=c(A1=1.5, A2=1.5, B1=1.2, B2=1.2, C1=.3, C2=.3,
  #   v1=3.32, v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1,
  #   sw1=1, sw2=1, rD=0.3, tD=.3, swt=0.5, t0=0.08)

  # Stage 1 LBA
  v1 <- rtnorm(n, pVec["v1"], pVec["sv1"], 0, Inf)[,1]
  v2 <- rtnorm(n, pVec["v2"], pVec["sv2"], 0, Inf)[,1]
  sp <- matrix(runif(2*n,0,pVec[c("A1","A2")]),nrow=2)

  # Calcualte thresholds
  b1 <- sum(pVec[c("A1","B1")])
  b2 <- sum(pVec[c("A2","B2")])
  c1 <- b1 + pVec[c("C1")]
  c2 <- b2 + pVec[c("C2")]

  # Race
  dt <- rbind((c(b1,b2)-sp[1,])/v1,(c(b1,b2)-sp[2,])/v2)
  # dt[dt<0] <- Inf
  choice <- apply(dt,2,which.min)
  rt <- dt[cbind(choice,1:n)]

  # Calculate effective switch times
  swt_b <- pVec["swt"] + pVec["tD"]
  swt_r <- pVec["swt"] + pVec["rD"]

  # Which switch is first
  swt <- pmin(swt_b,swt_r)
  if (swt_b==swt_r) {
    change <- "both"
  } else if (swt_r < swt_b) {
    change <- "rate"
  } else {
    change <- "threshold"
  }

  # Which are finished?
  done <- rt <= swt
  n2 <- sum(!done)

  # Stage 2 LBA

  # Distance left to travel for those not finished
  # threshold - distance already travelled
  if ( change=="rate" ) {
    B1 <- b1 - (sp[1,!done] + swt*v1[!done])
    B2 <- b2 - (sp[2,!done] + swt*v2[!done])
  } else {
    B1 <- c1 - (sp[1,!done] + swt*v1[!done])
    B2 <- c2 - (sp[2,!done] + swt*v2[!done])
  }


  # Change rates?
  if ( change=="threshold" ) {
    w1 <- v1[!done]; w2 <- v2[!done]
  } else {
    w1 <- rtnorm(n2, pVec["w1"], pVec["sw1"],0, Inf)[,1]
    w2 <- rtnorm(n2, pVec["w2"], pVec["sw2"],0, Inf)[,1]
  }

  # Race
  dt <- rbind(B1/w1,B2/w2)
  # dt[dt<0] <- Inf
  choice[!done] <- apply(dt,2,which.min)
  rt[!done] <- swt+dt[cbind(choice[!done],1:n2)]

  if ( change != "both" ) { # Stage 3 LBA

    if ( change=="threshold" ) swt1 <- swt_r else swt1 <- swt_b
    t2 <- swt1-swt

    # Which are finished?
    done1 <- rt[!done] < swt1
    n2 <- sum(!done1)

    if ( !all(done1) ) {

      # Distance left to travel for those not finished
      # Distance left at end of stage 1 - further travel
      B1 <- B1[!done1] - t2*w1[!done1]
      B2 <- B2[!done1] - t2*w2[!done1]

      if ( change=="threshold" ) {
        w1 <- rtnorm(n2,pVec["w1"],pVec["sw1"],0, Inf)[,1]
        w2 <- rtnorm(n2,pVec["w2"],pVec["sw2"],0, Inf)[,1]
      }  else {
        w1 <- w1[!done1];
        w2 <- w2[!done1]
        B1 <- B1 + pVec["C1"]
        B2 <- B2 + pVec["C2"]
      }

      # Race
      dt <- rbind(B1/w1,B2/w2)
      # dt[dt<0] <- Inf
      choice[!done][!done1] <- apply(dt,2,which.min)
      rt[!done][!done1] <- swt1+dt[cbind(choice[!done][!done1],1:n2)]
    }

  }

  # save results
  data.frame(R=choice, RT=pVec["t0"] + rt)
}

#' @rdname rlnr
#' @export
rlnrR <- function (n, meanlog, sdlog, t0, st0 = 0)
{
  n_acc <- ifelse(is.null(dim(meanlog)), length(meanlog), dim(meanlog)[1])
  dt    <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog),
               nrow = n_acc) + t0
  winner <- apply(dt, 2, which.min)
  if (st0[1] == 0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
    data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=winner)
}


#' Canonical Linear Ballistic Accumulation/Accumualtor Model
#'
#' \code{makeR} stands for making/generating/simulating responses from
#' a LBA model. \code{make_r} and \code{make.r} use C++ function. These
#' make \code{r}, \code{_r}, \code{.r} functions are essentially \code{rLBA},
#' including \code{rlba_norm}. They uses a LBA model with parameters, b, A,
#' mean_v, sd_v and t0 (no st0) to generate choice RT random deviates.
#'
#' \code{make_v} draws drift rate from normal or truncated normal distribution.
#' Each trial is stored as a row and each column is a drift rate for an
#' accumulator. You need to transpose drift rates generated by make_v for
#' \code{makeR}.
#'
#' \code{make.r} is a wrapper function of \code{make_r}. You may
#' need to use ":::" to call make.r, because of S3 method naming convention. If
#' you call \code{make_r} directly, beware it returns C index and is only a
#' numeric matrix. It does not carry a string vector for the column names, RTs
#' and responses. See timing test to see why it might be a good idea not to
#' return it as a data frame. \code{rlbaCnorm} is R version of correlated LBA
#' model.
#'
#' \code{rlba_norm} adds checks and calls \code{make_v} and \code{make_r}.
#' \code{rlba_norm} is only slightly quicker than \code{make_r}.
#'
#' \code{n1PDFfixedt0} is defective density function for the fisrt node LBA
#' model. Defective means its probability does not normally normalize to 1.
#' Only the probabilities from all nodes/accumulators add to 1.
#' \code{n1PDFfixedt0} is equation (3) on  page 159 in Brown and
#' Heathcote (2008).  This equation assumes st0 is 0.
#'
#' \code{fptcdf} and \code{fptpdf} are distribution and density functions with
#' four parameters A, b, v and sv, assuming t0 is zero. \code{fptcdf} and
#' \code{fptpdf} are respectively equation (1) and equation (2) on page 159 in
#' Brown and Heathcote (2008).
#'
#' @param rt response times. A vector.
#' @param drifts a n x n_v drift rate matrix. It can be a vector with 2 or more
#' elements. n is the numbers of observation. n_v is the numbers of
#' response/accumulator.
#' @param b decision threshold, a vector or a scalar.
#' @param A start point upper bound, a vector of a scalar.
#' @param mean_v mean drift rate. Must be a n_acc-element vector
#' @param sd_v standard deviation of the drift rate. a scalar or a vector with
#' n_v elements.
#' @param h bandwidth. Only for PDA.
#' @param n_v numbers of response/accumulator, an integer. Note n_v must match
#' the length/size of \code{drifts} vector.
#' @param t0 nondecision time, a vector or a scalar.
#' @param st0 nondecision time variation, a vector of a scalar. It is the upper
#' bound of a uniform distribution for t0 variability.
#' @param n numbers of observation/model simulations. This must be a scalar.
#' @param return.ttf a boolean switch indicating if return RTs for all
#' accumulators. When \code{return.ttf} is TRUE, a n_v x n ttf matrix is
#' returned.
#' @param posdrift a boolean switch indicating if trimming off negative drift
#' rates when drawing random drift rates.
#'
#' @return \code{make_r} gives either a time-to-finish (ttf) matrix or a n x 2
#' matrix, storing RTs (first column) and responses (second column). \code{n}
#' equals to number of model simulations. ttf is a n_v x n matrix with RTs from
#' all accumulators.
#
#' @examples
#' ## Basic test ----
#' mean_v <- matrix(c(2.4, 2.2)); mean_v
#' n <- 10
#' n_v <- 2
#' A <- 1.2
#' b <- 2.7
#' t0 <- .2
#' sd_v <- 1
#' st0 <- 0
#' posdrift <- TRUE
#' return.ttf <- FALSE
#'
#' drifts <- ggdmc::make_v(n, A, b, t0, mean_v, sd_v, st0, posdrift)
#' dat0 <- ggdmc::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#' dat1 <- ggdmc:::make.r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#' dat1 <- ggdmc::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#' dat1 <- ggdmc::makeR(t(drifts), b, A, n_v, t0, st0, n, return.ttf)
#' dat2 <- ggdmc::rlba_norm(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
#'
#' ## return time to finish ----
#' return.ttf <- TRUE
#' r1 <- ggdmc:::make.r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#' r2 <- ggdmc::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#'
#' ## negative mean drift rates  ----
#' return.ttf <- FALSE
#' posdrift <- FALSE
#' mean_v <- matrix(c(2.4, -2.2)); mean_v
#' drifts <- ggdmc::make_v(n, A, b, t0, mean_v, sd_v, st0, posdrift)
#'
#' ## 3 accumulators ----
#' mean_v <- matrix(c(2.4, 2.3, 2.6)); mean_v
#' drifts <- ggdmc::make_v(n, A, b, t0, mean_v, sd_v, st0, posdrift)
#'
#' posdrift <- TRUE
#' n_v <- 3
#' driftsR <- t(drifts)
#'
#' r1 <- ggdmc:::make.r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#' r2 <- ggdmc::makeR(driftsR, b, A, n_v, t0, st0, n, return.ttf)
#' r3 <- ggdmc::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
#'
#' ## benchmark ----
#' \dontrun{
#' library(rbenchmark)
#' res <- benchmark(r1 = ggdmc::makeR(driftsR, b, A, n_v, t0, st0, n, return.ttf),
#'                  r2 = ggdmc:::make.r(drifts, b, A, n_v, t0, st0, n, return.ttf),
#'                  r3 = ggdmc::make_r(drifts, b, A, n_v, t0, st0, n, return.ttf),
#'                  replications = 100)
#' print(res[,1:4])
#' ## test replications elapsed relative
#' ##   r1          100   0.023     11.5  ## data.frame plus R overhead
#' ##   r2          100   0.019      9.5  ## data.frame time
#' ##   r3          100   0.002      1.0  ## optimal speed
#' }
#' ## rlba_norm ----
#' p.vector  <- c(A= .75, B=.25, t0=.2, mean_v.true=2.5, mean_v.false= 1.5)
#' #'
#' n <- 2^10
#' mean_v <- matrix(c(2.5, 1.5)); mean_v
#' n_v <- 2
#' A <- .75
#' b <- 1
#' t0 <- .2
#' sd_v <- 1
#' st0 <- 0
#' posdrift <- TRUE
#' return.ttf <- FALSE
#'
#' dat0 <- clba::rlba_norm(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
#' dat1 <- clba:::rlba.norm(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
#' head(dat0)
#' head(dat1)
#' ylba1 <- sort(dat0[dat0[,2] == 0, 1]) ## rt1
#' ylba2 <- sort(dat0[dat0[,2] == 1, 1]) ## rt2
#'
#' \dontrun{
#' den0 <- rtdists::dLBA(dat0[, 1], dat0[, 2]+1, A=A, b=b, t0=t0, mean_v=mean_v[,1],
#'                       sd_v=rep(sd_v, 2))
#'
#' df0 <- cbind(dat0, den0)
#' head(df0)
#' df1 <- df0[df0[,2]==0,]
#' df2 <- df0[df0[,2]==1,]
#' head(df1)
#' head(df2)
#' denlba1 <- df1[order(df1[,1]),3]
#' denlba2 <- df2[order(df2[,1]),3]
#' lbs <- 2
#' tks <- 1.5
#' par(mar = c(4, 5.3, 0.82, 1))
#' plot(ylba1, denlba1, type = "l", xlab="RT", ylab="Density", cex.lab=lbs,
#'      cex.axis=tks, lwd = 2, bty = "n")
#' lines(ylba2, denlba2, lwd=2)
#' text(0.6, 2.6, "Choice 1",  cex = tks)
#' text(.7, .8, "Choice 2",  cex = tks)
#' text(1.0, 2.0, "LBA Model", cex = tks, pos=4)
#' }
#'
#' ## make_v example
#' mean_v <- matrix(c(2.4, 2.2)); mean_v
#' n    <- 10
#' n_v  <- 2
#' A    <- 1.2
#' b    <- 2.7
#' t0   <- .2
#' sd_v <- 1
#' st0  <- 0
#' posdrift   <- TRUE
#' return.ttf <- FALSE
#' make_v(n, A, b, t0, mean_v, sd_v, st0, posdrift, return.ttf)
#' ##          [,1]     [,2]
#' ## [1,] 3.322569 3.579262
#' ## [2,] 1.657905 1.181219
#' ## [3,] 2.493664 1.093876
#' ## ...
#'
#' x <- seq(0, 5, .01)
#' head(rtdists:::n1PDFfixedt0)
#' head(ggdmc::n1PDFfixedt0)
#' den1 <- rtdists::n1PDF(x, A = .25, b = .6, mean_v = c(1, .25), sd_v = c(1, 2), t0 = .2, silent = TRUE)
#' den2 <- ggdmc::n1PDFfixedt0(x, A = .25, b = .6, mean_v = c(1, .25), sd_v = c(1, 2), t0 = .2)
#' den3 <- ggdmc::n1PDFfixedt0(x, A = rep(.25, 2), b = rep(.6, 2), mean_v = c(1, .25), sd_v = c(1, 2), t0 = rep(.2, 2))
#'
#' all.equal(den1, den2[,1])
#' all.equal(den1, den3[,1])
#' ## all TRUE
#'
#' \dontrun{
#' require(rbenchmark)
#' res <- benchmark(r1 = rtdists::n1PDF(x, A = .25, b = .6, mean_v = c(1, .25), sd_v = c(1, 2), t0 = .2, silent = TRUE),
#'                  r2 = ggdmc::n1PDFfixedt0(x, A = .25, b = .6, mean_v = c(1, .25), sd_v = c(1, 2), t0 = .2),
#'                  r3 = ggdmc::n1PDFfixedt0(x, A = rep(.25, 2), b = rep(.6, 2), mean_v = c(1, .25), sd_v = c(1, 2), t0 = rep(.2, 2)),
#'                  replications = 100)
#'
#' print(res[,1:4])
#' ## test replications elapsed relative
#' ##   r1          100   0.068    4.000
#' ##   r2          100   0.017    1.000
#' ##   r3          100   0.021    1.235
#' }
#'
#' data <- seq(0, 3, length.out = 1e3);
#' den1 <- ggdmc::n1PDFfixedt0_pda(data, nsim=2^18, 1, .5, c(2.4, 1.6), c(1,1), .5,
#'       h_in = 0.001, k = .09, debug=T)
#' plot(data, den1, type="l")
#'
#' ######################
#' ## n1PDFfixedt0_pda
#' ######################
#' mean_v <- matrix(c(2.4, 2.2)); mean_v
#' n <- 2^14
#' n_v <- 2
#' A <- 1.2
#' b <- 2.7
#' t0 <- .2
#' sd_v <- c(1, 1)
#' h <- .01
#'
#' x <- seq(0, 3, .001);
#' den1 <- ggdmc::n1PDFfixedt0_pda(x, A, b, mean_v, sd_v, t0, n, h)
#' den2 <- ggdmc::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0)
#' den3 <- rtdists::n1PDF(x, A, b, mean_v = as.vector(mean_v), sd_v = sd_v, t0)
#' plot(x, den1, type="l")
#' lines(x, den2)
#' lines(x, den3, lwd= 3)
#'
#' library(rbenchmark)
#' res <- benchmark(r1 = ggdmc::n1PDFfixedt0_pda(x, A, b, mean_v, sd_v, t0, n, h),
#'                  r2 = rtdists::n1PDF(x, A, b, mean_v = as.vector(mean_v), sd_v = sd_v, t0, silent = TRUE),
#'                  r3 = ggdmc::n1PDFfixedt0(x, A, b, mean_v, sd_v, t0),
#'                  replications = 100)
#' print(res[,1:4])
#'
#' ## test replications elapsed relative
#' ##   r1          100   0.512    4.876
#' ##   r2          100   0.257    2.448
#' ##   r3          100   0.105    1.000
#'
#' ## SB's CDF R script
#' ## fptcdf <- function(z,x0max,chi,v,sdv) {
#' ##   if (x0max==0) return(pnorm(chi/z,mean=v,sd=sdv,lower.tail=F))
#' ##   zs=z*sdv ; zu=z*v ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
#' ##   chizu=chiminuszu/zs ; chizumax=xx/zs
#' ##   tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
#' ##   tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
#' ##   1+(tmp1+tmp2)/x0max
#' ## }
#'
#' dt   <- seq(.2, 2, by=0.05)          ## Decision times
#' pVec <- c(A=.5, b=1, v=2.4, sv=1)    ## parameter vector
#' Y    <- fptcdf(dt=dt, A=pVec[1], b=pVec[2], v=pVec[3], sv=pVec[4])
#'
#' dt   <- seq(.1, 2, 0.01)           ## Decision times
#' pVec <- c(A=.5, b=1, v=2.4, sv=1)  ## parameter vector
#'
#' Y1 <- fptcdf(z = dt, x0max = pVec[1], chi = pVec[2], v = pVec[3], sdv=pVec[4])
#' Y2 <- ggdmc:::fptcdf(rt = dt, A = pVec[1], b = pVec[2], mean_v = pVec[3], sd_v=pVec[4], t0 = 0, posdrift = FALSE)
#' Y3 <- rtdists:::plba_norm_core(dt, pVec[1], pVec[2], 0, pVec[3], pVec[4], FALSE, posdrift =  FALSE, length(dt))
#'
#' all.equal(Y1, Y2[,1])
#' all.equal(Y3, Y2[,1])
#' all.equal(Y1, Y3)
#' ## all TRUE
#'
#' \dontrun{
#' ## res <- benchmark(
#' ##   r1 = fptcdf(z=dt, x0max=pVec[1], chi=pVec[2], v = pVec[3], sdv=pVec[4]),
#' ##   r2 = ggdmc::fptcdf(rt=dt, A=pVec[1], b=pVec[2], mean_v = pVec[3], sd_v=pVec[4], 0),
#' ##   r3 = rtdists:::plba_norm_core(dt, pVec[1], pVec[2], 0, pVec[3], pVec[4], TRUE, FALSE, length(dt)),
#' ##   replications = 5e3)
#' ## print(res[,1:4])
#'
#' ## test replications elapsed relative
#' ## r1         5000   0.296    1.003
#' ## r2         5000   0.295    1.000
#' ## r3         5000   0.528    1.790
#' }
#'
#' ## SB's PDF R script
#' ## fptpdf <- function(z,x0max,chi,v,sdv) {
#' ## if (x0max==0) return( (chi/z^2)*dnorm(chi/z,mean=v,sd=sdv))
#' ##   zs=z*sdv ; zu=z*v ; chiminuszu=chi-zu
#' ##     chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
#' ##       (v*(pnorm(chizu)-pnorm(chizumax)) +
#' ##         sdv*(dnorm(chizumax)-dnorm(chizu)))/x0max
#' ## }
#'
#' y1 <- fptpdf(z=dt, x0max=pVec[1], chi=pVec[2], v = pVec[3], sdv=pVec[4])
#' y2 <- ggdmc::fptpdf(rt=dt, A=pVec[1], b=pVec[2], mean_v = pVec[3], sd_v=pVec[4], t0 = 0, FALSE)
#' y3 <- rtdists:::dlba_norm_core(dt, pVec[1], pVec[2], 0, pVec[3], pVec[4], FALSE, FALSE, length(dt))
#'
#' all.equal(y1, y2[,1])
#' all.equal(y3, y2[,1])
#' ## all TRUE
#'
#' \dontrun{
#' res <- benchmark(
#'      r1 = fptpdf(z=dt, x0max=pVec[1], chi=pVec[2], v = pVec[3], sdv=pVec[4]),
#'      r2 = ggdmc::fptpdf(rt=dt, A=pVec[1], b=pVec[2], mean_v = pVec[3], sd_v=pVec[4], 0),
#'      r3 = rtdists:::dlba_norm_core(dt, pVec[1], pVec[2], 0, pVec[3], pVec[4], TRUE, FALSE, length(dt)),
#'      replications = 5e3)
#' print(res[,1:4])
#'
#' ## test replications elapsed relative
#' ##   r1         5000   0.296    1.138
#' ##   r2         5000   0.260    1.000
#' ##   r3         5000   0.510    1.962
#' }
#' @export
makeR <- function (drifts, b, A, n_v, t0, st0 = 0, n, return.ttf = FALSE)
{
  drifts[drifts < 0] <- 0
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  ttf <- t0 + (b - starts)/drifts
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  rt <- ttf[cbind(resp,1:n)]
  if (st0[1]>0) rt <- rt + runif(min = 0, max = st0[1], n = n)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad), "infinite RTs removed and less than",
                  n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(RT = rt, R = resp)
}

#' @rdname makeR
make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n, return.ttf = FALSE)
{
  tmp <- make_r(drifts, b, A, n_v, t0, st0, n, return.ttf)
  if (return.ttf) {
    return(tmp)
  } else {
    return(data.frame(RT = tmp[,1], R = tmp[,2]))
  }
}


#' @importFrom rtdists rdiffusion
#' @export
random <- function(type, pmat, n, p.prior = NA)
{
  if (type == "rd") {
    out <- rtdists::rdiffusion(n, a = pmat$a[1], v = pmat$v[1],
      t0 = pmat$t0[1],
      z  = pmat$z[1]*pmat$a[1], # convert to absolute
      d  = pmat$d[1],
      sz = pmat$sz[1]*pmat$a[1],
      sv = pmat$sv[1], st0 = pmat$st0[1], stop_on_error = TRUE)

  } else if (type %in% c("norm", "norm_pda", "norm_pda_gpu")) {

    ## pmat: A b t0 mean_v sd_v st0
    out <- rlba_norm(n, pmat[, 1], pmat[, 2], matrix(pmat[, 4]), pmat[, 5],
      pmat[,3], pmat[1,6])

  } else if (type %in% c("plba0_gpu") ) {

    out <- rplba0(n, pmat[,1], pmat[,2], pmat[1,7], pmat[,3], pmat[,5],
      pmat[, 4], pmat[1,6], pmat[1,8])

  } else if (type %in% c("plba1", "plba1_gpu") ) {

    out <- rplba1(n, pmat[,1], pmat[,2], pmat[1,7], pmat[,3], pmat[,5],
      pmat[, 4], pmat[1,6], pmat[1,8])

  } else if (type == "plba2") {
    ## A   b mean_v mean_w sd_v sd_w  rD  t0 swt
    # arma::mat rplba2(int n, arma::vec A, arma::vec b,
    #   arma::vec mean_v, arma::vec mean_w, arma::vec sd_v, arma::vec sd_w,
    #   double rD, double swt, double t0)
    out <- rplba2(n, pmat[,1], pmat[,2], pmat[,3], pmat[,4],
      pmat[,5], pmat[,6], pmat[1, 7], pmat[1, 9], pmat[1, 8])

  } else if (type == "plba3") {
    ## A   B C mean_v mean_w sd_v sd_w  rD   tD   t0 swt
    # arma::mat rplba3(int n, arma::vec A, arma::vec B, arma::vec C,
    #   arma::vec mean_v, arma::vec mean_w, arma::vec sd_v, arma::vec sd_w,
    #   double rD, double tD, double swt, double t0)
    out <- rplba3(n, pmat[,1], pmat[,2], pmat[,3],
      pmat[,4], pmat[,5], pmat[,6], pmat[,7], pmat[1, 8],
      pmat[1, 9], pmat[1, 11], pmat[1, 10])
  } else if (type == "lnr") {
    out <- ggdmc:::rlnrDF(n, pmat[,1],  pmat[,2], pmat[,3], pmat[1,4])
  } else {
    stop("Model type yet created")
  }
  return(out)
}

##' Return ns-npar matrix
##'
##' Contructs a ns x npar matrix, indicating the true paramters
##' used to simualte data. Each row represents a set of parameters for a
##' participant. One should enter either a valid vector or matrix for
##' true parameters (i.e., ps) or a list of (parameter) prior distributions
##' (p.prior). When \code{p.prior} is supplied, true parameters are drawn
##' from prior distributions.
##'
##' @param model a model object
##' @param ns number of subjects.
##' @param p.prior a list of parameter prior distributions
##' @param ps a vector or matirx. Each row indicates a set of true parameters
##' for a participant.
##'
##' @examples
##' model <- ggdmc::BuildModel(
##' p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
##' match.map = list(M=list(s1="r1", s2="r2")),
##' factors   = list(S=c("s1", "s2")),
##' constants = c(st0=0, d=0),
##' responses = c("r1","r2"),
##' type      = "rd")
##'
##' p.prior <- ggdmc::BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "beta", "tnorm", "beta"),
##'   p1=c(a=1, v=0, z=1, sz=1, sv=1, t0=1),
##'   p2=c(a=1, v=2, z=1, sz=1, sv=1, t0=1),
##'   lower=c(0, -5, NA, NA, 0, NA),
##'   upper=c(2,  5, NA, NA, 2, NA))
##'
##' ## Example 1: Randomly generate 2 sets of true parameters from
##' ## parameter priors (p.prior)
##' ggdmc::GetParameterMatrix(model, 2, p.prior)
##' ##            a         v         z        sz       sv        t0
##' ## [1,] 1.963067  1.472940 0.9509158 0.5145047 1.344705 0.0850591
##' ## [2,] 1.512276 -1.995631 0.6981290 0.2626882 1.867853 0.1552828
##'
##' ## Example 2: Use a user-selected true parameters
##' true.vector  <- c(a=1, v=1, z=0.5, sz=0.2, sv=1, t0=.15)
##' ggdmc::GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##' ggdmc::GetParameterMatrix(model, 2, ps = true.vector)
##'
##' ## Example 3: When a user enter arbritary sequence of parameters.
##' ## Note sv is before sz. It should be sz before sv
##' ## See correct sequence, by entering "attr(model, 'p.vector')"
##' ## GetParameterMatrix will rearrange the sequence.
##' true.vector  <- c(a=1, v=1, z=0.5, sv=1, sz = .2, t0=.15)
##' ggdmc::GetParameterMatrix(model, 2, NA, true.vector)
##' ##      a v   z  sz sv   t0
##' ## [1,] 1 1 0.5 0.2  1 0.15
##' ## [2,] 1 1 0.5 0.2  1 0.15
##'
##' @export
GetParameterMatrix <- function(model, ns, p.prior = NA, ps = NA) {
  ## h.simulate.dmc <- function(model, ns=2, n=2, p.prior=NA, ps=NA)
  message1 <- "Parameters are incompatible with model"
  # message2 <- "Must supply either a list of p.prior or a parameter vector."
  # if(anyNA(p.prior) & anyNA(ps)) stop(message2)
  pnames <- names(attr(model, "p.vector"))

  if (anyNA(p.prior)) { ## use ps
    if (is.vector(ps)) {
      if (check_pvec(ps, model)) stop(message1)
      ps    <- ps[pnames]
      pss   <- rep(ps, each = ns)
      psmat <- matrix(pss, ns, dimnames = list(NULL, pnames))
    } else if (is.matrix(ps)) {
      psmat <- matrix(ps, ns, dimnames = list(NULL, pnames))
    } else {
      if ((ns != dim(ps)[1])) stop("ps matrix must have ns rows")
      if (check_pvec(ps[1,], model)) stop(message1)
    }

  } else {  ## use p.prior; random-effect model
    if (!all( pnames %in% names(p.prior))) stop(message1)
    psmat <- ggdmc::rprior(p.prior[pnames], ns)

    ## A nasty way to deal with MG's error gate; ie keep redrawing until we
    ## pass his checks.
    # if (attr(model, "type") == "rd") {
    #   facs <- ggdmc::createfacsdf(model)
    #
    #   for (i in 1:ns) {
    #     for (j in 1:nrow(facs)) {
    #       psmat_allpar <- p.df.dmc(psmat[i,], j, model, FALSE)
    #       psmat_allpar <- checkddm3(psmat_allpar, j, model, p.prior)
    #       psmat[i,] <- as.numeric(psmat_allpar[1, pnames])
    #     }
    #   }
    # }
  }

  rownames(psmat) <- 1:ns
  return(psmat)
}


##' @export
##' @rdname simulate.model
simulate_one <- function(model, n, ps, p.prior = NA) {
  if (check_pvec(ps, model)) stop("p.vector and model incompatible")
  resp <- attr(model, "responses")
  type <- attr(model, "type")
  levs <- attr(model, "factors")
  facs <- ggdmc::createfacsdf(model)
  nvec <- ggdmc::check_n(n, facs)
  dat  <- ggdmc::nadf(nvec, facs, levs)
  row1 <- 1

  ## random is set in the ggdmc_random.R
  for (i in 1:nrow(facs)) {
    pmat <- p.df.dmc(ps, i, model, FALSE) ## simulate use n1.order == FALSE
    rown <- row1 + nvec[i] - 1
    dat[row1:rown, c("RT", "R")] <- random(type, pmat, nvec[i], p.prior)
    row1 <- rown+1
  }

  dat$R <- factor(dat$R, levels = 1:length(resp), labels = resp)
  if (type == "rd") dat <- FlipResponse_rd(model, dat, facs)
  return(dat)
}

##' @rdname simulate.model
##' @export
simulate_many <- function(model, n, ns, p.prior, ps) {

  n  <- ggdmc::GetNsim(model, n, ns)
  ps <- ggdmc::GetParameterMatrix(model, ns, p.prior, ps)

  ismanysub <- ggdmc::ismanymodels(model, ns)
  if(ismanysub) modeli <- model[[1]] else modeli <- model

  ndatai <- cumsum(c(0, matrixStats::rowSums2(n))); ## index boundaries
  datr <- (ndatai[1] + 1):(ndatai[2]); ## Fist subj's trial index
  ## Simulate first subject; modeli should be 'model' class
  dat <- cbind(s = rep.int(1, length(datr)),
    simulate_one(modeli, n[1,], ps[1,], p.prior))

  if (ns > 1) {
    for (i in 2:ns) {
      if (ismanysub) modeli <- model[[i]] else modeli <- model
      datr <- (ndatai[i] + 1):(ndatai[i + 1]) ## continue to index trials
      dat  <- rbind(dat,
        cbind(s = rep.int(i, length(datr)),
          simulate_one(modeli, n[i,], ps[i,], p.prior)))
    }
  }

  dat$s <- factor(dat$s)
  attr(dat, "parameters") <- ps
  ## if ps is not auto-created by p.prior, save the user's p.prior in 'attribute'
  if(!anyNA(p.prior)) attr(dat, "p.prior") <- p.prior
  return(dat)
}



##' Simulate Choice RT Data
##'
##' Simulate stochastic responses either for one subject or multiple subjects.
##' The simulation is based on the \code{model} object. For one subject, the
##' user must supply true parameters, \code{p.vector} at \code{ps} argument.
##' For multiple subjects, the user can supply a matrix (or a row vector),
##' indicating true parameters for each subject, separately on each row
##' (via \code{ps} argument). This is the fixed-effect model. If the user
##' wants to simulate from a random-effect (i.e., hierarchical) model, in which
##' case p.prior must be supplied and ps will be ignored. Note in some cases,
##' a random-effect model may fail to draw data from the model, because
##' true parameters are drawn from \code{p.prior} and a specific model, like
##' DDM, may has certain ranges from different parameters.
##'
##' \code{ps} can be a row vector, in which case each subject has identical
##' parameters. It can also be a matrix with one row per subject, in which
##' case it must have \code{ns} rows. The true values will be saved as
##' "parameters" attribute.
##'
##' @param object a model object
##' @param nsim number of trials/responses. \code{n} can be a single number for a
##' balanced design or matrix for an unbalanced design, where rows are
##' subjects and columns are design cells. If the matrix has one row then all
##' subjects have the same \code{n} in each cell, if it has one column then all
##' cells have the same \code{n}; Otherwise each entry specifies the \code{n}
##' for a particular design subject x design cell combination.
##' @param nsub number of subjects
##' @param p.prior parameter priors. A list of distributions based on which
##' the true parameters fro each subject are drawn.  It is usually created by
##' \code{BuildPrior} and will be saved as "p.prior" attribute.
##' @param ps p.vector matrix. Each row represent a subject.
##' @return a data frame
##' @examples
##' model <- ggdmc::BuildModel(
##'   p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1",
##'   sv = "1", t0 = "1", st0 = "1"),
##'   match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'   factors   = list(S = c("s1", "s2")),
##'   constants = c(st0 = 0, d = 0),
##'   responses = c("r1", "r2"),
##'   type      = "rd")
##' @export
simulate.model <- function(object, nsim = NA, nsub = NA, p.prior = NA, ps = NA) {
  if (is.na(nsub)) {
    if (is.na(nsim)) stop("How many response you want to generate? Must supply n")
    if (anyNA(ps)) stop("Some true parameters missing")
    out <- simulate_one(object, nsim, ps)
  } else {
    message1 <- "Must supply either p.prior or p.vector."
    if (anyNA(p.prior) & anyNA(ps)) stop(message1)
    out <- simulate_many(object, nsim, nsub, p.prior, ps)
  }
  return(out)
}


