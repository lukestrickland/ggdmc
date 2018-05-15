#' Create an array of all factorial combinations of factor levels
#'
#' Take a list storing one or multiple string vectors and glues the strings in
#' these vectors with a dot symbol.
#'
#' @param factors a list storing one or multiple factor vector. Each vector
#' can have different numbers of level.  For example,
#' \code{f <- list(E = c("nc", "wc"), S = c("n", "w", "pn", "pw"))}
#' @return a table showing the combinations of factor levels.
#'
#' @examples
#' ## Example 1
#' factors <- list(S = c("s1", "s2"))
#' ggdmc:::make.level.array(factors)
#' ##   S
#' ##  s1   s2
#' ## "s1" "s2"
#'
#' factors <- list(S = c("left", "right"))
#' ggdmc:::make.level.array(factors)
#' ##      S
#' ##   left   right
#' ##  "left" "right"
#'
#' ## Example 2
#' factors <- list(A = c("00", "45", "90", "135", "180"),
#'                 S = c("mirror", "normal"))
#' ggdmc:::make.level.array(factors)
#' ##       S
#' ## -------------------------------
#' ##   A        mirror       normal
#' ## -------------------------------
#' ##  00   "00.mirror"  "00.normal"
#' ##  45   "45.mirror"  "45.normal"
#' ##  90   "90.mirror"  "90.normal"
#' ## 135  "135.mirror" "135.normal"
#' ## 180  "180.mirror" "180.normal"
#'
#' ## Example 3
#' factors <- list(E = c("nc", "wc"),
#'                 S = c("n", "w", "pn", "pw"))
#' ggdmc:::make.level.array(factors)
#'
#' ##         S
#' ## ---------------------------------
#' ## E       n      w      pn      pw
#' ## ---------------------------------
#' ## nc "nc.n" "nc.w" "nc.pn" "nc.pw"
#' ## wc "wc.n" "wc.w" "wc.pn" "wc.pw"
#'
#' @export
make.level.array <- function(factors = NA) {
  if (all(is.na(factors))) stop("Found no factors!")
  out <- factors[[1]]
  nf  <- length(factors)

  if ( nf > 1 ) {
    for (i in 2:nf) out <- outer(out, factors[[i]], "paste", sep = ".")
    dimnames(out) <- factors
  } else {
    out <- array(out, dim = length(factors[[1]]), dimnames = factors)
  }
  return(out)
}

# p.map = list(A="1",B=c("R"),t0="1",
#   mean_v=c("A","S","M"),sd_v="M",st0="1")
#   match.map = list(M=list(mirror="MIRROR",normal="NORMAL"))
#   factors=list(A=c("00","45","90","135","180"),S=c("mirror","normal"))
#   constants = c(sd_v.false=1,st0=0,
#     mean_v.00.mirror.false=0,mean_v.45.mirror.false=0,
#     mean_v.90.mirror.false=0,mean_v.135.mirror.false=0,
#     mean_v.180.mirror.false=0,mean_v.00.normal.false=0,
#     mean_v.45.normal.false=0,mean_v.90.normal.false=0,
#     mean_v.135.normal.false=0,mean_v.180.normal.false=0)
#   responses = c("NORMAL","MIRROR")
#   type="norm"

# p.map=list(t0="1",A="1",B=c("E","R"),mean_v=c("E","V"),sd_v="1",st0="1")
# factors=list(E=c("nc", "wc"),S=c("n","w","pn","pw"))
# responses=c("N","W","P")
# match.map=list(M=list(n="N",w="W",pw="P",pn="P"),V=map3)
# constants=c(sd_v=1,st0=0)
# type="norm"; posdrift=TRUE


######### Modeling checks
checknull    <- function(p.map, factors, responses) {
  if (is.null(p.map)) stop("Must supply p.map")
  if (is.null(factors)) stop("Must supply factors")
  if (is.null(responses)) stop("Must supply responses")
}

checkmap <- function(match.map, responses) {

  ## Extract essential information
  nmap      <- length(match.map)
  mapnames  <- names(match.map)
  message1 <- "match.map must be a list of lists"
  message2 <- "match.map has no list named M"
  message3 <- "match.map$M has index or name not in response names. If you use
  integers, it must be continuous."
  message4 <- "Not all response names are scored by match.map$M"
  message5 <- "Entries in match.map besides M must be factors"

  if (nmap < 1 || class(match.map[[1]]) != "list") stop(message1)
  if (!any(names(match.map) %in% "M") ) stop(message2)
  map.names  <- mapnames[mapnames != "M"]
  umap       <- unlist(match.map$M) ## unlisted match map

  ## We allow the user to enter integers (1, 2, ..., n) to represent
  ## "response 1", "response 2", ... "response n". If "is.numeric" evaluates to
  ## TRUE, it means the user enters integers as response types. Hence, we use
  ## lapply to pick sequentially 1st stimulus type matching to 1st response
  ## type etc. Response types are stored in another variable "responses"
  if (is.numeric(umap)) {
    match.map$M <- lapply(match.map$M, function(x){responses[x]})
    umap <- unlist(match.map$M)
  }
  if ( !all(umap %in% responses) ) stop(message3)
  if ( !(all(sort(responses) == sort(unique(umap)))) ) stop(message4)
  if ( nmap > 1 && !all(lapply(match.map[mapnames != "M"], class)=="factor") )
    stop(message5)
  if ( length(match.map)>1 &&
      !all(lapply(match.map[names(match.map)!="M"], class)=="factor") )
    stop("Entries in match.map besides M must be factors")

}
checkddm1 <- function(p.map, responses,  type) {
  nr       <- length(responses)
  message1 <- "DDM only applicable for two responses"
  message2 <- "Cannot use M or R in DDM p.map"
  if (type == "rd" & (nr != 2)) stop(message1)

  # Does the par use a match factor or a response factor?
  is.M <- unlist(lapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "M"))
  }))
  is.R <- unlist(lapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "R"))
  }))
  if (type =="rd"  & ( any(is.M) | any(is.R) )) stop(message2)
}

checkddm3 <- function(pmat, i, model, p.prior, debug = FALSE) {
  ## RT = .3, precision = 2.5 are dummies
  if (anyNA(p.prior)) stop("p.prior not found in checkddm3")
  ps <- c(pmat$a[1], pmat$v[1], pmat$z[1], pmat$d[1], pmat$sz[1], pmat$sv[1],
    pmat$t0[1], pmat$st0[1], .3, 2.5)
  isbadparameter <- ggdmc::checkddm2(ps)
  j <- 0

  repeat {
    if (isbadparameter) {
      j <- j+1
      ps <- rprior(p.prior)
      pmat <- p.df.dmc(ps, i, model, FALSE)
      ps <- c(pmat$a[1], pmat$v[1], pmat$z[1], pmat$d[1], pmat$sz[1], pmat$sv[1],
        pmat$t0[1], pmat$st0[1], .3, 2.5)
      isbadparameter <- checkddm2(ps)
      if (debug) message("ps changed")
    } else if (j >= 1e3) {
      stop("Fail to simulate DDM data");
      break
    } else {
      break
    }

  }
  return(pmat)
}


checkfacresp <- function(p.map, match.map, factors, responses) {


  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]
  map.levels <- unlist(map.names, levels)

  ufac <- unlist(factors)    ## unlisted factors
  uml  <- unlist(map.levels) ## unlisted map.levels
  reservednames <- c("1", "R", "R2", "M", "s")
  tfnames       <- c("true", "false")
  ufacmaplv     <- unlist(c(factors, map.levels))

  message1 <- "Do not use s, M, R, R2 or 1 as a factor name"
  message2 <- paste(map.names,"used in match.map, can not use as a factor name")
  message3 <- "All factors levels must be unqiue"
  message4 <- "All match.map levels must be unqiue"
  message5 <- "\"true\" and \"false\" cannot be used as factor levels"
  message6 <- "\"true\" and \"false\" cannot be used as match.map levels"
  message7 <- "Factor levels cannot overlap match.map levels"

  # Check factors and add responses
  if ( any(names(factors) %in% reservednames) ) { stop(message1) }
  if ( any(names(factors) %in% mapnames) ) {      stop(message2) }
  if ( length(ufac) != length(unique(ufac)) ) {   stop(message3) }
  if ( length(uml) != length(unique(uml)) ) {     stop(message4) }
  if ( any(ufac %in% tfnames) ) {                 stop(message5) }
  if ( any(map.levels %in% tfnames) ) {           stop(message6) }
  if ( length(ufacmaplv) != length(unique(ufacmaplv)) ) { stop(message7) }
}
checkpmap    <- function(p.map) {
  pmapnames <- names(p.map)
  has.dot   <- unlist(lapply(strsplit(pmapnames, "[.]"), length)) > 1
  message1 <- paste("Dots not allowed in p.map names, fix:", pmapnames[has.dot])
  if (any(has.dot)) stop(message1)

  # Check M and R are last
  if (any(unlist(lapply(p.map, function(x){any(x == "M") && x[length(x)]!="M"}))))
    stop("M factors must always be last")
  if (any(unlist(lapply(p.map, function(x){any(x == "R") && x[length(x)]!="R"}))))
    stop("R factors must always be last")
}
addmap2facs  <- function(match.map, factors, responses) {
  umap <- unlist(match.map$M)
  if (is.numeric(umap)) match.map$M <- lapply(match.map$M, function(x){responses[x]})
  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]

  factors$R     <- responses
  if (!is.null(match.map)) {
    factors$M <- c("true", "false")
    for (i in map.names) factors[[i]] <- levels(match.map[[i]])
  }
  return(factors)
}
getparnames  <- function(p.map, factors) {
  ## Make parameter names
  out <- character(0)
  for ( i in names(p.map) )
  {
    if ( length(p.map[[i]]) == 1 && p.map[[i]] == "1" ) {
      new.names <- i
    } else {
      new.names <- paste(i, factors[[p.map[[i]][1]]], sep=".")
      if ( length(p.map[[i]]) > 1 ) {
        for (j in 2:length(p.map[[i]])) {
          new.names <- as.vector(outer(
            new.names, factors[[p.map[[i]][j]]], "paste", sep="."))
        }
      }
    }
    out <- c(out, new.names)
  }

  return(out)
}

checkdesign  <- function(match.map, levelarray) {
  ## Check match.map and expand
  for (i in names(match.map$M)) {
    match.position <- grep(i, levelarray)
    if ( length(match.position)==0 )
      stop(paste(i, "in match.map is not in the design"))
  }
}

splitfacresp <- function(colfac, responses, maplevels) {
  res <- lapply(colfac, function(x){
    if ( length(x) == 0 ) out <- c(NA, NA)
    if ( length(x) == 1 ) {
      if ( x %in% c(responses, "true", "false", maplevels) )
        out <- c(NA, x) else out <- c(x, NA)
    }
    if ( length(x) > 1 )
      if ( x[length(x)] %in% c(responses, "true", "false", maplevels) )
        out <- c(paste(x[-length(x)], collapse="."), x[length(x)]) else
          out <- paste(x, collapse=".")
        out
  })
  return(res)
}

## Check reserved factor name (M, R, and other user indicated names in match.map)
checkrn <- function(p.map, match.map) {
  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]
  ## Reserved uppercase letters we don't want the user to use
  is.M <- unlist(lapply(p.map, function(x) {
    any(unlist(strsplit(x, "[.]") == "M"))
  }))

  ## Does the par use a response factor
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x, "[.]") == "R"))
  }))

  ## Does the par use a map factor, other than M
  is.map <- unlist(lapply(p.map,function(x) {
    any(unlist(strsplit(x, "[.]") %in% map.names))
  }))

  out <- cbind(is.M, is.R, is.map)

  if ( any(apply(out, 1, function(x){ sum(x) > 1 })) ) {
    message("M and R are reserved names. Also if you have used a name in")
    message("match.map other than M. It should not be used to name a paramter.")
    message("I will stop here. Please check your match.map or")
    stop("consult the package maintainer.")
  }
  # if ( any(is.map) ) {
  #   p.map.name <- lapply(p.map, function(x){
  #     unlist(strsplit(x,"[.]" ))[
  #       unlist(strsplit(x,"[.]")) %in% map.names]
  #   })
  #   nr <- length(responses)
  #   n  <- length(level.array)
  #   map.shuffle <- matrix(aperm(array(1:n, dim = c(n/nr, nr, nr)), c(1,3,2)), ncol=nr)
  # }

  if (any(is.map)) { stop("map.shuffle is pending.") }
  return(out)
}

makemodelarray <- function(match.map, responses, col.par, rnamemat, allpar,
  level.array) {

  umap <- unlist(match.map$M)
  if (is.numeric(umap)) match.map$M <- lapply(match.map$M, function(x){responses[x]})
  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]
  map.levels <- unlist(map.names, levels)


  modeldn  <- list(as.vector(level.array), allpar, responses)
  modeldim <- c(length(level.array), length(allpar), length(responses))
  use.par  <- array(NA, modeldim, modeldn)
  is.M   <- rnamemat[,1]
  is.R   <- rnamemat[,2]
  is.map <- rnamemat[,3]

  ## col.par = column parameter type (1st name)
  col.par <- strsplit(modeldn[[2]], "[.]")
  col.fac <- lapply(col.par, function(x){x[-1]})
  col.par <- unlist(lapply(col.par, function(x){x[1]}))
  ## split into fac and resp
  col.fac <- splitfacresp(col.fac, responses, map.levels)
  col.resp <- unlist(lapply(col.fac, function(x){x[2]}))
  col.fac  <- unlist(lapply(col.fac, function(x){x[1]}))

  row.fac <- strsplit(modeldn[[1]], "[.]")
  row.fac <- unlist(lapply(row.fac, function(x){
    paste(x[-length(x)], collapse=".")}))


  for ( p in unique(col.par) ) {
    is.col <- p==col.par
    ncols  <- sum(is.col)
    if ( ncols == 1 ) {
      use.par[, is.col,] <- TRUE
    } else {
      # there are parameter subtypes
      for ( i in 1:ncols ) {
        ## each parameter subtype
        ## select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE, dim(use.par)[1])
        if ( !is.na(tmp) ) is.fac.row[!grepl(tmp, row.fac)] <- FALSE
        ## set not applicable rows to false
        use.par[!is.fac.row, is.col, ][,i,] <- FALSE
        if ( is.M[p] ) {
          # has a match factor
          for ( j in names(match.map$M) ) {
            # response cell
            correct.response <- match.map$M[[j]]
            is.rcell <- is.fac.row & grepl(j, row.fac)
            for ( k in responses ) {
              # responses
              if ( k==correct.response ) {
                if ( grepl("true", col.resp[is.col][i]) )
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE
              } else {
                if ( grepl("false", col.resp[is.col][i]) )
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE
              }
            }
          }
        } else if ( is.R[p] ) {
          for ( k in responses )
            use.par[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
        }  else if ( is.map[p] ) {
          use.par[is.fac.row,is.col,][,i,] <-
            match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
        } else use.par[is.fac.row,is.col,][,i,] <- TRUE
      }
    }
  }

  if ( any(is.na(use.par)) ) stop("Some cells of the map were not assigned!")
  return(use.par)
}

getcolpar <- function(names.par) {
  col.par <- strsplit(names.par, "[.]")   ## dim2
  col.fac <- lapply(col.par, function(x){x[-1]})
  col.par <- unlist(lapply(col.par, function(x){x[1]}))
  return(col.par)
}

checkconst <- function(use.par, constants) {
  all.par <- use.par[1,,1]  ## First array, first row
  all.par[1:length(all.par)] <- NA
  if (length(constants)>0) {
    if ( !any(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }
  return(all.par)
}

addisr1 <- function(match.map, responses, use.par) {
  is.r1 <- rep(FALSE, length(row.names(use.par)))
  names(is.r1) <- row.names(use.par)
  is.r1[unlist(lapply(
    lapply(
      as.list(names(match.map$M)[match.map$M==responses[1]]),
      function(x)grepl(x, row.names(use.par))
    ),
    function(x)which(x==TRUE)))] <- TRUE
  attr(use.par, "is.r1") <- is.r1
  return(use.par)
}

getn1order <- function(responses, level.array, use.par) {
  resp <- unlist(lapply(strsplit(level.array,"[.]"),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp), ncol=nr)
  for (i in 1:length(resp))
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses], c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)
  return(n1.order)

}

getmatchcell <- function(match.map, n1.order) {
  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell))
    match.cell[match.num[match.num %in%
        grep(names(match.map$M[i]),names(match.cell))]] <- TRUE
  }
  return(match.cell)
}

##' Is a p.vector compatible with model?
##'
##' Check if the user supplies a p.vector, incompatiable with the p.vector
##' created by \code{model.dmc}.
##'
##' @param pvec p.vector. I use pvec to avoid naming confound with S3 method.
##' @param model a model object
##' @export
check_pvec <- function(pvec, model) {
  modpvec <- names(attr(model, "p.vector"))
  ism1 <- modpvec %in% names(pvec)
  ism2 <- names(pvec) %in% modpvec
  bad  <- any(!ism1)
  if (bad) warning(paste("Parameter", modpvec[!ism1],
    "in model not present in p.vector\n"))
  bad <- bad | any(!ism2)
  if (any(!ism2)) warning(paste("Parameter",
    names(pvec)[!ism2], "in p.vector not present in model\n"))
  invisible(bad)
}

check_cell <- function(cell, model) {

  dim1 <- dimnames(model)[[1]]
  if(is.numeric(cell)) {
    if ( cell > length(dim1) ) stop("cell out-of-range!")
    cell <- dim1[cell]
  }
  cell
}

check_rd <- function(type, model) {
  ## depreciate this
  if(type != "rd") { isr1 <- 0 } else { isr1 <- attr(model, "is.r1") }
  isr1
}

##' Get Response-Parameter Data Frame
##'
##' This is a wrapper function for C++ function, \code{pdfdmc}.
##' \code{p.df.dmc} extracts and arranges the values in a parameter vector,
##' which either is from a sampler or directly provided by the user, to a
##' response x parameter matrix. The matrix is usually used by a likelihood
##' function that picks a trial corresponding to its an experimental design cell
##' and calculates its probability density.
##'
##' @param p.vector a parameter vector
##' @param cell a string or a integer indicating a design cell, e.g.,
##' \code{s1.f1.r1} or 1. Note the integer cannot exceed the number of cell.
##' use \code{length(dimnames(model))} to check the upper bound.
##' @param model a model, a cognitive model usually created by \code{model.dmc}
##' @param n1order a boolean switch to use node 1 ordering. This is only for
##' LBA model and its n1PDF likelihood function.
##' @return each row corresponding to the model parameter for a response.
##' When \code{n1.order} is FALSE, p.df.dmc returns a martix in natural order,
##' which is used by \code{simulate.dmc}. By default \code{n1.order} is TRUE,
##' the returned matrix is usually used by n1PDF and \code{likelihood_norm}.
##' @export
##' @examples
##' model <- ggdmc:::model.dmc(
##'   p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2")),
##'   constants = c(st0 = 0, sd_v = 1),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1, mean_v.false = .25)
##' ggdmc:::print.cell.p(p.vector, model)
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
##'
##' for(i in dimnames(model)[[1]]) {
##'   tmp <- ggdmc::p_df(p.vector, i, pnames, allpar, parnames, model, type,
##'   dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
##'   isr1, n1idx, TRUE)
##'   cat(i, "\n")
##'   print(tmp)
##' }
##'
##' ## Older note
##' m1 <- model.dmc(
##'   p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="F",t0="1",st0="1"),
##'   match.map = list(M=list(s1="r1", s2="r2")),
##'   factors   = list(S=c("s1","s2"), F=c("f1","f2")),
##'   constants = c(st0=0,d=0),
##'   responses = c("r1","r2"),
##'   type      = "rd")
##'
##' m2 <- model.dmc(
##'   p.map = list(A = "1", B = "1", mean_v = "M", sd_v = "1",
##'   t0 = "1", st0 = "1"),
##'   constants = c(st0 = 0, sd_v = 1),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2")),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' pvec1 <- c(a=1.15, v.f1 = -0.10, v.f2 = 3, z=0.74, sz=1.23,
##'   sv.f1=0.11, sv.f2=0.21, t0=0.87)
##' pvec2 <- c(A = .75, B = .25, mean_v.true = 2.5, mean_v.false = 1.5,
##'   t0 = .2)
##'
##' print.cell.p(pvec1, m1)
##' print.cell.p(pvec2, m2)
##'
##' accMat1 <- ggdmc::p.df.dmc(pvec1, "s1.f1.r1", m1, FALSE)
##' accMat2 <- ggdmc::p.df.dmc(pvec2, "s1.r1",    m2, FALSE)
##'
##' ##    a    v   t0    z d   sz   sv st0
##' ## 1.15 -0.1 0.87 0.26 0 1.23 0.11   0
##' ## 1.15 -0.1 0.87 0.26 0 1.23 0.11   0
##'
##' ##    A b  t0 mean_v sd_v st0
##' ## 0.75 1 0.2    2.5    1   0
##' ## 0.75 1 0.2    1.5    1   0
##'
##' \dontrun{
##' library(microbenchmark)
##' res <- microbenchmark(ggdmc::p.df.dmc(pvec2, "s1.r2", m2, FALSE),
##'  p.df.dmc(pvec2, "s1.r2", m2, FALSE),
##'  ggdmc::p_df(pvec2, "s1.r2", pnames, allpar, parnames, m2, type,
##'    dim1, dim2, dim3, isr1, n1, FALSE), times = 1e2)
##' res
##' }
##'
##' ## Unit: microseconds
##' ##                        expr
##' ## ggdmc::p.df.dmc(pvec2, ...)
##' ## p.df.dmc(pvec2, ... )
##' ## ggdmc::p_df(pvec2, ... )
##' ##     min       lq      mean   median       uq     max neval cld
##' ## 119.359 125.5400 132.10528 131.1270 138.2160 153.581   100   b
##' ## 164.826 171.7745 178.52793 176.3145 180.4355 349.346   100   c
##' ##  14.389  17.4265  24.17361  25.5630  28.7055  64.953   100   a
##'
##' ## Test internal function, p_df()
##' ## cell <- 1
##' ## n1order <- TRUE
##' ## pnames   <- names(attr(model, "p.vector"))
##' ## allpar   <- attr(model, "all.par")
##' ## parnames <- attr(model, "par.names")
##' ## type     <- attr(model, "type")
##' ## n1       <- attr(model, "n1.order")
##' ## resp     <- attr(model, "responses")
##' ## dim1 <- dimnames(model)[[1]]
##' ## dim2 <- dimnames(model)[[2]]
##' ## dim3 <- dimnames(model)[[3]]
##' ## if(type == "rd") { isr1 <- attr(model, "is.r1") } else { isr1 <- 0 }
##' ## if(is.numeric(cell)) {
##' ##   if ( cell > length(dim1) ) stop("cell out-of-range!")
##' ##   cell <- dim1[cell]
##' ## }
##'
p.df.dmc <- function(p.vector, cell,  model, n1order = TRUE)
{
  pnames   <- names(attr(model, "p.vector"))
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1       <- attr(model, "n1.order")
  resp     <- attr(model, "responses")
  cell     <- check_cell(cell, model)
  isr1     <- check_rd(type, model)

  out <- as.data.frame(p_df(p.vector, cell, pnames, allpar, parnames,
                       model, type, dimnames(model)[[1]], dimnames(model)[[2]],
                       dimnames(model)[[3]], isr1, n1, n1order))

  if(type == "rd") {
    names(out) <- c("a","v","z","d","sz","sv","t0","st0")
    rownames(out) <- attr(model, "response")
  }

  if(type %in% c("norm", "norm_pda", "norm_pda_gpu")) {
    if (dim(out)[[2]] != 6) {
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0", "nacc")
    } else {
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0")
    }
  }

  if(type %in% c("plba0_gpu")) {
    names(out) <- c("A", "b", "mean_v", "sd_v", "mean_w", "rD", "t0","swt")
  }

  if(type %in% c("plba1", "plba1_gpu")) {
    names(out) <- c("A", "b", "mean_v", "sd_v", "mean_w", "rD", "t0","swt")
  }

  if(type %in% c("plba2")) {
    names(out) <- c("A","b", "mean_v","mean_w","sd_v", "sd_w", "rD","t0","swt")
  }

  if(type %in% c("plba3")) {
    names(out) <- c("A", "B", "C","mean_v","mean_w","sd_v","sd_w", "rD", "tD",
      "t0","swt")
  }

  if(type %in% c("lnr")) {
    names(out) <- c("meanlog", "sdlog", "t0", "st0")
    rownames(out) <- attr(model, "response")
  }

  return(out)
}


##' @export
##' @rdname simulate.model
createfacsdf <- function(model) {
  dim1 <- dimnames(model)[[1]]
  responses <- attr(model, "responses")
  levs      <- attr(model, "factors")
  nr <-  length(responses)

  facs  <- lapply(strsplit(dim1, "[.]"), function(x){x[-length(x)]})
  nf    <- length(facs)
  facs  <- facs[1:( nf/nr )]
  fnams <- names(levs)
  facs  <- data.frame(t(matrix(unlist(facs), length(fnams))))
  names(facs) <- fnams
  return(facs)
}

##' @importFrom data.table is.data.table
##' @export
##' @rdname simulate.model
check_n <- function(n, facs) {
  if ( length(n) == 1 ) n <- rep(n, dim(facs)[1])
  test_n <- ifelse(is.data.frame(n), (nrow(n) != dim(facs)[1]), (length(n) != dim(facs)[1]))
  if (test_n) stop(paste("n must either be length 1 or", dim(facs)[1], "for this model."))

  if (data.table::is.data.table(n)) {
    n <- n$N
  } else {
    if ( !is.null(dimnames(n)) ) n <- merge(facs, data.frame(n), sort = F)$Freq
    for (i in 1:dim(facs)[1]) {if (n[i] < 0) warning("n is less than 0")}
  }

  return(n)
}


##' @export
##' @rdname simulate.model
nadf <- function(n, facs, levs) {
  # create a data-frame container
  data <- data.frame(lapply(facs, rep, n))

  for (i in names(levs)) data[[i]] <- factor(as.character(data[[i]]), levels=levs[[i]])
  data$R <- NA
  data$RT <- NA
  return(data)
}

##' @export
##' @rdname simulate.model
FlipResponse_rd <- function(model, data, facs) {
  cell.names <- apply(data[, 1:length(names(facs)), drop=F], 1, paste, collapse=".")
  M <- attr(model, "match.map")$M
  R <- attr(model, "responses")
  for ( i in names(M) )
    if ( M[[i]] == R[1] )
      data[grep(i, cell.names),"R"] <- as.character(
        factor(as.character(data[grep(i, cell.names), "R"]),
          levels=R, labels=R[2:1]))
  return(data)
}

#' @export
#' @rdname BindDataModel
check_dmi <- function(data, model) {
  ## Remember to make CheckDMI and check_dmi consistent
  message1 <- "Model list is for multiple subjects\nNo s column is found in
 data model instance"
  message2 <- "List of models must be same length as number of subjects"
  message3 <- "data must be a data frame"

  if (is.list(model)) {
    if (!any(names(data)=="s")) stop(message1)
    if (length(model) != length(levels(data$s))) stop(message2)
    subject.models <- TRUE
    modeli <- model[[1]]
  } else {
    subject.models <- FALSE
    modeli <- model
  }

  fnams   <- names(attr(modeli, "factors"))
  factors <- attr(modeli, "factors")
  resps   <- attr(modeli, "responses")
  message4 <- paste("data must have columns named:", paste(fnams, collapse=" "),
    "R", "RT")

  if ( !is.data.frame(data) ) stop(message3)

  # check data
  if ( !all(c(fnams, "R", "RT") %in% names(data)) ) stop(message4)
  for ( i in fnams ) {
    factori <- factors[[i]]
    if ( !all(sort(factori) == sort(levels(data[,i]))) )
      stop(paste("Factor", i, "must have levels:", paste(factori, collapse=" ")))
  }

  if ( !all(sort(resps) == sort(levels(data[, "R"]))) )
    stop(paste("data$R must have levels:", paste(resps, collapse=" ")))
  if ( !is.numeric(data$RT) ) stop("data$RT must be of type numeric")

  list(issm=subject.models, model=modeli)
}

#' Bind Data and Models
#'
#' Binding a data frame and a DMC model description. This function also check
#' if they are compatible and adding a \code{cell.index} and many other
#' attributes to the data frame in order to speed likelihood computation.
#'
#' @param data a data frame stored choice-RT data
#' @param model a DMC model
#' @param npda number of model simulations
#' @param bw bandwidth
#' @param gpuid indicate using which GPU card. Default is 0.
#' @param debug if print debugging information
#' @examples
#' model <- BuildModel(
#'    p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'    constants = c(st0=0,d=0),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors   = list(S=c("s1","s2")),
#'    responses = c("r1","r2"),
#'    type      = "rd")
#'
#' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat <- simulate(model, p.vector, 128)
#' dmi <- BindDataModel(model, dat)
#'
#' @export
BindDataModel <- function(data, model, n.pda=16384, bw = 0.01, gpuid=0, debug=FALSE)
{
  res <- check_dmi(data, model)
  subject.models <- res$issm
  modeli <- res$model
  fnams <- names(attr(modeli, "factors"))

  if ( any(names(data)=="s") ) # more than one subject
  {
    dat <- data
    data <- vector(mode="list", length=length(levels(dat$s)))
    names(data) <- levels(dat$s)
    if (subject.models) names(model) <- levels(dat$s)
    is.sim <- !is.null(attr(dat, "parameters"))
    for (s in names(data))
    {
      if (subject.models) modeli <- model[[s]] else modeli <- model
      data[[s]] <- dat[dat$s==s, names(dat)!="s"]
      # add model and index attribute to data
      cells <- apply(data[[s]][, c(fnams,"R")], 1, paste, collapse=".")
      cell.index <- vector(mode="list", length=dim(modeli)[1])
      names(cell.index) <- row.names(modeli)
      for ( i in names(cell.index) ) cell.index[[i]] <- cells %in% i
      attr(data[[s]], "cell.index") <- cell.index
      attr(data[[s]], "cell.empty") <-
        unlist(lapply(cell.index,function(x){sum(x)}))==0
      attr(data[[s]], "model") <- modeli
      if (is.sim) attr(data[[s]], "parameters") <- attr(dat, "parameters")[s,]

      attr(data[[s]], "n.pda") <- n.pda
      attr(data[[s]], "bw")    <- bw
      attr(data[[s]], "gpuid") <- gpuid
      attr(data[[s]], "debug") <- debug
    }

  } else { # one subject
    # add model and index attribute to data
    cells      <- apply(data[,c(fnams, "R")], 1, paste, collapse=".")
    cell.index <- vector(mode="list", length=dim(model)[1])
    names(cell.index) <- row.names(model)
    for ( i in names(cell.index) )
      cell.index[[i]] <- cells %in% i
    attr(data, "cell.index") <- cell.index
    attr(data, "cell.empty") <-
      unlist(lapply(cell.index, function(x){sum(x)}))==0
    attr(data, "model") <- model
    attr(data, "n.pda") <- n.pda
    attr(data, "bw")    <- bw
    attr(data, "gpuid") <- gpuid
    attr(data, "debug") <- debug
  }
  data
}



######### CUSTOM MAP MAKING

#' @rdname BuildModel
#' @export
MakeEmptyMap <- function(FR, levels)
{
  ## This derives from DMC's empty.map
  level.array <- make.level.array(FR)
  map <- rep("", length = length(level.array))
  names(map) <- level.array
  factor(map, levels = levels)
}


#' @rdname BuildModel
#' @export
assignMap <- function(map, value = "", eq.list = list(),
                       funs=NULL, include=NA, match.values=NA)
{

  if (any(is.na(include))) ok <- rep(TRUE,length(map)) else
  {
    ok <- grepl(include[1],names(map))
    if (length(include)>1) for (i in 2:length(include))
      ok <- ok | grepl(include[i],names(map))
  }
  if ( length(eq.list)==0 ) # Fill in empty elements if included
    map[is.na(map) & ok] <- value else if (all(unlist(lapply(eq.list,length))==1))
  {
    if (all(is.na(match.values)) || length(match.values)!=length(eq.list))
      stop("If eq.list only has length one entries match.value must contain target for each entry")
    match.mat <- matrix(unlist(lapply(eq.list,function(x){
      (unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
        y[x]})))})),ncol=length(eq.list))
    map[apply(match.mat,1,function(x){all(x==match.values)})] <- value
  } else {
    if ( is.null(funs) ) funs <- lapply(eq.list,function(x){
      list("identity","identity")
    })
    if (length(funs)!=length(eq.list))
      stop("Must specify one function pair in funs for each entry in eq.list")
    if ( class(funs)!="list" || !all(unlist(lapply(funs,class))=="list") )
      stop("funs must be  list of lists")
    if ( !all(unlist(lapply(funs,length))==2) )
      stop("Each entry in funs must be length 2")
    funs <- lapply(funs,function(x){lapply(x,function(y){
      if ( is.character(y) && y=="" ) "identity" else y
    })})
    pair.list <- lapply(eq.list,function(x){
      matrix(unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
        y[x]})),nrow=2)})
    map.mat <- matrix(nrow=length(map),ncol=length(eq.list))
    for ( i in 1:length(eq.list) )
      map.mat[,i] <- apply(pair.list[[i]],2,function(x){
        do.call(funs[[i]][[1]],list(x[1])) ==
        do.call(funs[[i]][[2]],list(x[2]))
      })
      map[apply(map.mat,1,all) & ok] <- value
  }
  map
}

# map <- empty.map()
# map <- assign.map( map,value="true",eq.list=list(c(1,2)) )
# map <- assign.map(map,value="false")
