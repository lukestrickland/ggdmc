

#' @export
score <- function(dmi, digit = 2) {
  correct <- tolower(dmi$S)==tolower(dmi$R)
  cat("Accuracy: \n"); print(round(mean(correct), digit))
  mrt <- tapply(dmi$RT, list(correct), mean)
  iqr <- tapply(dmi$RT, list(correct), IQR)
  SD  <- tapply(dmi$RT, list(correct), sd)
  out <- data.frame(rbind(mrt, iqr, SD))
  names(out) <- c("error", "correct")
  print(round(out, digit))
  invisible(out)
}

#' Creat a Model Object
#'
#' \code{model.dmc} Creates an array object ("model") with a set of attributes
#' specifying a particular model and parameterisation. Call \pkg{coda} to
#' summarise the model parameters in a DMC samples with multiple participants
#' at the hyper level.
#'
#' @param p.map list factors and constants for parameters
#' @param factors Factor names and levels
#' @param responses Response (accumulator) names
#' @param match.map Scores responses
#' @param constants Parameters set to constant value
#' @param type model type. Should go to class in the future
#' @param posdrift only used by norm (ie LBA model)
#' @param verbose Print p.vector, constants and type
#' @export
#' @examples
#' m1 <- BuildModel(
#'    p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors=list(S=c("s1","s2"), F=c("f1","f2")),
#'    constants = c(st0=0,d=0),
#'    responses = c("r1","r2"),
#'    type = "rd")
BuildModel <- function(
  p.map,                        # list factors and constants for parameters
  responses,                    # Response (accumulator) names
  factors = list(dummy = "1"),  # Factor names and levels
  cvs = NULL,                   # Names of trial covariates (in data)
  responses2 = NULL,            # Second response name (multi-threshold models)
  match.map = NULL,             # Scores responses
  constants = numeric(0),       # Parameters set to constant value
  constant.prior = NULL,        # Parameter sampled from a fixed prior
  type = "norm",                # model type
  posdrift = TRUE,              # only used by norm
  verbose = TRUE,               # Print p.vector, constants and type
  check.model = TRUE)
{

  grepl.dot <- function(pattern,x) {
    # Splits up pattern at ".", to get n parts. Matches each part to x and
    # returens TRUE for each x where exactly n matches in any order.
    grepl.exact <- function(xdot,pattern) {
      xvec <- strsplit(xdot, ".", fixed=TRUE)[[1]]
      any(pattern==xvec)
    }

    ps <- strsplit(pattern,".",fixed=TRUE)[[1]]
    out <- sapply(x,grepl.exact,pattern=ps[1])
    if (length(ps)>1) {
      for (i in ps[-1])
        out <- rbind(out,sapply(x,grepl.exact,pattern=i))
      apply(out,2,sum)==length(ps)
    } else out
  }

  # Check requried inputs supplied
  if (is.null(p.map)) stop("Must supply p.map")
  if (is.null(factors)) stop("Must supply factors")
  if (is.null(responses)) stop("Must supply responses")

  # Check factors
  if ( length(unlist(factors)) != length(unique(unlist(factors))) )
    stop("All factors levels must be unqiue")
  if ( any(names(factors) %in% c("1","R","R2","s")) )
    stop("Do not use s, R, R2, or 1 as a factor name")
  # Check no parameter names have a dot
  has.dot <- unlist(lapply(strsplit(names(p.map),".",fixed=TRUE),length))>1
  if ( any(has.dot) )
    stop(paste("Dots not allowed in p.map names, fix:",paste(names(p.map)[has.dot])))
  # Check R last if used
  if (any(unlist(lapply(p.map,function(x){any(x=="R") && x[length(x)]!="R"}))))
    stop("R factors must always be last")

  # Check responnses
  if ( type =="rd" ) {
    if (is.null(match.map))
      stop("Must specify supply a match.map for the DDM")
    if ( length(responses)!=2 )
      stop("DDM only applicable for two responses")
  }
  if (!is.null(responses2)) if (!is.character(responses2) || length(responses2)<2)
    stop("responses2 must be a character vector of length 2 or greater")

  # Check match.map (if supplied)
  if (!is.null(match.map)) {
    # Check structure
    if ( length(match.map)<1 || class(match.map[[1]]) != "list" )
      stop("match.map must be a list of lists")
    # Check match.map contains at least name M
    if ( !any(names(match.map) %in% "M") )
      stop("match.map must have a list named M")
    map.names <- names(match.map)[names(match.map)!="M"]
    map.levels <- unlist(lapply(match.map[names(match.map)!="M"],levels))
    # convert match.map$M to responses and check
    if ( is.numeric(unlist(match.map$M)) )
      match.map$M <- lapply(match.map$M,function(x){responses[x]})
    if ( !all(unlist(match.map$M) %in% responses) )
      stop("match.map$M has index or name not in response names")
    if ( !(all(sort(responses)==sort(unique(unlist(match.map$M))))) )
      stop("Not all response names are scored by match.map$M")
    if ( length(match.map)>1 &&
         !all(lapply(match.map[names(match.map)!="M"],class)=="factor") )
      stop("Entries in match.map besides M must be factors")
    if ( length(unlist(map.levels)) != length(unique(unlist(map.levels))) )
      stop("All match.map levels must be unqiue")
    # Check factors
    if ( any(names(factors) == "M") )
      stop("Do not use M as a factor name")
    if ( any(names(factors) %in% names(match.map)) )
      stop(paste(map.names,"used in match.map, can not use as a factor name"))
    if ( any(unlist(factors) %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as factor levels")
    if ( any(map.levels %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as match.map levels")
    if ( length(unlist(c(factors,map.levels))) !=
         length(unique(unlist(c(factors,map.levels)))) )
      stop("Factor levels cannot overlap match.map levels")
    # Check M and R are last
    if (any(unlist(lapply(p.map,function(x){any(x=="M") && x[length(x)]!="M"}))))
      stop("M factors must always be last")

  }

  factors.short <- factors
  factors$R <- responses
  if (!is.null(match.map)) factors$M <- c("true","false")

  # protect againt grep problems
  for (i in unlist(factors)) if ( length(grep(i,unlist(factors)))!=1 )
    stop("Factor, response or map level is not unique or is substring of another
         level or of \"true\" or \"false\"!" )

  # Add in extra match.map factors (if any)
  if ( !is.null(match.map) ) for (i in map.names)
    factors[[i]] <- levels(match.map[[i]])


  # Make parameter names
  names.par <- character(0)
  for ( i in names(p.map) )
  {
    if ( length(p.map[[i]])==1 && p.map[[i]] == "1" ) new.names <- i else
    {
      new.names <- paste(i,factors[[p.map[[i]][1]]],sep=".")
      if ( length(p.map[[i]])>1 ) for ( j in 2:length(p.map[[i]]) )
        new.names <- as.vector(outer(
          new.names,factors[[p.map[[i]][j]]],"paste",sep="."
        ))
    }
    names.par <- c(names.par,new.names)
  }

  # Make level array for manifest design and accumulators
  level.array <- make.level.array(factors[1:(length(factors.short)+1)])

  # Check match.map
  if ( !is.null(match.map) ) for ( i in names(match.map$M) ) {
    match.position <- grep(i,level.array)
    if ( length(match.position)==0 )
      stop(paste(i,"in match.map is not in the design"))
  }

  # Does the par use the match factor?
  is.M <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="M"))
  }))

  # Does the par use a response factor
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="R"))
  }))

  if ( type =="rd"  & ( any(is.M) | any(is.R) ) )
    stop("Cannot use M or R in DDM p.map")

  # Does the par use a map factor (i.e., in match map but not M)
  if ( !is.null(match.map) ) {
    is.map <- unlist(lapply(p.map,function(x){
      any(unlist(strsplit(x,".",fixed=T) %in% map.names))
    }))
  } else {
    is.map <- logical(length(p.map))
    names(is.map) <- names(p.map)
  }

  if ( any(is.map) ) {
    p.map.name <- lapply(p.map,function(x){
      unlist(strsplit(x,".",fixed=T))[
        unlist(strsplit(x,".",fixed=T)) %in% map.names]
    })
    nr <- length(responses)
    n <- length(level.array)
    map.shuffle <- matrix(aperm(array(1:n,dim=c(n/nr,nr,nr)),c(1,3,2)),ncol=nr)
  }

  if ( any(apply(cbind(is.M,is.R,is.map),1,function(x){sum(x)>1})) )
    stop("Parameters cannot have more than one of match.map and R factors")

  # use.par = boolean matrix for parameter use, cells x pars x resposnes
  use.par <- array(NA,
                   dim=c(length(level.array),length(names.par),length(responses)))
  dimnames(use.par) <-
    list(as.vector(level.array),names.par,responses)

  # col.par = column parameter type (1st name)
  if ( is.null(match.map) )
    col.par.levels <- responses else
      col.par.levels <- c(responses,"true","false",map.levels)

  col.par <- strsplit(dimnames(use.par)[[2]],".",fixed=T)
  col.fac <- lapply(col.par,function(x){x[-1]})
  col.par <- unlist(lapply(col.par,function(x){x[1]}))
  # split into fac and resp
  col.fac <- lapply(col.fac,function(x){
    if ( length(x)==0 ) out <- c(NA,NA)
    if ( length(x)==1 ) {
      if ( x %in% col.par.levels )
        out <- c(NA,x) else out <- c(x,NA)
    }
    if ( length(x)>1 )
      if ( x[length(x)] %in% col.par.levels )
        out <- c(paste(x[-length(x)],collapse="."),x[length(x)]) else
          out <- paste(x,collapse=".")
        out
  })
  col.resp <- unlist(lapply(col.fac,function(x){x[2]}))
  col.fac <- unlist(lapply(col.fac,function(x){x[1]}))

  row.fac <- strsplit(dimnames(use.par)[[1]],".",fixed=T)
  #  row.resp <- unlist(lapply(row.fac,function(x){x[length(x)]}))
  row.fac <- unlist(lapply(row.fac,function(x){
    paste(x[-length(x)],collapse=".")}))

  # Fill use.par array
  for ( p in unique(col.par) )
  { # parameters
    is.col <- p==col.par
    ncols <- sum(is.col)
    if ( ncols==1 ) use.par[,is.col,] <- TRUE else
    { # there are parameter subtypes
      for ( i in 1:ncols )
      { # each parameter subtype
        # select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE,dim(use.par)[1])
        if ( !is.na(tmp) ) is.fac.row[!grepl.dot(tmp,row.fac)] <- FALSE
        # set not applicable rows to false
        use.par[!is.fac.row,is.col,][,i,] <- FALSE
        if ( is.M[p] )
        { # has a match factor
          for ( j in names(match.map$M) )
          { # response cell
            correct.response <- match.map$M[[j]]
            is.rcell <- is.fac.row & grepl.dot(j,row.fac)
            for ( k in responses )
            { # responses
              if ( k==correct.response )
              {
                if ( grepl("true", col.resp[is.col][i]) )
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE

              } else {
                if ( grepl("false",col.resp[is.col][i]) )
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

  if ( any(is.na(use.par)) )
    stop("Some cells of the map were not assigned!")

  # add in constants
  all.par <- use.par[1,,1]
  all.par[1:length(all.par)] <- NA
  if ( length(constants)>0 ) {
    if ( !all(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }
  if ( !is.null(constant.prior) ) {
    if (!all(names(constant.prior) %in% names(all.par)))
      stop("Name(s) in constant prior not in p.map")
    all.par[names(constant.prior)] <- rprior(constant.prior)
  }

  attr(use.par, "all.par")  <- all.par
  attr(use.par, "p.vector") <- all.par[is.na(all.par)]

  if (length(attr(use.par,"p.vector"))<2)
    stop("DMC cant handle models with less than two parameters")

  if (verbose) {
    cat("\nParameter vector names (unordered) are: ( see attr(,\"p.vector\") )\n")
    print(names(all.par[is.na(all.par)]))
    cat("\nConstants are (see attr(,\"constants\") ):\n")
    print(constants)
    if (!is.null(constant.prior)) cat(paste("\nConstant prior for parameters:",
                                            paste(names(constant.prior),collapse=" "),"\n"))
    mod <- paste("\nModel type =", type)
    if (type == "norm") mod <- paste(mod, "(posdrift =", posdrift,")")
    cat(paste(mod,"\n\n"))
  }

  # parameters x cells x accumulators, used by p.list.dmc
  attr(use.par,"pca") <- aperm(use.par, c(2,1,3))

  # Names of parameter types (cannot have a period)
  attr(use.par,"par.names") <- unique(col.par)

  attr(use.par,"type") <- type
  attr(use.par,"factors") <- factors.short
  attr(use.par,"cvs") <- cvs
  attr(use.par,"responses") <- responses
  attr(use.par,"responses2") <- responses2
  attr(use.par,"constants") <- constants
  if (!is.null(constant.prior))
    attr(use.par,"constant.prior") <- constant.prior
  attr(use.par,"posdrift") <- posdrift

  par.df <- matrix(nrow=0,ncol=length(p.map))
  dimnames(par.df) <- list(NULL,names(p.map))
  par.df <- data.frame(par.df)
  if (!is.null(cvs))
    attr(par.df,"cvs") <- data.frame(matrix(rep(1,length(cvs)),
                                            nrow=1,dimnames=list(NULL,cvs)))

  # par.df <- try(transform.dmc(par.df),silent=TRUE)
  # if ( class(par.df)=="try-error" ) { # check list version
  #   par.df <- vector(mode="list",length=length(p.map))
  #   names(par.df) <- names(p.map)
  #   for (i in names(p.map)) par.df[[i]] <- matrix(nrow=0,ncol=length(responses))
  #   if (!is.null(cvs))
  #     attr(par.df,"cvs") <- data.frame(matrix(rep(1,length(cvs)),
  #       nrow=1,dimnames=list(NULL,cvs)))
  #   par.df <- try(transform.dmc(par.df),silent=TRUE)
  # }
  # if ( class(par.df)=="try-error")
  #   stop(paste("p.map must have names for each of the possible external parameter
  #     name (see top of model file for definitions)"))
  # attr(use.par,"internal.par.names") <- names(par.df)

  # save "n1" orders
  resp <- unlist(lapply(strsplit(level.array,".",fixed=TRUE),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp),ncol=nr)
  for (i in 1:length(resp))
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses],c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)

  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  if ( !is.null(match.map) ) for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell))
    match.cell[match.num[match.num %in%
                           grep(names(match.map$M[i]),names(match.cell))]] <- TRUE
  }

  attr(use.par, "n1.order") <- n1.order
  attr(use.par, "match.cell") <- match.cell
  if ( !is.null(match.map) ) attr(use.par, "match.map") <- match.map

  if (type=="rd") # r1 cells have z flipped
  {
    is.r1 <- rep(FALSE,length(row.names(use.par)))
    names(is.r1) <- row.names(use.par)
    is.r1[unlist(lapply(
      lapply(
        as.list(names(match.map$M)[match.map$M==responses[1]]),
        function(x)grepl(x,row.names(use.par))
      ),
      function(x) which(x==TRUE)))] <- TRUE
    attr(use.par,"is.r1") <- is.r1

    # add bound attributes
    bound <- rep("lower",dim(use.par)[1])
    bound[as.vector(sapply(
      paste("",names(attr(use.par,"match.map")$M),
            attr(use.par,"match.map")$M,sep="*"),
      function(x){grep(glob2rx(x),row.names(use.par))}))] <- "upper"
    names(bound) <- row.names(use.par)
    attr(use.par,"bound") <- bound
  }

  # CHECK NOT COMPATIBLE WITH p.list.dmc ONLY p.df.dmc
  #   if ( check.model ) { # Check model is good
  #     p.vector <- attr(use.par,"p.vector")
  #     p.vector[1:length(p.vector)] <- 1:length(p.vector)
  #     if (class(try(print.cell.p(p.vector,use.par,verbose=FALSE),silent=TRUE))=="try-error")
  #       stop("There is something wrong with the model that is not handled by the
  #            checks implemented so far, please send a bug report to
  #            andrew.heathcote@utas.edu.au.")
  #   }

  class(use.par) <- "model"
  return(use.par)
}

#' @rdname BuildModel
#' @export
print.model <- function(model, p.vector = NULL) {

  if (is.null(p.vector)) {
    nr <- length(attr(model, "response"))
    for (i in 1:nr) {
      dim3 <- dimnames(model)[[3]]
      cat(dim3[i], "\n")
      print(model[,, i])
    }
    message("model has the following attributes: ")
    print(names(attributes(model)))
    return(invisible(model))
  } else {
    dim1 <- dimnames(model)[[1]]

    out <- lapply(dim1, function(x) {
      print(x)
      print(p.df.dmc(p.vector, x, model))
    })
    return(invisible(out))
  }
}

