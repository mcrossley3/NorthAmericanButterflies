
setupF <- function(fit, xvar, call.env, data) {
  CALL <- if (isS4(fit)) fit@call else fit$call
  if (!is.null(data)) {
    Data <- data
  } else if (!is.null(CALL) && ('data' %in% names(CALL)) && exists(as.character(CALL$data), call.env)) {
    env <- call.env
    Data <- eval(CALL$data, envir=env)
  } else if (isS4(fit)) {
    FRAME <- try(fit@frame, silent=TRUE)
    DATA <- try(fit@data, silent=TRUE)
    if (!inherits(DATA, 'try-error')) {
      Data <- DATA
    } else if (!inherits(FRAME, 'try-error')) {
      Data <- FRAME
    } else {
      stop("visreg cannot find the data set used to fit your model; supply it using the 'data=' option", call.=FALSE)
    }
  } else {
    ENV <- environment(fit$terms)
    if ("data" %in% names(fit) && is.data.frame(fit$data)) {
      Data <- fit$data
      env <- NULL
    } else if (is.null(CALL$data)) {
      env <- NULL
      Data <- NULL
    } else if (exists(as.character(CALL$data), call.env)) {
      env <- call.env
      Data <- eval(CALL$data, envir=env)
    } else if (exists(as.character(CALL$data), ENV)) {
      env <- ENV
      Data <- eval(CALL$data, envir=ENV)
    } else {
      stop("visreg cannot find the data set used to fit your model; supply it using the 'data=' option", call.=FALSE)
    }
  }
  form <- formula(fit)
  if (!is.null(Data)) names(Data) <- gsub('offset\\((.*)\\)', '\\1', names(Data))
  if (inherits(fit, 'mlm') && fit$terms[[2L]] != 'call') {
    ff <- form
    ff[[2]] <- NULL
    av <- get_all_vars(ff, Data)      # If mlm with matrix as Y, outside of data frame framework
  } else {
    av <- get_all_vars(form, Data)    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14905
  }
  f <- as.data.frame(av)

  if (inherits(CALL$random, "call")) {
    rf <- as.data.frame(as.list(get_all_vars(CALL$random, Data)))
    rf <- rf[, setdiff(names(rf), names(f)), drop=FALSE]
    f <- cbind(f, rf)
  }
  if ("subset" %in% names(CALL) & !(inherits(fit, 'averaging'))) {
    s <- CALL$subset
    subset <- eval(substitute(s), Data, env)
    f <- f[which(subset==TRUE), , drop=FALSE]
  }
  suppressWarnings(f <- f[!apply(is.na(f), 1, any), , drop=FALSE])

  # Handle some variable type issues
  needsUpdate <- FALSE
  f <- droplevels(f)
  frameClasses <- sapply(f, class)
  if (any(frameClasses=="Surv")) needsUpdate <- TRUE
  if (any(frameClasses=="character")) {
    needsUpdate <- TRUE
    for (j in 1:ncol(f)) if (typeof(f[,j])=="character") f[,j] <- factor(f[,j])
  }
  if (any(frameClasses=="logical")) {
    needsUpdate <- TRUE
    for (j in 1:ncol(f)) if (typeof(f[,j])=="logical") f[,j] <- as.double(f[,j])
  }
  if (missing(xvar)) {
    all_x <- strsplit(parseFormula(formula(fit)[3]), ' + ', fixed=TRUE)[[1]]
    inModel <- sapply(names(f), function(x) x %in% all_x)
    const <- sapply(f, function(x) all(x==x[1]))
    xvar <- names(f)[!const & inModel]
  }
  if (length(xvar)==0) stop("The model has no predictors; visreg has nothing to plot.", call.=FALSE)
  for (i in 1:length(xvar)){if (!is.element(xvar[i], names(f))) stop(paste(xvar[i], "not in model"), call.=FALSE)}

  attr(f, "needsUpdate") <- needsUpdate
  attr(f, "xvar") <- xvar
  f
}


##################################
setupCond <- function(cond, f, by, breaks) {
  for (i in seq_along(cond)) if(!is.character(cond[[i]]) & is.factor(f[, names(cond)[i]])) cond[[i]] <- as.character(cond[[i]])
  
  if (missing(by)) {
    cond <- list(cond)
  } else {
    cond.orig <- cond
    if(is.numeric(f[, by])) {
      if (length(breaks)==1) {
        unique.by <- unique(f[, by])
        if (breaks >= length(unique.by)) {
          lev <- sort(unique.by)
        } else {
          a <- 1/5/2^(breaks-2)
          lev <- as.double(quantile(f[, by], seq(a, 1-a, length=breaks), type=1))
        }
      } else {
        lev <- breaks
      }
      n.by <- length(lev)
    } else {
      if (is.factor(breaks) || is.character(breaks)) {
        if (!all(breaks %in% levels(f[, by]))) stop("'breaks' does not match levels of 'by' variable", call.=FALSE)
        lev <- breaks
      } else {
        lev <- levels(f[, by])
      }
      n.by <- length(lev)
    }
    
    cond <- vector("list", n.by)
    for (i in 1:n.by) {
      a <- if (is.factor(lev)) as.character(lev[i]) else lev[i]
      cond[[i]] <- c(a, cond.orig)
      names(cond[[i]])[1] <- by
      cond[[i]] <- as.list(cond[[i]])
    }
    attr(cond, "lev") <- lev
  }
  cond
}


################################
makeYName <- function(fit, scale, trans, type) {
  if (scale=="response" | (class(fit)[1] %in% c("lm", "aov") & identical(trans, I))) {
    if (type=="contrast") {
      yName <- as.expression(substitute(list(Delta) * x, list(x=as.character(formula(fit)[2]))))
    } else {
      yName <- as.character(formula(fit)[2])
    }
  } else if (inherits(fit, "mlm")) {
    if (type=="contrast") {
      yName <- sapply(colnames(fit$residuals), function(y) {as.expression(substitute(list(Delta) * x, list(x=y)))})
    } else {
      yName <- colnames(fit$residuals)
    }
  } else if (inherits(fit, "randomForest")) {
    if (fit$type=="regression") yName <- as.character(formula(fit)[2])
    if (fit$type=="classification") yName <- paste0("Pr(", as.character(formula(fit)[2]), ")")
  } else {
    yName <- NULL
  }
  yName
}


###################################
# v is a list of three elements: fit, res, and meta
# alternatively (class "visregList"), a list of visreg elements

setupV <- function(fit, f, xvar, nn, cond, type, trans, alpha=0.05, by, yName) {
  
  # Initial setup
  if (length(xvar) > 1 & length(cond) > 1) stop("Cannot specify 'by' and multiple x variables simultaneously", call.=FALSE)
  J <- max(length(xvar), length(cond))
  Attempt <- try(max(attr(terms(as.formula(formula(fit))), "order")) > 1, silent=TRUE)
  hasInteraction <- ifelse(inherits(Attempt, 'try-error'), FALSE, Attempt)
  lev <- attr(cond, "lev")

  # Get xy list
  xy <- vector("list", J)
  for (j in 1:J) {
    cond.j <- if (length(cond) > 1) cond[[j]] else cond[[1]]
    name <- if (length(xvar) > 1) xvar[j] else xvar
    xy[[j]] <- getXY(fit, f, name, nn, cond.j, type, trans, alpha=0.05)
  }
  if (!missing(by)) xy <- subsetV(xy, f, by, lev, type)

  # Format
  meta <- list(x=xvar, y=xy[[1]]$y$name, hasInteraction=hasInteraction, yName=yName, trans=trans, class=class(fit))
  K <- xy[[1]]$y$n
  if (K==1) {
    if (!missing(by)) {
      meta$by <- by
      v <- list(fit=NULL, res=NULL, meta=meta)
      for (j in 1:length(xy)) {
        fit.j <- data.frame(xy[[j]]$x$DD, visregFit=xy[[j]]$y$fit, visregLwr=xy[[j]]$y$lwr, visregUpr=xy[[j]]$y$upr)
        res.j <- data.frame(xy[[j]]$x$D, visregRes=xy[[j]]$y$r, visregPos=xy[[j]]$y$pos)
        fit.j[, xvar] <- xy[[j]]$x$xx
        res.j[, xvar] <- xy[[j]]$x$x
        v$fit <- rbind(v$fit, fit.j)
        v$res <- rbind(v$res, res.j)
      }
      class(v) <- "visreg"
    } else {
      v <- vector("list", J)
      for (j in 1:J) {
        meta.j <- meta
        meta.j$x <- xvar[j]
        v[[j]] <- list(fit=data.frame(xy[[j]]$x$DD, visregFit=xy[[j]]$y$fit, visregLwr=xy[[j]]$y$lwr, visregUpr=xy[[j]]$y$upr),
                       res=data.frame(xy[[j]]$x$D, visregRes=xy[[j]]$y$r, visregPos=xy[[j]]$y$pos),
                       meta=meta.j)
        v[[j]]$fit[, xvar[j]] <- xy[[j]]$x$xx
        v[[j]]$res[, xvar[j]] <- xy[[j]]$x$x
        class(v[[j]]) <- "visreg"
      }
      if (J==1) {
        v <- v[[1]]
      } else {
        class(v) <- "visregList"
      }
    }
  } else {
    if (!missing(by)) {
      meta$by <- by
      v <- vector("list", K)
      for (k in 1:K) {
        meta.k <- meta
        meta.k$y <- meta$y[k]
        meta.k$yName <- meta$yName[k]
        v[[k]] <- list(fit=NULL, res=NULL, meta=meta.k)
        for (j in 1:J) {
          fit.jk <- data.frame(xy[[j]]$x$DD, visregFit=xy[[j]]$y$fit[,k], visregLwr=xy[[j]]$y$lwr[,k], visregUpr=xy[[j]]$y$upr[,k])
          res.jk <- data.frame(xy[[j]]$x$D, visregRes=xy[[j]]$y$r[,k], visregPos=xy[[j]]$y$pos[,k])
          fit.jk[, xvar] <- xy[[j]]$x$xx
          res.jk[, xvar] <- xy[[j]]$x$x
          v[[k]]$fit <- rbind(v[[k]]$fit, fit.jk)
          v[[k]]$res <- rbind(v[[k]]$res, res.jk)
        }
        class(v[[k]]) <- "visreg"
      }
      class(v) <- "visregList"
    } else {
      v <- vector("list", J*K)

      for (j in 1:J) {
        for (k in 1:K) {
          meta.jk <- meta
          meta.jk$x <- meta$x[j]
          meta.jk$y <- meta$y[k]
          meta.jk$yName <- meta$yName[k]
          l <- (j-1)*K + k
          v[[l]] <- list(fit=data.frame(xy[[j]]$x$DD, visregFit=xy[[j]]$y$fit[,k], visregLwr=xy[[j]]$y$lwr[,k], visregUpr=xy[[j]]$y$upr[,k]),
                         res=data.frame(xy[[j]]$x$D, visregRes=xy[[j]]$y$r[,k], visregPos=xy[[j]]$y$pos[,k]),
                         meta=meta.jk)
          v[[l]]$fit[, xvar[j]] <- xy[[j]]$x$xx
          v[[l]]$res[, xvar[j]] <- xy[[j]]$x$x
          class(v[[l]]) <- "visreg"
        }
      }
      class(v) <- "visregList"
    }
  }
  v
}


###################################
getXY <- function(fit, f, name, nn, cond, type, trans, alpha=0.05) {
  if (type=="conditional") {
    x <- setupD(fit, f, name, nn, cond)
    y <- Response(fit, x, trans, alpha)
  } else if (type=="contrast") {
    x <- setupX(fit, f, name, nn, cond)
    y <- Terms(fit, f, x, trans, alpha)
    x <- setupD(fit, f, name, nn, cond)
  } 
  list(x=x, y=y)
}


##################################
setupD <- function(fit, f, name, nn, cond, whitespace) {
  ## Set up n-row data frame for residuals
  x <- f[, name]
  xdf <- data.frame(x)
  names(xdf) <- name
  df <- fillFrame(f, xdf, cond)

  rhs_form <- formula(fit)
  rhs_form[2] <- NULL
  simple_rhs <- paste(all.vars(rhs_form), collapse = ' + ')
  simple_form <- as.formula(paste("~", simple_rhs))
  D <- model.frame(simple_form, df)
  condNames <- setdiff(names(D), name)
  condNames <- intersect(condNames, names(df))
  D <- cbind(D, df[, setdiff(names(df), names(D)), drop=FALSE])

  ## Set up nn-row data frame for prediction
  dots <- list()
  if (is.factor(x)) {
    xx <- factor(levels(x), levels=levels(x))
  } else {
    xx <- seq(min(x), max(x), length=nn)
  }
  xxdf <- data.frame(xx)
  names(xxdf) <- name
  df <- fillFrame(f, xxdf, cond)
  DD <- model.frame(simple_form, df)
  DD <- cbind(DD, df[, setdiff(names(df), names(DD)), drop=FALSE])

  list(x=x, xx=xx, D=D, DD=DD, factor=is.factor(x), name=name, cond=D[1, condNames, drop=FALSE])
}


##############################################
fillFrame <- function(f, x, cond) {
  ## x  = data frame of x variable(s) being changed
  ## x2 = variables being filled by median
  ## x3 = variables specified by cond
  if (missing(cond)) cond <- NULL
  for (j in 1:ncol(f)) {
    if (is.factor(f[,j]) && !is.element(names(f)[j], names(cond)) && !is.element(names(f)[j], names(x))) {
      mode = names(sort(-table(f[j])))[1]
      eval(parse(text=c('cond=c(cond, list(', names(f)[j],'=factor(mode, levels=levels(f[, names(f)[j]]))))')))
    }
  }
  exclude <- c(names(x), names(cond), names(which(sapply(f, function(x) inherits(x, "Surv")))))
  x2 <- lapply(as.data.frame(f[, setdiff(names(f), exclude)]), median)
  names(x2) <- setdiff(names(f), exclude)
  x3 <- cond
  for (j in seq_along(x3)) {
    if (is.character(x3[[j]])) x3[[j]] <- factor(x3[[j]], levels=levels(f[, names(x3)[j]]))
  }

  if (length(x2)>0 & length(x3)>0) newdf <- data.frame(x, x2, x3, check.names=FALSE)
  else if (length(x2)>0 & length(x3)==0) newdf <- data.frame(x, x2, check.names=FALSE)
  else if (length(x2)==0 & length(x3)>0) newdf <- data.frame(x, x3, check.names=FALSE)
  else newdf <- x

  newdf
}


#####################################
Response <- function(fit, x, trans, alpha=0.05) {

  ## Calculate partial residuals
  rr <- visregResid(fit)
  nr <- if (is.matrix(rr)) nrow(rr) else length(rr)
  if (nr>0 && nrow(x$D) != nr) warning("Residuals do not match data; have you changed the original data set?  If so, visreg is probably not displaying the residuals for the data set that was actually used to fit the model.")
  y <- visregPred(fit, x$D)
  if (is.null(rr)) {
    r <- NULL
  } else {
    r <- y + rr
  }

  # Calculate predictions
  p <- visregPred(fit, x$DD, se.fit=TRUE)

  ## Format output
  if (inherits(p, "svystat")) {
    p <- list(fit=as.double(p), se.fit=sqrt(attr(p,"var")))
  } else if (inherits(fit, "rq")) {
    p <- list(fit=as.double(p[,1]), se.fit=as.double(p[,3]-p[,2])/(2*qnorm(.975)))
  } else if (inherits(fit, "rms")) {
    p$fit <- p$linear.predictors
  } else if (is.double(p)) {
    p <- list(fit=p, se.fit=NA)
  }
  m <- ifelse(identical(class(fit), "lm"), qt(1-alpha/2, fit$df.residual), qnorm(1-alpha/2))
  upr <- p$fit + m*p$se.fit
  lwr <- p$fit - m*p$se.fit
  if (is.matrix(p$fit)) {
    if (length(r)==0) {
      R <- matrix(NA, nrow(x$D), ncol=ncol(p$fit))
    } else {
      R <- matrix(trans(r), ncol=ncol(p$fit))
    }
    val <- list(fit=matrix(trans(p$fit), ncol=ncol(p$fit)), lwr=matrix(trans(lwr), ncol=ncol(p$fit)), upr=matrix(trans(upr), ncol=ncol(p$fit)), r=R)
    val$name <- colnames(val$fit) <- colnames(p$fit)
  } else {
    if (length(r)==0) r <- rep(NA_real_, nrow(x$D))
    val <- list(fit=as.double(trans(p$fit)), lwr=as.double(trans(lwr)), upr=as.double(trans(upr)), r=as.double(trans(r)), name=as.character(formula(fit)[2]))
  }
  val$pos <- rr>0
  if (length(val$pos)==0) {
    if (is.matrix(p$fit)) {
      val$pos <- matrix(NA, nrow(x$D), ncol(p$fit))
    } else {
      val$pos <- rep(NA_real_, nrow(x$D))
    }
  }
  val$n <- if (is.matrix(p$fit)) ncol(p$fit) else 1
  val
}


######################################
visregResid <- function(fit) {
  if (inherits(fit, "randomForest")) {
    if (fit$type=="regression") rr <- fit$y - fit$predicted
    if (fit$type=="classification") {
      P <- predict(fit, type="prob")
      rr <- (fit$y==colnames(P)[2]) - P[,2]
    }
  } else if (inherits(fit, 'coxph')) {
    rr <- residuals(fit, type='deviance')
  } else if (inherits(fit, 'gamlss')) {
    rr <- residuals(fit, what='mu')
  } else {
    rr <- residuals(fit)
  }
  if (!is.matrix(rr) & length(rr)>0) rr <- rr[!is.na(rr)]
  rr
}


##############################
visregPred <- function(fit, Data, se.fit=FALSE) {
  predict.args <- list(object=fit, newdata=Data)
  if (inherits(fit, "lme")) predict.args$level <- 0
  if (inherits(fit, "merMod")) predict.args$re.form <- NA
  if (inherits(fit, "rq")) predict.args$interval <- "confidence"
  if (inherits(fit, "svm")) predict.args$probability <- TRUE
  if (inherits(fit, "multinom") | inherits(fit, "polr")) predict.args$type <- "probs"
  if (inherits(fit, "gbm")) predict.args$n.trees <- length(fit$trees)
  if (inherits(fit, "betareg")) predict.args$type <- "link"
  dots <- list()
  if (length(dots)) predict.args[names(dots)] <- dots

  if (se.fit) {
    if (inherits(fit, "mlm")) {
      p <- list(fit = suppressWarnings(do.call("predict", predict.args)), se.fit = se.mlm(fit, newdata=Data))
    } else if (inherits(fit, "randomForest") && fit$type=="classification") {
      predict.args$type <- "prob"
      P <- suppressWarnings(do.call("predict", predict.args))
      p <- list(fit=P[,2], se.fit=NA)
    } else if (inherits(fit, "loess")) {
      predict.args$se <- TRUE
      p <- suppressWarnings(do.call("predict", predict.args))
    } else {
      predict.args$se.fit <- TRUE
      p <- suppressWarnings(do.call("predict", predict.args))
    }
  } else {
    if (inherits(fit, "randomForest") && fit$type=="classification") {
      p <- predict(fit, type="prob")[,2]
    } else if (inherits(fit, 'rq')) {
      p <- suppressWarnings(do.call("predict", predict.args))[,1]
    } else {
      p <- suppressWarnings(do.call("predict", predict.args))
    }
  }
  if (inherits(fit, "svm") && fit$type < 3) p <- attr(p, "probabilities")
  p
}


################################
## Subsets xy so that residuals appear only once
## Should probably be renamed subsetXY
subsetV <- function(v, f, by, lev, type) {
  ## Calculate distance
  if (is.numeric(f[, by])) {
    D <- matrix(NA, nrow(f), length(v))
    for (i in 1:length(v)) {
      D[,i] <- (f[, by]-lev[i])^2
    }
  }
  for (i in 1:length(v)) {
    if (is.factor(f[, by])) {
      ind <- as.character(f[, by])==as.character(lev[i])
    } else {
      ind <- (apply(D, 1, which.min)==i)
    }
    v[[i]]$x$D <- v[[i]]$x$D[ind,]
    v[[i]]$x$x <- v[[i]]$x$x[ind]
    v[[i]]$y$r <- if (v[[i]]$y$n == 1) v[[i]]$y$r[ind] else v[[i]]$y$r[ind,]
    v[[i]]$y$pos <- if (v[[i]]$y$n == 1) v[[i]]$y$pos[ind] else v[[i]]$y$pos[ind,]
  }
  v
}


###################################
# setupV for visreg2d
# Returns a list of x, y, and z for plotting
setupV2 <- function(fit, f, xvar, yvar, nn, cond, type, scale, trans) {
  n.z <- if (inherits(fit, "mlm")) ncol(coef(fit)) else 1
  form <- parseFormula(formula(fit)[3])
  x <- f[, xvar]
  y <- f[, yvar]
  xx <- if(is.factor(x)) factor(levels(x), levels=levels(x)) else seq(min(x), max(x), length=nn)
  yy <- if(is.factor(y)) factor(levels(y), levels=levels(y)) else seq(min(y), max(y), length=nn)
  xydf <- as.data.frame(expand.grid(xx, yy))
  names(xydf) <- c(xvar, yvar)

  if (type=="conditional") {
    df <- fillFrame(f, xydf, cond)
    DD <- model.frame(as.formula(paste("~", form)), df)
    DD <- cbind(DD, df[, setdiff(names(df), names(DD)), drop=FALSE])
    P <- predict(fit, newdata=DD)
    if (inherits(fit, "mlm")) {
      z <- vector("list", n.z)
      for (i in 1:n.z) z[[i]] <- matrix(trans(P[,i]), nrow=length(xx), ncol=length(yy))
    } else {
      z <- matrix(trans(P), nrow=length(xx), ncol=length(yy))
    }
  } else if (type=="contrast") {
    xref <- if(is.factor(x)) xx[1] else xref <- mean(x)
    yref <- if(is.factor(y)) yy[1] else yref <- mean(y)
    xydf <- rbind(c(xref, yref), xydf)
    df <- fillFrame(f, xydf, cond)
    DD <- rbind(f, df)
    if (inherits(fit, "mlm")) {
      ind <- apply(is.finite(coef(fit)), 1, all)
      if (!identical(ind, apply(is.finite(coef(fit)), 1, any))) stop("Inconsistent NA/NaN coefficients across outcomes", call.=FALSE)
    } else ind <- is.finite(coef(fit))
    XX. <- model.matrix(as.formula(paste("~", formula(fit)[3])), DD)[-(1:nrow(f)), ind]
    XX <- t(t(XX.[-1,])-XX.[1,])
    if (inherits(fit, "mlm")) {
      z <- vector("list", n.z)
      for (i in 1:n.z) z[[i]] <- matrix(trans(XX%*%coef(fit)[ind, i]), nrow=length(xx), ncol=length(yy))
    } else {
      z <- matrix(trans(XX%*%coef(fit)[ind]), nrow=length(xx), ncol=length(yy))
    }
  }
  zname <- makeYName(fit, scale, trans, type)
  D <- model.frame(as.formula(paste("~", form)), df)
  condNames <- setdiff(names(D), c(xvar, yvar))
  condNames <- intersect(condNames, names(df))
  baseMeta <- list(x=xvar, y=yvar, trans=trans, class=class(fit), cond=D[1, condNames, drop=FALSE])

  if (n.z > 1) {
    v <- vector("list", n.z)
    for (i in 1:n.z) {
      meta <- baseMeta
      meta$z <- zname[i]
      v[[i]] <- list(x=xx, y=yy, z=z[[i]], meta=meta)
      class(v[[i]]) <- 'visreg2d'
    }
    class(v) <- 'visregList'
  } else {
    meta <- baseMeta
    meta$z <- zname
    v <- list(x=xx, y=yy, z=z, meta=meta)
    class(v) <- 'visreg2d'
  }
  v
}


#############################
parseFormula <- function(form) {
  form <- gsub("\\|", "\\+", as.character(form))
  form <- gsub("\\*", "\\+", as.character(form))
  f <- if (grepl("\\+", form)) unlist(strsplit(form, "\\+")) else form
  n.f <- length(f)
  for (i in 1:n.f) {
    f[i] <- gsub(" ", "", f[i])
    if (substr(f[i], 1, 3) %in% c("te(", "ti(", "lp(") || substr(f[i], 1, 2) %in% c("s(")) {
      matched <- gregexpr("\\((?>[^()]|(?R))*\\)", f[i], perl = T)
      fi <- substring(f[i], matched[[1]]+1, matched[[1]] + attr(matched[[1]], "match.length")-2)
      fi <- gsub("\\([^\\)]+.*\\)", "", fi)
      fi <- unlist(strsplit(fi, ","))
      fi <- fi[!grepl("=", fi) | grepl("by\\s*=", fi)]
      fi <- gsub("by\\s*=", "", fi)
      f[i] <- paste(fi, collapse=" + ")
    }
    f[i] <- gsub("\\b[^\\(]*\\(([^,]+).*\\)", "\\1", f[i])
  }
  f <- f[f!=""]
  val <- paste(f, collapse=" + ")
  val <- gsub("()", "", val)
  val
}
