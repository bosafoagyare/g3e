#%%%%%%%%%%%%%%%%%%%%%%% Function to be Exported %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @name g3e
#'
#' @title Generalized Expectile Estimating Equations
#'
#' @description \code{ge3} obtains the parameter estimates for generalized
#'   expectile estimating equations for repeated measure data
#'
#' @param id 
#' @param y 
#' @param x 
#' @param expectile 
#' @param intercept 
#' @param corstr 
#' @param scaled 
#'
#' @return
#' @export
#'
#@examples
g3e <- function(id, y, x, expectile, intercept='FALSE', corstr, scaled){

  ## remove NAs if exist
  complete_index <- !is.na(y) 
  id <- id[complete_index]
  y  <- y[complete_index]
  x  <- x[complete_index, ]
  
  
  if (intercept == 'TRUE') x = cbind(1, x)

  ltps = split(1:length(id), id)
  tpslengvec = unlist(lapply(ltps, length))
  sampsize = length(ltps)
  xncol = ncol(x)
  expectleng = length(expectile)

  listbeta = list()
  leres = list()
  lsandCovhat = list()
  lmodCovhat = list()
  lrho = list()
  lsig = list()
  lQlike = list()
  lQIC = list()

  if (corstr=="ind"){

    for (k in 1:expectleng){
      # elmfit = lm(y~x-1)
      elmfit <- geepack::geeglm(y~x[, -1], id=id, corstr=corstr)
      eres = elmfit$residuals
      betaOld = elmfit$coef
      deltaParam = 1
      iter = 1
      while( deltaParam > 1e-05) {
        firstDeriv = 0
        scoreFun = 0
        for(i in 1:sampsize ){
          itps = ltps[[i]]
          ires = eres[itps]
          xi = x[itps, ]
          iphi = asymweight(ires, expectile[k])
          txiPhi = t(xi) %*% diag(iphi)
          firstDeriv = firstDeriv + txiPhi %*% xi
          scoreFun = scoreFun + txiPhi %*% ires
        }
        betaNew = betaOld + MASS::ginv(firstDeriv) %*% scoreFun
        eres = y - x %*% betaNew
        betaDiff = abs(betaNew-betaOld)
        deltaParam = max(abs((betaNew - betaOld)/(betaOld + .Machine$double.eps)))  #  max(betaDiff)
        betaOld = betaNew
        iter = iter + 1
      }

      listbeta[[k]] = c(betaNew)
      leres[[k]] = c(eres)
      sigparamhat = sigparam(eres, expectile[k], xncol)
      rhoparam = 0
      lsig[[k]] = sigparamhat
      VarCovlist = listCov(x, eres, expectile[k], sigparamhat, rhoparam, ltps, sampsize, corstr, scaled)
      lsandCovhat[[k]] = VarCovlist$sandCovhat
      lmodCovhat[[k]] = VarCovlist$ModCovhat
    }

    out = list(id=id, y=y, x=x, expectile=expectile, listbeta=listbeta, lmodCovhat=lmodCovhat,
               lsandCovhat=lsandCovhat, lsig=lsig, lrho=lrho, leres=leres)


  } else {

    vdof = dof(tpslengvec,  xncol, corstr)

    for (k in 1:expectleng){

      # elmfit = lm(y~x-1)
      elmfit <- geepack::geeglm(y~x[ ,-1], id=id, corstr=corstr)
      eres = elmfit$residuals
      sigparamhat = 1 # sum(eres^2)/(length(eres) - xncol)
      rhoparamOld = 0.5 # (length(eres) - xncol)/vdof
      betaOld = elmfit$coef
      deltaParam = 1
      iter = 1

      while( deltaParam > 1e-05) {
        firstDeriv = 0
        scoreFun = 0

        for(i in 1:sampsize ){
          itps = ltps[[i]]
          ires = eres[itps]
          xi = x[itps, ]
          iphi = asymweight(ires, expectile[k])
          ivcovmat = ivcov(sigparamhat, rhoparamOld, length(itps), corstr, scaled)
          ivcovmat = MASS::ginv(ivcovmat)
          corEprx = t(xi) %*% (ivcovmat %*% diag(iphi))
          firstDeriv = firstDeriv + corEprx %*% xi
          scoreFun = scoreFun + corEprx %*% ires
        }

        betaNew = betaOld + MASS::ginv(firstDeriv) %*% scoreFun
        eres = y - x %*% betaNew
        sigparamhat = sigparam(eres, expectile[k], xncol)
        rhoparamNew = sum(unlist(lapply(1:sampsize, function(i)irhoparam(eres[ltps[[i]]], expectile[k], corstr))))/(vdof * sigparamhat)
        betaDiff = abs((betaNew - betaOld)/(betaOld + .Machine$double.eps))
        rhoDiff = abs((rhoparamNew-rhoparamOld)/(rhoparamOld + .Machine$double.eps))
        deltaParam = max(c(betaDiff, rhoDiff))  #  max(betaDiff) #max(c(betaDiff, rhoDiff))
        betaOld = betaNew
        rhoparamOld = rhoparamNew
        iter = iter + 1

      }
      lrho[[k]] = rhoparamOld
      lsig[[k]] = sigparamhat
      listbeta[[k]] = c(betaNew)
      leres[[k]] = c(eres)
      VarCovlist = listCov(x, eres, expectile[k], sigparamhat, rhoparamNew, ltps, sampsize, corstr, scaled)
      lsandCovhat[[k]] = VarCovlist$sandCovhat
      lmodCovhat[[k]] = VarCovlist$ModCovhat
    }

    out = list(id=id, y=y, x=x, expectile=expectile, listbeta=listbeta, lmodCovhat=lmodCovhat,
               lsandCovhat=lsandCovhat, lsig=lsig, lrho=lrho, leres=leres)

  }


  out

}
