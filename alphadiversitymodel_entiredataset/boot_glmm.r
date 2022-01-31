###R script written by Roger Mundry for calculating confidence intervals.

#' @title Parametric bootstraps of confidence intervals of model estimates and potentially fitted values for
#' (G)LMMs fitted with lmer or glmmTMB functions.
#' @description Derive Bootstrapped confidence intervals of model estimates and potentially fitted values for
#' (G)LMMs fitted with lmer or glmmTMB functions and using a parametric bootstrap.
#' @param m fitted model object as obtained from \code{lmer}, \code{glmer}, \code{glmer.nb}, or \code{glmmTMB}.
#' @param data data frame comprising all variables needed to fit the model (only \code{boot.glmmtmb}).
#' @param discard.warnings logical (function \code{boot.lmer}); determines whether bootstraps that failed
#' with an error or issued a warning should be discarded (\code{T}; a 'dangerous' option; see datails) or not (\code{F}; default).
#' @param discard.non.conv logical (function \code{boot.glmmtmb}); determines whether bootstraps that failed
#' with an error should be discarded (\code{T}; a 'dangerous' option; see datails) or not (\code{F}; default).
#' @param nboots number of bootstraps to be conducetd (default: 1000).
#' @param para logical; determines whether bootstraps should be parallellized (\code{T}) or not (\code{F};
#' default; see details).
#' @param n.cores number cores to be used when \code{para=T}. Can be one of "all-1" (the default), "all",
#' or an integer.
#' @param resol resolution with which confidence limits of fitted values for the effects of covariates should be
#' determined (default: 1000; see details).
#' @param level confidence level (default: 0.95).
#' @param use, character vector or list of character vectors indicating the names of the predictors for which
#' confidence intervals of fitted values should be obtained. Defaults to \code{NULL} in which case no
#' confidence intervals of fitted values are returned (see details).
#' @param circ.var.name character; name of a potential circular variable present in the model (only relevant
#' when \code{use} is not left at its default).
#' @param circ.var numeric vector with the values of the circular variable as handed over to the model fitting
#' function.
#' @param use.u logical; determines whether bootstraps should be conditional on the particular levels of the
#' random effects present in the data \code{T} or not (\code{F}; default; see \code{\link{bootMer}} for details).
#' @param save.path path in which results of individual bootstraps should be saved (defaults to \code{NULL} in
#' which case results of individual bootstraps aren't saved).
#' @param load.lib logical determining whether the package needed to conduct the bootstrap should be loaded
#' (\code{T}) or not (\code{F}; default).
#' @param lib.loc path from where the package should be loaded (defaults to \code{.libPaths()}; usually one
#' dosn't need to bother about; see details).
#' @param set.all.effects.2.zero don't know what this is doing (just leave it at its default; \code{F}).
#' @details
#'  Both functions are actually wrappers of the functions \code{bootMer{lme4}} and
#' \code{simulate.glmmTMB{glmmTMB}}, respectively, making their use more convenient, allowing to obtain
#' bootstrapped confidence intervals of fitted values, keep track of warnings issued, and avoid boostraps
#' failing with an error.
#'
#' Note that these functions try to fit many models (at least as many as indicated to the argument
#' \code{n.boots}). As a consequence, it can take a while until they finish. For instance, for a complex GLMM
#' with, e.g., binomial, negative binomial, Poisson, or beta error distribution, it can easily take days
#' for the function to complete, even when it runs parallelized on on, say, 16 or 64 cores.
#'
#' Use parallelisation (i.e., \code{para=T}) on dedicated machines, and I suggest to not use it on laptops to avoid
#' over heating. The default for the number cores (\code{n.cores="all-1"}) keeps some capacity for the OS,
#' allowing to easier have a look at what's going on.
#'
#' Excluding bootstraps which issued a warning or failed with an error (\code{discard.warnings=T} or
#' \code{discard.non.conv=T}, respectively) is pretty 'dangerous' since in case this likely happens many models will
#' tried to be fitted and then discarded. Hence, getting, say, 1,000 bootstraps may take ages...
#'
#' When the argument \code{use} is not left at its default, the function will return confidence intervals
#' of fitted values. When a vector is handed over to \code{use} which includes a single covariate and the model doesn't
#' include any factors, the function will return a data frame with one column for each term in the model with all
#' except the intercept and the particular covariate(s) set to zero, as well as columns with the fitted value
#' ("fitted") and its confidence limits ("lwr", "upr"). The number rows will be equal to \code{resol} whereby the
#' covariate will range from its minimum to its maximum with as many values as \code{resol}, equally spaced.
#' When the vector handed over to \code{use} comprises two or more covariates, confidence intervals of fitted
#' values will be determined for each combination of their values (i.e., all will range from their minimum to
#' their maximum with \code{resol} values and the resulting data frame will have \code{resol}^n rows, with n
#' being the numbe of covariates). When the model
#' comprises factors, these will automatically be included in each vector handed over to use. For instance, if a model
#' comprises the factors 'sex' and 'species' and a covariate 'age' and one runs the function with use="sex", then one
#' will get fitted values with confidence intervals for each combination of sex, species, and age.
#'
#' When a list of vectors is handed over to \code{use} then a data frame ass described in the previous para is generated
#' for each element of the list.
#'
#' When the model comprises a circular variable fitted by including its sine and cosine into the model
#' (Stolwijk et al. 1999) and one wants to get confidence intervals of fitted values, then it is important
#' that the vector with the circular variable (as handed over to the model fitting function is handed over to
#' \code{circ.var} and its name (as present in the model to \code{to circ.var.name}.
#'
#' Caution, the function \code{boot.glmmtmb} is somewhat premature. For instance, it likely fails if the fitted
#' model doesn't contain random effects...
#'
#' @return Returns a list with the following objects:
#' \item{ci.estimates}{bootsrapped confidence intervals of model estimates (together with original estimate). In case of
#' the function \code{boot.lmer} this is a date frame, and in case of the function \code{boot.glmmtmb} this
#' is a named list comprising confidence intervals separately for fixed effects (\code{ci.estimates$fe}) and
#' random effects (\code{ci.estimates$re}).}
#' \item{ci.fitted}{confidence intervals of fitted model with regard to what was handed over to \code{use}. When \code{use} was
#' a single vector it will be a single data frame in case of the function \code{boot.lmer} with one column for
#' each term in the model (with all covariates not present in \code{use} set to zero), the fitted value, and
#' its confidence intervals (columns headed "lwr", and "upr"). In case of the function \code{boot.glmmtmb}, it
#' will be three such data frames, one for each of the conditional, the dispersion, and the zero-inflation model
#' (given these were present in the model). In case one handed over a list to \code{use}, \code{ci.fitted}
#' will be a list with one object as described so far for each entry in the list.}
#'
#' \item{warnings}{list with the warnings issued by the models fitted to the individual bootstraps (only
#' present when \code{boot.lmer} was used).}
#'
#' \item{all.boots}{Data frame with model estimates for each individual bootstrap.}
#'
#' @seealso
#' \code{\link[lme4]{bootMer}}, \code{\link[glmmTMB]{simulate.glmmTMB}}, \code{\link[kyotil]{keepWarnings}} (which I borrowed)
#' @references
#' Stolwijk, AM, Straatman H, & Zielhuis GA. 1999. Studying seasonality by using sine and cosine functions in regression analysis. J Epidemiol Community Health, 53:235–238.
#' @examples
#' ##examples for \code{boot.lmer}
#' \donttest{
#' ##minimal examlpe with 10 boostraps for testing:
#' library(lme4)
#' data(sleepstudy)
#' m = lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' boot.m=boot.lmer(m=m, nboots=10)
#' boot.m$ci.estimates
#' #obtain confidence intervals of fitted values (low resolution for testing):
#' boot.m=boot.lmer(m=m, nboots=10, use="Days", resol=10)
#' boot.m$ci.fitted
#' }
#' \donttest{
#' ##with factor in model:
#' data(Orthodont, package="nlme")
#' m=lmer(distance ~ age + Sex + (1|Subject), data=Orthodont)
#' boot.m=boot.lmer(m=m, nboots=10, use="age", resol=10)
#' boot.m$ci.fitted
#' #comprises a range of age values for each of females and males
#' boot.m=boot.lmer(m=m, nboots=10, use="Sex", resol=10)
#' boot.m$ci.fitted
#' ##age is set to zero in the data frame for which confidence intervals of fitted values are determined
#'
#' ##with a list handed over to argument use:
#' boot.m=boot.lmer(m=m, nboots=10, use=list("Sex", "age"), resol=10)
#' names(boot.m$ci.fitted)
#' #comprises two data frames; have a look:
#' boot.m$ci.fitted$Sex
#' boot.m$ci.fitted$age
#' }
#'
#' ##examples for boot.glmmtmb
#' \donttest{
#' library(glmmTMB)
#' m = glmmTMB(count ~ mined + (1|site), zi=~mined, family=poisson, data=Salamanders)
#' boot.m=boot.glmmtmb(m=m, nboots=10, resol=10, data=Salamanders, para=T)
#' boot.m$ci.estimates$fe
#' boot.m$ci.estimates$re
#'
#' ##bootstrap fitted values, too:
#' boot.m=boot.glmmtmb(m=m, nboots=10, resol=10, data=Salamanders, use="mined")
#' names(boot.m$ci.fitted)
#' boot.m$ci.fitted$cond
#' boot.m$ci.fitted$zi
#' boot.m$ci.fitted$disp##null because there is no dispersion part
#'
#' ##bootstrap fitted values for combinations of pedictors and zero-inflated part:
#' m = glmmTMB(count ~ spp + mined + DOP + (1|site), zi=~mined, family=poisson, data=Salamanders)
#' boot.m=boot.glmmtmb(m=m, nboots=10, resol=10, data=Salamanders, use=c("mined", "DOP"))
#' boot.m$ci.estimates$fe
#' head(boot.m$ci.fitted$cond)
#' head(boot.m$ci.fitted$zi)
#' boot.m$ci.fitted$disp##null because there is no dispersion part
#'
#' ##with use being a list:
#' boot.m=boot.glmmtmb(m=m, nboots=10, resol=10, data=Salamanders, use=list("mined", "DOP"))
#' names(boot.m$ci.fitted)
#' boot.m$ci.fitted$mined$cond
#' boot.m$ci.fitted$mined$zi
#' }
#'
#' @rdname boot.lmer
#' @export
boot.lmer<-function(m, discard.warnings=F, nboots=1000, para=F, resol=1000, level=0.95,
  use=NULL, circ.var.name=NULL, circ.var=NULL, use.u=F,
  n.cores=c("all-1", "all"), save.path=NULL, load.lib=F, lib.loc=.libPaths(), set.all.effects.2.zero=F){
  if(load.lib){library(lme4, lib.loc=lib.loc)}
  n.cores=n.cores[1]
  if(!is.null(use)){
    if(is.list(use)){use.list=use}else{use.list=list(use)}
  }else{
    use.list=NULL
  }
  #keepWarnings<-function(expr){##from package kyotil
  #  localWarnings <- list()
  #  value <- withCallingHandlers(expr,
	#															 warning = function(w) {
	#															   localWarnings[[length(localWarnings)+1]] <<- w
	#																 #invokeRestart(r="muffleWarnings")
	#															 }
  #  )
  #  list(value=value, warnings=localWarnings)
  #}
  ##define function extracting all estimated coefficients (fixed and random effects) and also the model summary wrt random effects:
  extract.all<-function(mres, use.u.=use.u){
    ##extract random effects model summary:
    vc.mat=as.data.frame(summary(mres)$varcor)##... and prepare/extract variance covariance matrix from the model handed over
    xx=lapply(summary(mres)$varcor, function(x){attr(x, "stddev")})##append residual variance
    ##create vector with names of the terms in the model
    xnames=c(names(fixef(mres)), paste(rep(names(xx), unlist(lapply(xx, length))), unlist(lapply(xx, names)), sep="@"))
    if(class(mres)[1]=="lmerMod"){xnames=c(xnames, "Residual")}##append "Residual" to xnames in case of Gaussian model
    if(class(mres)[1]=="lmerMod"){res.sd=vc.mat[vc.mat$grp=="Residual", "sdcor"]}##extract residual sd i case of Gaussian model
    #if(vc.mat$grp[nrow(vc.mat)]=="Residual"){vc.mat=vc.mat[-nrow(vc.mat), ]}##and drop residuals-row from vc.mat in case of Gaussisn model
    ##deal with names in vc.mat which aren't exactly the name of the resp. random effect
    r.icpt.names=names(ranef(mres))##extract names of random intercepts...
    not.r.icpt.names=setdiff(vc.mat$grp, r.icpt.names)##... and names of the random effects having random slopes
    not.r.icpt.names=unlist(lapply(strsplit(not.r.icpt.names, split=".", fixed=T), function(x){paste(x[1:(length(x)-1)], collapse=".")}))
    vc.mat$grp[!vc.mat$grp%in%r.icpt.names]=not.r.icpt.names
    if(vc.mat$grp[nrow(vc.mat)]=="Residual"){vc.mat$var1[nrow(vc.mat)]=""}
    xnames=paste(vc.mat$grp, vc.mat$var1, sep="@")
    re.summary=unlist(vc.mat$sdcor)
    names(re.summary)=xnames
    ranef(mres)
    ##extract random effects model summary: done
    ##extract random effects details:
    if(use.u.){
			re.detail=ranef(mres)
			xnames=paste(
				rep(x=names(re.detail), times=unlist(lapply(re.detail, function(x){nrow(x)*ncol(x)}))),
				unlist(lapply(re.detail, function(x){rep(colnames(x), each=nrow(x))})),
				unlist(lapply(re.detail, function(x){rep(rownames(x), times=ncol(x))})),
				sep="@")
			re.detail=unlist(lapply(re.detail, function(x){
				return(unlist(c(x)))
			}))
			names(re.detail)=xnames
		}else{
			xnames=NULL; re.detail=NULL
		}
    #deal with negative binomial model to extract theta:
    xx=as.character(summary(mres)$call)
    if(any(grepl(x=xx, pattern="negative.binomial"))){
      xx=xx[grepl(x=xx, pattern="negative.binomial")]
      xx=gsub(x=xx, pattern="MASS::negative.binomial(theta = ", replacement="", fixed=T)
      theta=as.numeric(gsub(x=xx, pattern=")", replacement="", fixed=T))
      #re.detail=c(re.detail, as.numeric(xx))
      #xnames=c(xnames, "theta")
    }else{
			theta=NULL
		}
    ns=c(length(fixef(mres)), length(re.summary), length(re.detail), length(theta))
    names(ns)=c("n.fixef", "n.re.summary", "n.blups", "n.theta")
    return(c(fixef(mres), re.summary, re.detail, theta=theta, ns))
  }
  if(discard.warnings){
    boot.fun<-function(x, m., keepWarnings., use.u., save.path.){
      xdone=F
      while(!xdone){
        i.res=keepWarnings.(bootMer(x=m., FUN=extract.all, nsim=1, use.u=use.u.)$t)
        if(length(unlist(i.res$warnings))==0){
          xdone=T
        }
      }
      est.effects=i.res
      i.warnings=NULL
      if(length(save.path.)>0){save(file=paste(c(save.path., "/b_", x, ".RData"), collapse=""), list=c("est.effects", "i.warnings"))}
      return(list(ests=i.res$value, warns=i.warnings))
    }
  }else{
    boot.fun<-function(y, m., keepWarnings., use.u., save.path.){
      #keepWarnings.(bootMer(x=m., FUN=fixef, nsim=1)$t)
      i.res=keepWarnings.(bootMer(x=m., FUN=extract.all, nsim=1, use.u=use.u.)$t)
      if(length(save.path.)>0){
        est.effects=i.res$t
        i.warnings=i.res$warnings
        save(file=paste(c(save.path., "/b_", y, ".RData"), collapse=""), list=c("est.effects", "i.warnings"))
      }
      return(list(ests=i.res$value, warns=unlist(i.res$warnings)))
    }
  }
  #browser()
  if(para){
    on.exit(expr = parLapply(cl=cl, X=1:length(cl), fun=function(x){rm(list=ls())}), add = FALSE)
    on.exit(expr = stopCluster(cl), add = T)
    library(parallel)
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    if(n.cores!="all"){
      if(n.cores!="all-1"){n.cores=length(cl)-1}
      if(n.cores<length(cl)){
        cl=cl[1:n.cores]
      }
    }
    parLapply(cl=cl, 1:length(cl), fun=function(x, .lib.loc=lib.loc){
      library(lme4, lib.loc=.lib.loc)
      return(invisible(""))
    })
    all.res=parLapply(cl=cl, X=1:nboots, fun=boot.fun, m.=m, keepWarnings.=keepWarnings, use.u.=use.u, save.path.=save.path)
  }else{
    all.res=lapply(X=1:nboots, FUN=boot.fun, m.=m, keepWarnings.=keepWarnings, use.u.=use.u, save.path.=save.path)
  }
  if(!discard.warnings){
    all.warns=lapply(all.res, function(x){
      xxx=unlist(x$warns)
      if(length(xxx)==0){xxx=""}
      return(xxx)
    })
  }else{
    all.warns=NULL
  }
	all.res=lapply(all.res, function(x){x$ests})
  if(!is.null(use.list)){
    #extract fixed effects terms from the model:
    xcall=as.character(m@call)[2]
    model.terms=attr(terms(as.formula(xcall)), "term.labels")
    REs=names(ranef(m))
    model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
    model=paste(model.terms, collapse="+")
    #build model wrt to the fixed effects:
    #exclude interactions and squared terms from model.terms:
    model.terms=unique(unlist(strsplit(x=model.terms, split=":", fixed=T)))
    model.terms=model.terms[!grepl(x=model.terms, pattern="I(", fixed=T)]
    model.terms=model.terms[!grepl(x=model.terms, pattern="^2)", fixed=T)]
    #create new data to be used to determine fitted values:
    ii.data=m@frame
    all.fitted=lapply(use.list, function(use){
      if(length(circ.var.name)==1){
        set.circ.var.to.zero=sum(circ.var.name%in%use)==0
      }else{
        set.circ.var.to.zero=F
      }

      if(length(use)==0){use=model.terms}
      new.data=vector("list", length(model.terms))
      usel=model.terms%in%use
      #if(length(use)>0)
      for(i in 1:length(model.terms)){
        if(is.factor(ii.data[, model.terms[i]])){
          new.data[[i]]=levels(ii.data[, model.terms[i]])
        }else if(!is.factor(ii.data[, model.terms[i]]) & usel[i] & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
          new.data[[i]]=seq(from=min(ii.data[, model.terms[i]]), to=max(ii.data[, model.terms[i]]), length.out=resol)
        }else  if(!is.factor(ii.data[, model.terms[i]]) & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
          new.data[[i]]=mean(ii.data[, model.terms[i]])
        }
      }
      names(new.data)=model.terms
      if(length(circ.var.name)==1){
        new.data=new.data[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
        if(sum(grepl(pattern=circ.var.name, x=use))>0){
          new.data=c(new.data, list(seq(min(circ.var, na.rm=T), max(circ.var, na.rm=T), length.out=resol)))
          names(new.data)[length(new.data)]=circ.var.name
        }else{
          new.data=c(new.data, list(0))
        }
        model.terms=model.terms[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
      }
      xnames=names(new.data)
      #browser()
      new.data=data.frame(expand.grid(new.data))
      names(new.data)=xnames
      if(length(circ.var.name)==1){
        names(new.data)[ncol(new.data)]=circ.var.name
      }
      if(set.all.effects.2.zero){
        for(iterm in setdiff(colnames(new.data), c("(Intercept)", use))){
          new.data[, iterm]=0
        }
      }
      m.mat=model.matrix(object=as.formula(paste(c("~", model), collapse="")), data=new.data)
      if(set.circ.var.to.zero){
        m.mat[,paste(c("sin(", circ.var.name, ")"), collapse="")]=0
        m.mat[,paste(c("cos(", circ.var.name, ")"), collapse="")]=0
      }
      #get CIs for fitted values:
      ci=lapply(all.res, function(x){
        #return(apply(m.mat[names(fixef(m)), , drop=F]*as.vector(x[, names(fixef(m))]), 2, sum))
        return(as.vector(m.mat[, names(fixef(m)), drop=F]%*%as.vector(x[, names(fixef(m))])))
      })
      ci=matrix(unlist(ci), ncol=nboots, byrow=F)
      ci=t(apply(ci, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
      colnames(ci)=c("lwr", "upr")
      fv=m.mat[, names(fixef(m))]%*%fixef(m)
      if(class(m)[[1]]!="lmerMod"){
        if(m@resp$family$family=="binomial"){
          ci=exp(ci)/(1+exp(ci))
          fv=exp(fv)/(1+exp(fv))
        }else if(m@resp$family$family=="poisson" | substr(x=m@resp$family$family, start=1, stop=17)=="Negative Binomial"){
          ci=exp(ci)
          fv=exp(fv)
        }
      }
      return(data.frame(new.data, fitted=fv, ci))
    })
    result=all.fitted
    if(length(use.list)==1){
			result=as.data.frame(result[[1]])
		}else{
			names(result)=unlist(lapply(use.list, paste, collapse="@"))
		}
  }else{
    result=NULL
  }
  all.boots=matrix(unlist(all.res), nrow=nboots, byrow=T)
  colnames(all.boots)=colnames(all.res[[1]])
  ci.est=apply(all.boots, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
  if(length(fixef(m))>1){
    ci.est=data.frame(orig=fixef(m), t(ci.est)[1:length(fixef(m)), ])
  }else{
    ci.est=data.frame(orig=fixef(m), t(t(ci.est)[1:length(fixef(m)), ]))
  }
  return(list(ci.estimates=ci.est, ci.fitted=result, warnings=all.warns, all.boots=all.boots))
}

#' @rdname boot.lmer
#' @export
boot.glmmtmb<-function(m, data, discard.non.conv=F, nboots=1000, para=F, resol=100, level=0.95,
                       use=NULL, circ.var.name=NULL, circ.var=NULL,
                       n.cores=c("all-1", "all"), save.path=NULL, load.lib=T, lib.loc=.libPaths(), set.all.effects.2.zero=F){
  if(!summary(m)$link%in%c("log", "logit", "identity")){
    stop("link functions other than log, logit, or identity aren't supported yet; please contact Roger")
  }
  if(load.lib){library(lme4, lib.loc=lib.loc)}
  if(!is.null(use)){
    if(is.list(use)){use.list=use}else{use.list=list(use)}
  }else{
    use.list=NULL
  }
  print("doesn't account for circular variables in the zero-inflated model")
  print("haven't tested it yet with a cbind response")
  n.cores=n.cores[1]
  #if(is.null(getCall(m)$control)){
  #  contr=glmmTMBControl()
  #}
  ##does the model comprise any random effects?
  has.RE=any(!unlist(lapply(summary(m)$varcor, is.null)))
  extract.BLUPs.from.glmmTB<-function(m){
    to.do=lapply(c("cond", "zi"), function(x){
      if(length(ranef(m)[[x]])>0){
        lapply(ranef(m)[[x]], function(y){
          return(list(
            what=rep(x, prod(dim(as.matrix(y)))),
            blup=unlist(c(as.matrix(y))),
            level=rep(rownames(y), times=ncol(y)),
            Name=rep(colnames(y), each=nrow(y))
          ))
        })
      }else{
        return(list(list(what=NULL, blup=NULL, level=NULL, Name=NULL)))
      }
    })
    xc=sapply(lapply(to.do[[1]], "[[", "what"), length)
    xz=sapply(lapply(to.do[[2]], "[[", "what"), length)
    res=data.frame(
      what=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "what")}))),
      Grp=c(rep(names(xc), times=xc), rep(names(xz), times=xz)),
      Name=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "Name")}))),
      level=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "level")}))),
      blup=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "blup")})))
    )
    res=res[order(res$what, res$Grp, res$Name, res$level), ]
    return(res)
  }
  extract.all<-function(mres){
    if(has.RE){
			##extract random effects model summary:
			m=extract.ranef.glmmTMB(m=mres)
			##extract random effects, BLUPs:
			BLUPs=extract.BLUPs.from.glmmTB(m=mres)
		}else{
			m=NULL
			BLUPs=NULL
		}
    ##dispersion parameter:
    dp=summary(mres)$sigma
    ##fixed effects
    xx=fixef(mres)
    fe=data.frame(
      what=rep(names(xx), times=sapply(xx, length)),
      term=unlist(lapply(xx, names)),
      est=unlist(xx)
    )
    fe=fe[order(fe$what, fe$term), ]
    return(list(fe=fe, sigma=dp, m=m, BLUPs=BLUPs))
  }
  ##prepare for new calls:
  xcall=as.character(m$call)
  names(xcall)=names(m$call)
  xcall=list(cond.form=xcall["formula"], zi.form=xcall["ziformula"], disp.form=xcall["dispformula"], xfam=xcall["family"])
  if(grepl(xcall[["xfam"]], pattern="(", fixed=T)){
    xfam=xcall[["xfam"]]
    xfam=unlist(strsplit(xfam, split="(", fixed=T))
    xfam[2]=gsub(x=xfam[2], pattern=")", replacement="", fixed=T)
    xfam[2]=gsub(x=xfam[2], pattern="link = ", replacement="", fixed=T)
    xfam[2]=gsub(x=xfam[2], pattern="\"", replacement="", fixed=T)
    if(substr(xcall[["xfam"]], start=1, stop=5)!="Gamma"){
      xcall[["xfam"]]=get(xfam[1])(xfam[2])
    }else{
      if(xfam[2]=="log"){
        xcall[["xfam"]]=Gamma(link="log")
      }else if(xfam[2]=="inverse"){
        xcall[["xfam"]]=Gamma(link="inverse")
      }else if(xfam[2]=="identity"){
        xcall[["xfam"]]=Gamma(link="identity")
      }else{
        stop("Error: family not supported")
				}
    }
  }

  #xweights=try(m$frame[, "(weights)"], silent=T)
  #if(class(xweights)[[1]]=="try-error"){
  #  xweights=rep(1, nrow(m$frame))
  #}
  #data$xweights.X.=xweights
  rv.name=gsub(unlist(strsplit(xcall[["cond.form"]], split="~"))[1], pattern=" ", replacement="", fixed=T)
  ##define function doing the bootstrap:
  boot.fun<-function(x, xcall., data., rv.name., m., discard.non.conv., save.path., extract.all.){
    xdone=F
    while(!xdone){
      done2=F
      while(!done2){
        data.[, rv.name.]=simulate(object=m.)[, 1]
        if(xcall[["xfam"]][["family"]]!="beta" |
           (xcall[["xfam"]][["family"]]=="beta" & min(data.[, rv.name.])>0 & max(data.[, rv.name.])<1)){
          done2=T
        }
      }
      i.res=try(update(m., data=data.), silent=T)
      #return(update(m., data=data.))
      #i.res=try(glmmTMB(formula=as.formula(xcall.[["cond.form"]]), ziformula=as.formula(xcall.[["zi.form"]]), dispformula=as.formula(xcall.[["disp.form"]]),
      #              data=data., control=contr., family=xcall[["xfam"]], weights=xweights.X.), silent=T)
      if(class(i.res)[[1]]!="try-error"){
				if(i.res$sdr$pdHess | !discard.non.conv.){
					xdone=T
				}
			}
    }
		if(class(i.res)[[1]]!="try-error"){
			conv=i.res$sdr$pdHess
			i.res=extract.all.(mres=i.res)
			if(length(save.path.)>0){save(file=paste(c(save.path., "/b_", x, ".RData"), collapse=""), list=c("i.res", "conv"))}
			return(list(boot.res=i.res, conv=conv))
		}else{
			return(NULL)
		}
  }
  ##do bootstrap:
  if(para){
    library(parallel)
    cl <- makeCluster(getOption("cl.cores", detectCores()))
    if(n.cores!="all"){
      if(n.cores=="all-1"){n.cores=length(cl)-1}
      if(n.cores<length(cl)){
        cl=cl[1:n.cores]
      }
    }
    clusterExport(cl=cl, varlist="extract.ranef.glmmTMB")
    parLapply(cl=cl, 1:length(cl), fun=function(x, .lib.loc=lib.loc){
      library(glmmTMB, lib.loc=.lib.loc)
      return(invisible(""))
    })
    all.res=parLapply(cl=cl, X=1:nboots, fun=boot.fun, m.=m, xcall.=xcall, data.=data, rv.name.=rv.name, 
                      discard.non.conv.=discard.non.conv, save.path.=save.path, extract.all.=extract.all)
  }else{
    all.res=lapply(X=1:nboots, FUN=boot.fun, m.=m, xcall.=xcall, data.=data, rv.name.=rv.name, 
                   discard.non.conv.=discard.non.conv, save.path.=save.path, extract.all.=extract.all)
  }
  
  ##extract results:
  all.res=all.res[!sapply(all.res, is.null)]
  ##convergence:
  all.conv=sapply(all.res, "[[", "conv")
  ##sigma:
  all.sigma=sapply(all.res, function(x){x$boot.res$sigma})
  ##estimates, fixef effects:
  all.fe=lapply(all.res, function(x){x$boot.res$fe$est})
  all.fe=matrix(unlist(all.fe), nrow=length(all.fe), byrow=T)
  colnames(all.fe)=apply(all.res[[1]]$boot.res$fe[, c("what", "term")], 1, paste, collapse="@")
  if(has.RE){
		##estimates, random effects:
		all.re=lapply(all.res, function(x){x$boot.res$m$sdcor})
		all.re=matrix(unlist(all.re), nrow=length(all.re), byrow=T)
		colnames(all.re)=apply(all.res[[1]]$boot.res$m[, c("part", "grp", "var1", "var2")], 1, paste, collapse="@")
		##estimates, BLUPs (well, it doesn't make any sense  to consider them):
		all.BLUPs=lapply(all.res, function(x){x$boot.res$BLUPs$blup})
		all.BLUPs=matrix(unlist(all.BLUPs), nrow=length(all.BLUPs), byrow=T)
		colnames(all.BLUPs)=apply(all.res[[1]]$boot.res$BLUPs[, c("what", "Grp", "Name", "level")], 1, paste, collapse="@")
	}else{
		all.re=NULL
		all.BLUPs=NULL
	}
  ##store results:
  all.ind.boots=list(fe=all.fe, re=all.re, sigma=all.sigma, conv=all.conv)
  ##get confidence intervals and merge them with original values:
  ##fixed effects:
  orig=extract.all(m)
  xx=orig$fe
  xx$comb=apply(xx[, c("what", "term")], 1, paste, collapse="@")
  #xx$comb=paste(rownames(xx), xx$comb, sep="@")
  ci=apply(all.fe, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
  #times=unlist(lapply(summary(m)$coefficients, nrow))
  #comb=paste(rep(c("cond", "zi", "disp")[1:length(times)], times=times), unlist(lapply(summary(m)$coefficients, rownames)), sep="@")
  ci=ci[, match(xx$comb, colnames(ci))]
  #xx=xx[match(comb, xx$comb), ]
  ci.fe=data.frame(orig=xx$est, t(ci))
  ##order according to original:
  ci.fe=ci.fe[paste(rep(c("cond", "zi", "disp"), times=unlist(lapply(fixef(m), length))), unlist(lapply(fixef(m), names)), sep="@"), ]
  ##random effects:
  if(has.RE){
		xx=orig$m
		xx$comb=apply(xx[, c("part", "grp", "var1", "var1")], 1, paste, collapse="@")
		ci=apply(all.re, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
		ci.re=data.frame(orig=xx$sdcor, t(ci))
	}else{
		ci.re=NULL
	}
  ##store results:
  ci.estimates=list(fe=ci.fe, re=ci.re)

  if(length(use.list)>0){
    ci.fitted.cond=NULL
    ci.fitted.zi=NULL
    ci.fitted.disp=NULL
    #extract fixed effects terms from the model:
    model.terms=c(attr(terms(as.formula(xcall[["cond.form"]])), "term.labels"),
                  attr(terms(as.formula(xcall[["zi.form"]])), "term.labels"),
                  attr(terms(as.formula(xcall[["disp.form"]])), "term.labels"))
    #exclude random effects:
    model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
    #exclude random effects, interactions and squared terms from model.terms:
    model.terms=model.terms[!grepl(x=model.terms, pattern="I(", fixed=T)]
    model.terms=model.terms[!grepl(x=model.terms, pattern=":", fixed=T)]
    model.terms=unique(model.terms)

    all.fitted=lapply(use.list, function(use){
      #create new data to be used to determine fitted values:
      if(length(use)==0){use=model.terms}
      new.data=vector("list", length(model.terms))
      if(length(circ.var.name)==1){
        set.circ.var.to.zero=sum(circ.var.name%in%use)==0
      }else{
        set.circ.var.to.zero=F
      }
      usel=model.terms%in%use
      #if(length(use)>0)
      for(i in 1:length(model.terms)){
        if(is.factor(data[, model.terms[i]])){
          new.data[[i]]=levels(data[, model.terms[i]])
        }else if(!is.factor(data[, model.terms[i]]) & usel[i] & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
          new.data[[i]]=seq(from=min(data[, model.terms[i]]), to=max(data[, model.terms[i]]), length.out=resol)
        }else  if(!is.factor(data[, model.terms[i]]) & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
          new.data[[i]]=mean(data[, model.terms[i]])
        }
      }
      names(new.data)=model.terms

      if(length(circ.var.name)==1){
        new.data=new.data[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
        if(sum(grepl(pattern=circ.var.name, x=use))>0){
          new.data=c(new.data, list(seq(min(circ.var, na.rm=T), max(circ.var, na.rm=T), length.out=resol)))
          names(new.data)[length(new.data)]=circ.var.name
        }else{
          new.data=c(new.data, list(0))
        }
        model.terms=model.terms[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
      }
      xnames=names(new.data)
      #browser()
      new.data=data.frame(expand.grid(new.data))
      names(new.data)=xnames
      #browser()
      if(length(circ.var.name)==1){
        names(new.data)[ncol(new.data)]=circ.var.name
      }
      if(set.all.effects.2.zero){
        for(iterm in setdiff(colnames(new.data), c("(Intercept)", use))){
          new.data[, iterm]=0
        }
      }
      ##prepare model frame for conditional model for prediction:
      model=attr(terms(as.formula(xcall[["cond.form"]])), "term.labels")
      model=model[!grepl(x=model, pattern="|", fixed=T)]
      if(length(model)==0){model="1"}
      cond.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
      if(set.circ.var.to.zero){
        cond.m.mat[,paste(c("sin(", circ.var.name, ")"), collapse="")]=0
        cond.m.mat[,paste(c("cos(", circ.var.name, ")"), collapse="")]=0
      }
      ##get bootstrapped fitted values for the conditional part:
      ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=5)=="cond@", drop=F]
      est.names=gsub(colnames(ests), pattern="cond@", replacement="")
      ci.cond=lapply(1:nrow(ests), function(x){
        return(cond.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
      })
      ci.cond=matrix(unlist(ci.cond), ncol=length(all.res), byrow=F)
      ci.cond=t(apply(ci.cond, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
      colnames(ci.cond)=c("lwr", "upr")
      ##get fitted value for the conditional part:
      ests=ci.estimates$fe
      ests=ests[substr(rownames(ests), start=1, stop=5)=="cond@", , drop=F]
      est.names=gsub(rownames(ests), pattern="cond@", replacement="")
      fv.cond=cond.m.mat[, est.names, drop=F]%*%ests$orig
      if(summary(m)$link=="log"){
        ci.fitted.cond=data.frame(fitted=exp(fv.cond), exp(ci.cond))
      }else if(summary(m)$link=="logit"){
        ci.fitted.cond=data.frame(fitted=plogis(fv.cond), plogis(ci.cond))
      }else if(summary(m)$link=="identity"){
        ci.fitted.cond=data.frame(fitted=fv.cond, ci.cond)
      }
      ci.fitted.cond=data.frame(new.data, ci.fitted.cond)

      ##prepare model frame for ZI model for prediction:
      if(!is.null(summary(m)$coefficients$zi)){
        model=attr(terms(as.formula(xcall[["zi.form"]])), "term.labels")
        model=model[!grepl(x=model, pattern="|", fixed=T)]
        if(length(model)==0){model="1"}
        zi.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
        ##get bootstrapped fitted values for the conditional part:
        ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=3)=="zi@", drop=F]
        est.names=gsub(colnames(ests), pattern="zi@", replacement="")
        ci.zi=lapply(1:nrow(ests), function(x){
          return(zi.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
        })
        ci.zi=matrix(unlist(ci.zi), ncol=length(all.res), byrow=F)
        ci.zi=t(apply(ci.zi, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
        colnames(ci.zi)=c("lwr", "upr")
        ##get fitted value for the conditional part:
        ests=ci.estimates$fe
        ests=ests[substr(rownames(ests), start=1, stop=3)=="zi@", , drop=F]
        est.names=gsub(rownames(ests), pattern="zi@", replacement="")
        fv.zi=zi.m.mat[, est.names, drop=F]%*%ests$orig
        ci.fitted.zi=data.frame(fitted=plogis(fv.zi), plogis(ci.zi))
        #ci.fitted.zi=ci.fitted*(1-ci.fitted.zi)
        ci.fitted.zi=data.frame(new.data, ci.fitted.zi)
      }

      ##prepare model frame for disp model for prediction:
      if(!is.null(summary(m)$coefficients$disp)){
        model=attr(terms(as.formula(xcall[["disp.form"]])), "term.labels")
        model=model[!grepl(x=model, pattern="|", fixed=T)]
        if(length(model)==0){model="1"}
        disp.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
        ##get bootstrapped fitted values for the conditional part:
        ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=5)=="disp@", drop=F]
        est.names=gsub(colnames(ests), pattern="disp@", replacement="")
        ci.disp=lapply(1:nrow(ests), function(x){
          return(disp.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
        })
        ci.disp=matrix(unlist(ci.disp), ncol=length(all.res), byrow=F)
        ci.disp=t(apply(ci.disp, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
        colnames(ci.disp)=c("lwr", "upr")
        ##get fitted value for the dispersion part:
        ests=ci.estimates$fe
        ests=ests[substr(rownames(ests), start=1, stop=5)=="disp@", , drop=F]
        est.names=gsub(rownames(ests), pattern="disp@", replacement="")
        fv.disp=disp.m.mat[, est.names, drop=F]%*%ests$orig
        ci.fitted.disp=data.frame(fitted=exp(fv.disp), exp(ci.disp))
        #ci.fitted.disp=ci.fitted*(1-ci.fitted.disp)
        ci.fitted.disp=data.frame(new.data, ci.fitted.disp)
      }
      ci.fitted=list(
        cond=ci.fitted.cond,
        zi=ci.fitted.zi,
        disp=ci.fitted.disp
      )
      return(ci.fitted)
    })
    if(length(use.list)>1){
      names(all.fitted)=unlist(lapply(use.list, paste, collapse="@"))
    }else{
      all.fitted=all.fitted[[1]]
    }
  }else{
    all.fitted=NULL
  }
  return(list(ci.estimates=ci.estimates, ci.fitted=all.fitted, all.boots=all.ind.boots))
}
