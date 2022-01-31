##function plotting a quantitative response against a covariate and adding the model line
##written by Roger Mundry
no.int.plot<-function(
	plot.data, covariate, coefs=NULL, response, link=c("identity", "logit", "log"),
  grid.resol=11, 
	xlab="", ylab="", x.at=NULL, x.labels=NULL, y.at=NULL, y.labels=NULL, xlim=NULL, ylim=NULL,
	cex.lab=1, cex.axis=1, size.fac=NA, col="black", pch=1, conf.int=NULL, ci.type=c("lines", "area"), lty=2, lwd=1, pt.lwd=1,
	my.par=list(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1), quiet=T, reset.par=F, mean.fun=mean,
	weights=NULL){
	#last modified: 2018 AUG 20, works also with model without intercept
  old.par = par(no.readonly = TRUE)
  #version from april 15th 2016
	if(length(ci.type)>1){
		ci.type=ci.type[1]
	}
	#browser()
  link=link[1]
	#extract coefficients needed:
	if(!is.null(coefs)){
		coefs=coefs[names(coefs)%in%c("(Intercept)", covariate, paste(c("I(", covariate, "^2)"), collapse=""))]
		#coefs=coefs[c(which(names(coefs)=="(Intercept)"), grep(x=names(coefs), pattern=covariate))]
		#coefs=coefs[!grepl(x=names(coefs), pattern=":", fixed=T)]
	}
	if(length(dim(response))==2){
		response=response[, 1]/apply(response, 1, sum)
	}
  if(is.null(weights)){
		weights=rep(1, nrow(plot.data))
	}
	sel=!is.na(weights) & !is.na(response) & !is.na(plot.data[, covariate])
	weights=weights[sel]
	response=response[sel]
	plot.data=subset(plot.data, sel)
  if(!is.na(grid.resol)){
		xvar=seq(min(plot.data[, covariate]), max(plot.data[, covariate]), length.out=grid.resol)
		bin.x=cut(x=plot.data[, covariate[1]], breaks=xvar, labels=F, include.lowest=T)
		bin.x=min(xvar)+diff(xvar)[1]/2+((bin.x-1)*diff(xvar)[1])
  }else{
    bin.x=plot.data[, covariate]
  }

  if(!is.na(size.fac) & !is.na(grid.resol)){
		obs.mean=aggregate(response, list(bin.x), mean.fun)
    #get sample size per 'cell':
    N=aggregate(weights, by=list(bin.x), sum)
    #... and turn the result into a data frame:
		N=N$x
  }else if(!is.na(size.fac) & is.na(grid.resol)){
		obs.mean=aggregate(response, list(bin.x, response), mean.fun)
    N=aggregate(weights, by=list(bin.x, response), sum)
		N=N$x
  }else{
		obs.mean=aggregate(response, list(bin.x), mean.fun)
    N=rep(1, nrow(obs.mean))
    size.fac=1
  }
  names(obs.mean)[1]=covariate
	names(obs.mean)[2]="Freq"
  obs.mean=data.frame(obs.mean, N)
  obs.mean=subset(obs.mean, N>0)
	#get predicted values:
	xvals=seq(min(plot.data[, covariate]), max(plot.data[, covariate]), length.out=100)
  if(any(names(coefs)=="(Intercept)")){
		x.pred.data=data.frame(1, xvals)
		names(x.pred.data)=c("(Intercept)", covariate)
	}else{
		x.pred.data=data.frame(xvals)
		names(x.pred.data)=covariate
	}
	
	if(!is.null(coefs)){
		model.terms=names(coefs)
		if(length(intersect(model.terms, "(Intercept)"))>0){model.terms[model.terms=="(Intercept)"]="1"}
		rr=runif(n=nrow(x.pred.data))
		pvs=model.matrix(object=as.formula(paste(c("rr", paste(model.terms, collapse="+")), collapse="~")), data=x.pred.data)
		yvals=apply(t(coefs*t(pvs[, names(coefs)])), 1, sum)
		#transform according to link function
		if(link=="log"){yvals=exp(yvals)}
		if(link=="logit"){yvals=exp(yvals)/(1+exp(yvals))}
	}

	#prepare axis ticks and labels
  if(length(xlim)==0){
    xlim=range(plot.data[,covariate], na.rm=T)
  }
  if(length(ylim)==0){
		if(!is.null(coefs)){
			ylim=range(c(response, yvals))
		}else{
			ylim=range(response)
		}
	}
  if(length(x.at)==0){
    x.at=pretty(plot.data[, covariate])
    #x.at=x.at[x.at>=min(plot.data[, covariate])]
    #x.at=x.at[x.at<=max(plot.data[, covariate])]
    x.labels=x.at
  }
  if(length(y.at)==0){
    y.at=pretty(ylim)
    #y.at=y.at[y.at>=min(ylim)]
    #y.at=y.at[y.at<=max(ylim)]
    y.labels=as.character(y.at)
  }

  #plot:
  par(my.par)
	plot(obs.mean[, covariate], obs.mean$Freq, cex=size.fac*sqrt(obs.mean$N), type="n", 
		xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, tcl=-0.25, cex.lab=cex.lab)
	#add axes:
	axis(at=x.at, labels=x.labels, tcl=-0.2, side=1, tcl=-0.25, cex.axis=cex.axis)
	#mtext(text=xaxt.labels, at=xaxt.at, side=1, line=0.3, cex=0.7)
	axis(at=y.at, tcl=-0.2, side=2, las=1, labels=y.labels, tcl=-0.25, cex.axis=cex.axis)
	#add CI:
	if(length(conf.int)>0){
		if(ci.type=="lines"){
			ci.yvals=conf.int[, "lwr"]
			lines(x=conf.int[, covariate], y=conf.int[, "lwr"], lty=3, lwd=lwd)
			lines(x=conf.int[, covariate], y=conf.int[, "upr"], lty=3, lwd=lwd)
		}else if(ci.type=="area"){
			polygon(x=c(conf.int[, covariate], rev(conf.int[, covariate])), 
				y=c(conf.int[, "lwr"], rev(conf.int[, "upr"])), border=NA, col="grey")
		}
	}
	points(obs.mean[, covariate], obs.mean$Freq, cex=size.fac*sqrt(obs.mean$N), pch=pch, col=col, lwd=pt.lwd)
	if(!is.null(coefs)){
		lines(xvals, yvals, lty=lty, lwd=lwd)
	}
	if(reset.par){par(old.par)}
	if(!quiet){return(obs.mean)}
}
