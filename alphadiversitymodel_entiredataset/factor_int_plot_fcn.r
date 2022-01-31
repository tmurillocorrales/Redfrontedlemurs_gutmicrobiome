#This function was written by Roger Mundry for plotting variables that are factors in linear models.

factor.int.plot<-function(
  plot.data, factors, coefs=NULL, response, link=c("identity", "logit", "log"),
  conf.int=NULL, yaxt.labels=NULL, yaxt.at=NULL, ylab, ylim=NULL, size.fac=1,
	factor.seqs=NULL, factor.labels=NULL, to.show=c("bubbles", "data", "bars", "boxes"), weights=NULL, which.q=c(5, 3),
	pch=1, pt.col=par("fg"), rect.col=NULL, border.col=NULL, est.ci.col=par("fg"), est.ci.lty=1, quiet=T, average.response=F, log.y=F, 
	bg="white", median.col="black", median.lwd=2, fitted.lwd=2, quant.col=NULL,
	cex.lab=1, cex.lab.x=1, cex.axis=1, reset.par=T, percentile.col=NULL, percentile.lwd=NULL, my.par=NULL, cex.axis.lab=1, xlim=NULL,
	add.range=F, range.pch=4){
	#last updated: 2015 Nov 26
	#arguments:
		#plot.data: data frame comprising all variables involved
		#factors: character, vector with the names of the factors to be depicted (the hierarchical order in which they are displayed is the reverse of factors, i.e.,
			#the first factor will be the upper most at the x-axis and the last the lowermost)
		#coefs: named vector with the model coefficients to be shown (default NULL in which case no model coefficients are shown)
		#response: vector with the values of the response
		#link: character, one of "identity" (the default), "logit", or "log"; potentially leads to the model as inferred from coefs to be transformed
		#conf.int: numeric data frame or matrix with columns to be named "lwr" and "upr", respectively, and also columns for the resp. factors 
			#(precisely, the dummmy coded levels, named as in model output) indicating the CIs; default is NULL in which no CIs are depicted
		#yaxt.labels: numeric, labels to be displayed along the y-axis (default, NULL in which labels are determined from the response)
		#yaxt.at: numeric, values along the y-axis where ticks are to be placed (default, NULL in the values labels are determined from the response)NULL
		#ylab: as  usual
		#ylim: numeric, as usual (default: NULL, in which case ylim is the range of tthe response)
		#size.fac: numeric, determines the size of the points depicting the data in case to.show comprises "bubbles"
			#(default 1; if points are too large make it smaller, otherwise larger)
		#factor.seqs: character list, determines the order with which the levels of the factors are depicted along the x-axis
			#(default is NULL in which it corresponds to their levels)
		#factor.labels: character list, determines the labels which are displayed for the levels of the factors along the x-axis
			#(default is NULL in which it corresponds to their levels)
		#to.show: character vector, can be one of "bubbles", "bars", or "boxes" or a combination of he former with one of the latter two;
			#determines whether the data are shown (to.show comprises "bubbles") and whether a summary in the form of bars depicting medians ((to.show comprises "bars") or
			#boxes etc. (to.show comprises "boxes") are depicted (see also argument 'which.q')
		#which.q: numeric, one of 5 or 3; if 5 (the default) medians, quartiles and quantiles (2.5 and 97.5%) are shown, if 3 medians and quartiles are shown
		#pch, numeric, one of 1:25 (default is 1); determines the 'point character' (has the same interpretation as usual)
		#pt.col: character, color with which points/bubbles are depicted (default: par("fg"))
		#rect.col: character, color with which bars or boxes are depicted (default: NA  in which case they are shown with a border in the default foreground color and no filling)
		#est.ci.col: character, color with which CIs are depicted (default: par("fg"))
		#quiet: logical, determines whether the function returns the number of observations per bubble (only relevant if to.show comprieses "bubbles"
	old.par = par(no.readonly = TRUE)
	#if response is cbind:
	if(length(dim(response))==2){
		response=response[, 1]/apply(response, 1, sum)
	}
	which.q=which.q[1]
	#get rid of NAs:
	plot.data=droplevels(as.data.frame(na.omit(data.frame(plot.data[, factors, drop=F], response))))
	link=link[1]#get default link (in case not specified)
	#create lists factor.seqs and factor.seqs in case they are empty:
	if(length(factor.seqs)==0){
		factor.seqs=lapply(factors, function(x){levels(plot.data[, x])})
		names(factor.seqs)=factors
	}else{
		lapply(1:length(factors), function(x){
			if(length(levels(plot.data[, factors[x]]))!=length(factor.seqs[[x]])){
				stop("error: order of factors doesn't match in factors and factor.seqs")
			}else{
				if(sum(sort(levels(plot.data[, factors[x]]))!=sort(factor.seqs[[x]]))>0){
					stop("error: factor levels don't match in factors and factor.seqs")
				}
			}
		})
	}
	names(factor.seqs)=factors
	if(length(factor.labels)==0){factor.labels=factor.seqs}
	#set ylim in case it is empty:
	if(length(ylim)==0){
		ylim=range(plot.data[, "response"])
	}
	if(length(yaxt.labels)==0){
		yaxt.labels=pretty(ylim)
		yaxt.at=yaxt.labels
	}
	#create data matrix (needed for something taking place later (guess it has to do with the order of the factor levels)):
	pred.data=data.frame(expand.grid(lapply(factors, function(x){levels(plot.data[, x])})))
	names(pred.data)=factors
	pred.data$all=apply(pred.data[, factors, drop=F], 1, paste, collapse="@@@")
	pred.data=pred.data[match(apply(data.frame(expand.grid(factor.seqs)), 1, paste, collapse="@@@"), pred.data$all), ]
	##create corresponding matrix for x-axis labels
	labels.mat=data.frame(expand.grid(factor.labels))
	###subset such that only level combinations existing in the data are considered:
	to.keep=unique(apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@"))
	labels.mat=subset(labels.mat, pred.data$all%in%to.keep)
	pred.data=subset(pred.data, all%in%to.keep)
	if(length(coefs)>0){
		#extract coefficients needed:
		#paste factors and their levels together:
		to.search=lapply(factors, function(x){levels(plot.data[, x])[-1]})
		to.search=paste(rep(factors, times=unlist(lapply(to.search, length))), unlist(to.search), sep="")
		tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
			length(cf)==sum(unlist(lapply(to.search, function(v){
				cf==v
			})))
		}))
		if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
		coefs=coefs[tk]
		#determine fitted values:
		#browser()
		m.mat=model.matrix(as.formula(paste(c("~", paste(factors, collapse="*")), collapse="")), data=pred.data)
		#get fitted values:
		colnames(m.mat)=unlist(lapply(strsplit(colnames(m.mat), split=":", fixed=T), function(x){paste(sort(x), collapse=":")}))
		names(coefs)=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(x){paste(sort(x), collapse=":")}))
		xfitted=apply(t(t(m.mat[, names(coefs)])*coefs), 1, sum)
		if(link=="log"){
			xfitted=exp(xfitted)
		}else if(link=="logit"){
			xfitted=exp(xfitted)/(1+exp(xfitted))
		}
		###########################################################################################
		##probably needs to be adjusted to order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		###########################################################################################
	}
	if(log.y & min(plot.data[, "response"])>0){
		plot.data[, "response"]=log(plot.data[, "response"])
	}else if(log.y & min(plot.data[, "response"])==0){
		plot.data[, "response"]=log(plot.data[, "response"]+1)
	}
	if(sum("bubbles"%in%to.show)==1 & length(ylim)==0){
		ylim=range(plot.data[, "response"])
	}
	if(sum(c("bars", "boxes")%in%to.show)==1){
		plot.data$all=apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@")
		q.mat=tapply(plot.data[, "response"], plot.data$all, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
		xnames=names(q.mat)
		q.mat=as.data.frame(matrix(unlist(q.mat), nrow=length(q.mat), byrow=T))
		colnames(q.mat)=c("x0", "x025", "x25", "x5", "x75", "x975", "x1")
		rownames(q.mat)=xnames
		q.mat$all=rownames(q.mat)
		q.mat=q.mat[match(pred.data$all, q.mat$all), ]
	}
	
	hbw=0.25
	u.all.fac=pred.data[, "all"]
	#begin with plot:
  if(is.null(my.par)){
	  par(mar=c(length(factors)+0.2*length(factors)+0.6, 3.4, 0.2, 0.2), mgp=c(1.8, 0.35, 0), tcl=-0.2, bg=bg, mex=cex.lab*0.7, las=1)
  }else{
    par(my.par)
  }
	if(sum("bubbles"%in%to.show)==1 & !average.response){
		obs.mean=as.data.frame(table(plot.data[, c("response", factors)]))
		obs.mean$response=as.numeric(as.character(obs.mean$response))
		all.fac=apply(obs.mean[, factors, drop=F], 1, paste, collapse="@@@")
		obs.mean$all.fac=all.fac
		obs.mean=subset(obs.mean, Freq>0)
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
	}else if(sum("bubbles"%in%to.show)==1 & average.response){
		obs.mean=aggregate(plot.data$response, plot.data[, factors], mean)
		names(obs.mean)[ncol(obs.mean)]="response"
		obs.mean$Freq=aggregate(plot.data$response, plot.data[, factors], length)$x
		all.fac=apply(obs.mean[, factors, drop=F], 1, paste, collapse="@@@")
		obs.mean$all.fac=all.fac
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
	}else if(sum("data"%in%to.show)==1){
		plot.data$all.fac=apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@")
		if(is.null(weights)){weights=rep(1, nrow(plot.data))}
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(plot.data$all.fac, u.all.fac), y=plot.data$response, cex=size.fac*sqrt(weights), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
		obs.mean=data.frame(plot.data, Freq=weights)
	}else if(sum("bars"%in%to.show)==1){
		if(is.null(xlim)){xlim=xlim}
		plot(x=1, y=1, xlim=c(0.5, nrow(pred.data)+0.5), xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, type="n", yaxs="i", 
			yaxt="n", cex.lab=cex.axis.lab)
	}else{
		if(is.null(xlim)){xlim=c(0.5, nrow(pred.data)+0.5)}
		plot(x=1, y=1, xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, type="n", yaxt="n", 
      cex.lab=cex.axis.lab)
	}
	axis(side=2, las=1, at=yaxt.at, labels=yaxt.labels, cex.axis=cex.axis, cex.lab=cex.lab, las=my.par[["las"]])
	#add points in case specified;
	if(sum("bubbles"%in%to.show)==1){
		points(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), 
			pch=pch, col=pt.col)
	}else if(sum("data"%in%to.show)==1){
		points(x=match(plot.data$all.fac, u.all.fac), y=plot.data$response, cex=size.fac*sqrt(weights), pch=pch, col=pt.col)
	}
	#add bars or boxes etc. in case specified;
	if(is.null(rect.col)){rect.col=par("bg")}
	if(is.null(border.col)){border.col=par("fg")}
	if(length(quant.col)==0){quant.col=median.col}
	if(sum("bars"%in%to.show)==1){
		if(average.response){
			obs.mean.2=aggregate(plot.data$response, plot.data[, factors], mean)
			obs.mean.2$all=apply(obs.mean.2[, factors, drop=F], 1, paste, collapse="@@@")
			obs.mean.2=obs.mean.2[match(obs.mean.2$all, pred.data$all), "x"]
		}else{
			obs.mean.2=q.mat[match(q.mat$all, pred.data$all), "x5"]
		}
		rect(xleft=(1:nrow(pred.data))-hbw, xright=(1:nrow(pred.data))+hbw, ybottom=min(ylim), ytop=obs.mean.2, border=border.col, col=rect.col)
	}else if(sum("boxes"%in%to.show)==1){
		rect(xleft=(1:nrow(pred.data))-hbw, xright=(1:nrow(pred.data))+hbw, ybottom=q.mat[match(q.mat$all, pred.data$all), "x25"], ytop=q.mat[match(q.mat$all, pred.data$all), "x75"], col=rect.col, border=border.col)
		segments(x0=(1:nrow(pred.data))-hbw, x1=(1:nrow(pred.data))+hbw, y0=q.mat[match(q.mat$all, pred.data$all), "x5"], y1=q.mat[match(q.mat$all, pred.data$all), "x5"], lwd=median.lwd, col=median.col, lend=1)
		if(which.q==5){
			segments(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=q.mat[match(q.mat$all, pred.data$all), "x75"], y1=q.mat[match(q.mat$all, pred.data$all), "x975"], col=quant.col)
			segments(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=q.mat[match(q.mat$all, pred.data$all), "x25"], y1=q.mat[match(q.mat$all, pred.data$all), "x025"], col=quant.col)
		}
		if(add.range){
			points((1:nrow(pred.data)), y=q.mat[match(q.mat$all, pred.data$all), "x0"], pch=range.pch)
			points((1:nrow(pred.data)), y=q.mat[match(q.mat$all, pred.data$all), "x1"], pch=range.pch)
		}
	}
	#add text at x-axis:
	#u.all.fac=u.all.fac[match(u.all.fac, pred.data$all)]
	#u.all.fac.mat=matrix(unlist(strsplit(u.all.fac, split="@@@", fixed=T)), ncol=length(factors), byrow=T)
	#u.all.fac.mat=u.all.fac.mat[, ncol(u.all.fac.mat):1, drop=F]
	mtext(text=labels.mat[, 1], side=1, line=1-1+0.3*1, 
		at=1:nrow(labels.mat), cex=cex.lab.x)
	if(ncol(labels.mat)>1){
		for(i in 2:ncol(labels.mat)){
			to.add.1=c(T, labels.mat[-nrow(labels.mat), i]!=labels.mat[-1, i])
			to.add.2=cumsum(to.add.1)
			to.add.2=tapply(1:nrow(labels.mat), to.add.2, mean)
			mtext(text=labels.mat[to.add.1, i], side=1, line=i-1+0.3*i, 
				at=to.add.2,
				#at=sort(tapply(1:nrow(u.all.fac.mat), apply(u.all.fac.mat[, 1:(length(factors)+1-i), drop=F], 1, paste, collapse="@@@"), mean)),
				cex=cex.lab.x)
		}
	}
	#add model coefficients:
	if(length(coefs)>0){
		if(log.y & min(plot.data[, "response"])>0){
			xfitted=log(xfitted)
		}else if(log.y & min(plot.data[, "response"])==0){
			xfitted=log(xfitted+1)
		}
		hll=0.15
		segments(x0=(1:nrow(pred.data))-hll, x1=(1:nrow(pred.data))+hll, y0=xfitted, y1=xfitted, lwd=fitted.lwd, col=est.ci.col, lty=est.ci.lty, lend=1)
	}
	#add CIs:
	if(length(conf.int)>0){
		conf.int$all=apply(conf.int[, factors, drop=F], 1, paste, collapse="@@@")
		conf.int=conf.int[match(u.all.fac, conf.int$all), ]
		arrows(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=conf.int[, "lwr"], y1=conf.int[, "upr"], lwd=fitted.lwd, length=0.05, code=3, angle=90, col=est.ci.col, lty=est.ci.lty)
	}
  if(reset.par){par(old.par)}
	#if data are to be shown and function shouldn't be quiet
	if((sum(to.show%in%"bubbles")==1 | sum(to.show%in%"data")==1) & !quiet){
		obs.mean=obs.mean[, names(obs.mean)!="all.fac"]
		return(obs.mean)#return object telling sample size per bubble
	}
}
