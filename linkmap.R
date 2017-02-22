linkmap <- function(object, chr, chr.space = 2, m.cex = 0.6, ...){
  # map should be a list, the same format as  map <- pull.map(object)
  if ("data.frame" %in% class(object)){ # transform to a list
	  pos = object$pos
	  names(pos) = object$mar
	  map = split(pos,object$chr)
  } else if ("cross" %in% class(object)){
  	map <- pull.map(object)
  } else map = object # a map object or list
  dots <- list(...) # extra parameters inputed by user
  old.xpd <- par("xpd")
  par(xpd = TRUE)
  on.exit(par(xpd = old.xpd))
  
  if(!missing(chr)) {
    if(any(is.na(pmatch(chr, names(map)))))
      stop("Some names of chromosome(s) subset do not match names of map.")
    map <- map[chr]
  }
  n.chr <- length(map)
  mt <- list()

  maxlen <- max(unlist(lapply(map, max)))
  minlen <- min(unlist(lapply(map, min)))
  omap <- map
  if(!is.na(pmatch("cex", names(dots))))
      dots$cex <- NULL
  #    else cex <- par("cex")
  chrpos <- seq(1, n.chr * chr.space, by = chr.space)
  thelim <- range(chrpos) + c(-1.3, 1.3)

  for(i in 1:n.chr){
      mt[[i]] <- map[[i]]
      if(length(mt[[i]]) > 1){
          conv <- par("pin")[2]/maxlen
          for(j in 1:(length(mt[[i]]) - 1)){
              ch <- mt[[i]][j + 1]*conv - (mt[[i]][j]*conv + 10*par("csi")*m.cex/9)
              if(ch < 0){
                  temp <- mt[[i]][j + 1]*conv + abs(ch)
                  mt[[i]][j + 1] <- temp/conv
              }
          }
      }
  }
  
  maxlen <- max(c(unlist(lapply(omap, max)),unlist(lapply(mt, max))))
  names(mt) <- names(map)

  plot(0, 0, type = "n", ylim = c(maxlen, minlen), xlim = thelim,
       xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome",
       axes = FALSE, ...)
  axis(side = 2,  ylim = c(maxlen, minlen))

  for(i in 1:n.chr) {
    alis <- list(x = chrpos[i] + 0.50, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
    do.call("text", c(alis, dots))
    segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 0.3, map[[i]])
    segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, mt[[i]])
    segments(chrpos[i] + 0.40, mt[[i]], chrpos[i] + 0.45, mt[[i]])
	# JZ: add distance on the left
    alisL <- list(x = chrpos[i] - 0.50, y =  mt[[i]], labels = format(round(map[[i]], 1),nsmall=1),adj = c(1, 0.5), cex = m.cex)
	do.call("text", c(alisL, dots))
    segments(chrpos[i] - 0.25, map[[i]], chrpos[i] - 0.3, map[[i]])
    segments(chrpos[i] - 0.3, map[[i]], chrpos[i] - 0.4, mt[[i]])
    segments(chrpos[i] - 0.40, mt[[i]], chrpos[i] - 0.45, mt[[i]])
    
    map[[i]] <- omap[[i]]
    barl <- chrpos[i] - 0.07
    barr <- chrpos[i] + 0.07
    segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 3)
    segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 3)
    segments(barl - 0.13, map[[i]], barr + 0.13, map[[i]]) #0.2 from chrpos on each side
    # attempt to put curves at ends of chromosomes
    rs <- seq(0,pi,len=100)
	r <- (barr - barl)/2 # radius
	xunit = par("pin")[1]/abs(par("xaxp")[2] - par("xaxp")[1])
	yunit = par("pin")[2]/abs(par("yaxp")[2] - par("yaxp")[1])
    xseq <- r*cos(rs) 
    yseq <- r*sin(rs)*(xunit/yunit)
    lines(xseq + chrpos[i], min(map[[i]]) - yseq, lwd=3)
    lines(xseq + chrpos[i], max(map[[i]]) + yseq, lwd=3)
  }
  axis(side = 1, at = chrpos, labels = names(map), tick = F)
  if(is.na(pmatch("main", names(dots))) & !as.logical(sys.parent()))
    title("Genetic Map")
  invisible(list(mt = mt, map = map, chrpos = chrpos))
}
