linkmap <- function(object, chr, chr.space = 2, m.cex = 0.6, interval = FALSE, ...){
  # VERSION: 1.0.1
  # object: a "cross" object from R/qtl, or a "map" class from the output of "pull.map" in R/qtl, or a data frame with marker column, chromosme column and position column named as "mar", "chr" and "pos", respectively
  # chr: a vector of chromosome names that need to be drawn.
  # chr.space: space between each chromosomes
  # m.cex: font size
  # interval: NULL/TRUE/FALSE: plot no distance/marker interval/absolute distance. Default is absolute distance.
  # ...: other plot parameters
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
  #mt <- list()

  maxlen <- max(unlist(lapply(map, max)))
  minlen <- min(unlist(lapply(map, min)))
  omap <- map # omap is original map
  if(!is.na(pmatch("cex", names(dots))))
      dots$cex <- NULL
  #    else cex <- par("cex")
  chrpos <- seq(1, n.chr * chr.space, by = chr.space)
  thelim <- range(chrpos) + c(-1.3, 1.3)
  
  # function to get actual plotting postion
  finalPos = function(pos, maxlen) {#pos is a vector of positions of A genetic map
      posnew = pos
      if (length(posnew) > 1){
          conv <- par("pin")[2]/maxlen
          for(j in 1:(length(pos) - 1)){
              ch <- posnew[j + 1]*conv - (posnew[j]*conv + 10*par("csi")*m.cex/9)
              if(ch < 0){
                  temp <- posnew[j + 1]*conv + abs(ch)
                  posnew[j + 1] <- temp/conv
              }
          }
      }
      return(posnew)
  }
  
  # function to get interval
  getInterval = function(pos){#pos is a vector of positions of A genetic map
      pos2 = unique(round(pos,1))
      ll = length(pos2)
      pos3 = (pos2[1:(ll-1)] + pos2[2:ll])/2 # map position
      pos4 = pos2[2:ll] - pos2[1:(ll-1)] # intervals
      return(list(pos3,pos4))
  }
  
  # plot left side
  mt = lapply(map, finalPos, maxlen=maxlen) # final label postion on the right
  if (!is.null(interval)) {
      if (interval) {
          map2 = lapply(map, function(x) getInterval(x)[[1]]) # chromsome position for the left
          map3 = lapply(map, function(x) getInterval(x)[[2]]) # intervals
          mt2 = lapply(map2, finalPos, maxlen=maxlen) # distance plotting position on the left
      } else {
          map2 = map # left and right are the same
          map3 = map
          mt2 = mt # left and right are the same
      }
  }
  
  maxlen <- max(c(unlist(lapply(omap, max)),unlist(lapply(mt, max))))
  names(mt) <- names(map)

  plot(0, 0, type = "n", ylim = c(maxlen, minlen), xlim = thelim,
       xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome",
       axes = FALSE, ...)
  axis(side = 2,  ylim = c(maxlen, minlen))

  for(i in 1:n.chr) {
    # for the right side plotting
    alis <- list(x = chrpos[i] + 0.50, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
    do.call("text", c(alis, dots))
    segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 0.3, map[[i]])
    segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, mt[[i]])
    segments(chrpos[i] + 0.40, mt[[i]], chrpos[i] + 0.45, mt[[i]])
	# JZ: add distance on the left
    if (!is.null(interval)){
      alisL <- list(x = chrpos[i] - 0.50, y =  mt2[[i]], labels = format(round(map3[[i]], 1),nsmall=1),adj = c(1, 0.5), cex = m.cex)
	  do.call("text", c(alisL, dots))
      segments(chrpos[i] - 0.25, map2[[i]], chrpos[i] - 0.3, map2[[i]])
      segments(chrpos[i] - 0.3, map2[[i]], chrpos[i] - 0.4, mt2[[i]])
      segments(chrpos[i] - 0.40, mt2[[i]], chrpos[i] - 0.45, mt2[[i]])
    }
    # draw chromosome bar
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
