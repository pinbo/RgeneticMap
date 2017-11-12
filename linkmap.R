linkmap <- function(object, chr, chr.space = 2, m.cex = 0.6, interval = FALSE, ruler = FALSE, ...){
  # VERSION: 1.1.0
  # object: a "cross" object from R/qtl, or a "map" class from the output of "pull.map" in R/qtl, or a data frame with marker column, chromosme column and position column named as "mar", "chr" and "pos", respectively
  # chr: a vector of chromosome names that need to be drawn.
  # chr.space: space between each chromosomes
  # m.cex: font size
  # interval: NULL/TRUE/FALSE: plot no distance/marker interval/absolute distance. Default is absolute distance.
  # ruler: whether to draw a common left ruler on the left
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
  thelim <- range(chrpos) + c(-1.0, 1.0)
  
  # function to get actual plotting postion
  finalPos = function(pos, maxlen) {#pos is a vector of positions of A genetic map; maxlen is the length of chromosome segments
      posnew = pos
      if (length(posnew) > 1){
          conv <- par("pin")[2]/maxlen # pin is The current plot dimensions, (width, height), in inches.
          for(j in 1:(length(pos) - 1)){
              ch <- posnew[j + 1]*conv - (posnew[j]*conv + 10*par("csi")*m.cex/9) # csi is the default font height
              #cat(ch)
              if (ch < 0){
                  temp <- posnew[j]*conv + 10*par("csi")*m.cex/9
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
  total_len = maxlen - minlen # in case chromosome did not start from 0
  mt = lapply(map, finalPos, maxlen= total_len) # final label postion on the right
  if (!is.null(interval)) {
      if (interval) {
          map2 = lapply(map, function(x) getInterval(x)[[1]]) # chromsome position for the left
          map3 = lapply(map, function(x) getInterval(x)[[2]]) # intervals
          mt2 = lapply(map2, finalPos, maxlen=total_len) # distance plotting position on the left
      } else {
          map2 = map # left and right are the same
          map3 = map
          mt2 = mt # left and right are the same
      }
  }
  
  maxlen <- max(c(unlist(lapply(omap, max)),unlist(lapply(mt, max))))
  names(mt) <- names(map)
  
  par(mar=c(0.6 ,1.1 ,2.6 ,1.1))
  if (ruler) par(mar=c(0.6, 2.1, 2.6, 1.1))

  plot(0, 0, type = "n", ylim = c(maxlen, minlen), xlim = thelim,
       xaxs = "i", ylab = "", xlab = "", axes = FALSE, ...)
  if (ruler) axis(side = 2,  ylim = c(maxlen, minlen))
  else chrpos = chrpos - 0.3
  
  ## Draw chromosomes
  barwidth = 0.14
  seglen = 0.06

  for(i in 1:n.chr) {
    # for the right side plotting
    Rstart = chrpos[i] + barwidth # start point for legs on the RIGHT
    segments(Rstart, map[[i]], Rstart + seglen, map[[i]])
    segments(Rstart + seglen, map[[i]], Rstart + seglen*3, mt[[i]])
    segments(Rstart + seglen*3, mt[[i]], Rstart + seglen*4, mt[[i]])
    alis <- list(x = Rstart + seglen*4 + 0.05, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
    do.call("text", c(alis, dots)) # draw right labels
	# JZ: add distance on the left
    if (!is.null(interval)){
      Lstart = chrpos[i] - barwidth # start point for legs on the LEFT
	    segments(Lstart, map2[[i]], Lstart - seglen, map2[[i]])
      segments(Lstart - seglen, map2[[i]], Lstart - seglen*3, mt2[[i]])
      segments(Lstart - seglen*3, mt2[[i]], Lstart - seglen*4, mt2[[i]])
      alisL <- list(x = Lstart - seglen*4 - 0.05, y =  mt2[[i]], labels = format(round(map3[[i]], 1),nsmall=1),adj = c(1, 0.5), cex = m.cex)
      do.call("text", c(alisL, dots)) # draw left labels
    }
    # draw chromosome bar
    map[[i]] <- omap[[i]]
    barl <- chrpos[i] - barwidth/2
    barr <- chrpos[i] + barwidth/2
    segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 3)
    segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 3)
    segments(barl - barwidth/2, map[[i]], barr + barwidth/2, map[[i]]) # bar ribs
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
  axis(side = 3, at = chrpos, labels = names(map), tick = F, cex.axis=1.5) # side = 1 for bottom, side=3 for top
  #if(is.na(pmatch("main", names(dots))) & !as.logical(sys.parent()))
  #  title("Genetic Map")
  invisible(list(mt = mt, map = map, chrpos = chrpos))
}
