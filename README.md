# RgeneticMap
A R function to draw genetic maps (linkage map), and can be used with R/qtl directly.  
You can export the maps as svg format and edit them with any svg editor (Inkscape is a good FREE svg editor).

# Example

```r
library(qtl)
data("hyper")
summary(hyper)

source("linkmap.R")
linkmap(hyper,chr=c(1,2,3))

pdf(file="geneticMap.pdf",paper="USr")
linkmap(hyper,chr=c(1,2))
linkmap(hyper,chr=c(1,2,3))
linkmap(hyper,chr=c(1,2,3,4))
dev.off()

# Absolute positions
svg(file="GeneticMap.svg", width=7, height=7, bg=NA)
linkmap(hyper,chr=c(1,2,3), ruler = T, m.cex = 0.8)
dev.off()

# Relative positions
svg(file="GeneticMap-interval.svg", width=7, height=7, bg=NA)
linkmap(hyper,chr=c(1,2,3), ruler = T, m.cex = 0.8, interval=T)
dev.off()
```
# SVG format output example

### Absolute positions
![map](GeneticMap.svg)


### Relative positions
![map2](GeneticMap-interval.svg)

# Parameters
**linkmap <- function(object, chr, chr.space = 2, m.cex = 0.6, interval = FALSE, ruler = FALSE, ...){...}**

Input of this function:
- **object**:
  + a "cross" object from R/qtl
  + a "map" class from the output of "pull.map" in R/qtl
  + a data frame with marker column, chromosme column and position column named as "**mar**", "**chr**" and "**pos**", respectively
- **chr**: a vector of chromosome names that need to be drawn.
- **chr.space**: space between each chromosomes
- **m.cex**: font size
- **interval**: NULL/TRUE/FALSE: plot no distance/marker interval/absolute distance. Default is absolute distance.
- **ruler**: whether to draw ruler on the left
- **...**: other plot parameters

# Credits
This script is a modification of function "link.map.cross" in R package "wgaim" (https://cran.r-project.org/web/packages/wgaim/index.html). I added marker positions on the left side of the chromosomes and make it look better.

# Other genetic map packages in R
- **LinkageMapView: Plot Linkage Group Maps with Quantitative Trait Loci** (https://cran.r-project.org/web/packages/LinkageMapView/index.html) [Right now it can only export pdf format]
