#!/usr/bin/env Rscript
library("rgl")
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}


stuff <- read.table('test.txt',header=FALSE,colClasses=c("numeric","numeric","numeric","factor"))
x <- stuff[[1]]
y <- stuff[[2]]
z <- stuff[[3]]
fcolors <- get_colors(stuff[[4]])
rgl.open()
rgl.points(x,y,z,color=fcolors)
