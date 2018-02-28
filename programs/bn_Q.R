library(ggplot2)
library(viridis)
library(directlabels)

ET2 <- function(bb, cc){
  return(1 - exp(-bb) + cc*exp(-bb))
}

ET3 <- function(bb, cc){
  return( (1/6)*exp(-3*bb)*(1 - cc) - (3/2)*exp(-bb)*(1 - cc) + 4/3 )
}

ET4 <- function(bb, cc){
  return( (-1/30)*exp(-6*bb) + (1/3)*exp(-3*bb) - (9/5)*exp(-bb) + 3/2 + 
            cc*((1/30)*exp(-6*bb) - (1/3)*exp(-3*bb) + (9/5)*exp(-bb)) )
}

Qfactor <- function(bb, cc){
  return( (4*ET2(bb, cc)  - 6*ET3(bb, cc) + 3*ET4(bb, cc))/ET2(bb, cc) )
}

Qfactor.alt <- function(bb, cc){
  return(0.1*(5 - 4*exp(-bb)*(1 - cc) - exp(-6*bb)*(1-cc) )/(1 - exp(-bb) + cc*exp(-bb)))
}

test.vals <- expand.grid(bb = 10^seq(-2, 1, length.out = 2e2),
                         cc = 10^seq(-2, 2, length.out = 2e2))
N0 <- 1e4
theta <- 2e-8
L <- 1e4
test.vals$Q <- apply(test.vals, 1, FUN=function(xx) Qfactor.alt(xx[1], xx[2]))

pdf("Q_land.pdf", height=3.5, width=7*(3.5/5))
foo <- ggplot(data=test.vals, aes(x=bb, y=cc, z=2*Q)) + 
  geom_raster(data=test.vals, aes(fill=2*Q), show.legend=TRUE) + 
  scale_x_log10(breaks=c(.01,.1,1,10,100)) + scale_y_log10(breaks=c(.01,.1,1,10,100)) + 
  geom_contour(aes(color=..level..), show.legend = FALSE, breaks=c(0.85, 0.95, 1.0, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55))+
  scale_fill_viridis(option="B", direction = -1) + 
  labs(x="Time to pop size change", y="Pop size change", fill="Q") + theme_bw()
foo2 <- direct.label.ggplot(foo, method = list("far.from.others.borders"))
foo2
dev.off()
