d.sel.1 <- function(G, s, Ve){
    return(exp(s*G) * exp(s^2*Ve/2))
}

d.sel.2 <- function(G, s, k, Ve){
    p1 <- exp((2*G*(G*k+s)+s^2*Ve) / (2-4*k*Ve))
    return(sqrt(1-2*k*Ve)^(-1)*p1)
}

G.vals <- seq(.01, 3, by=.01)
tt.1 <- unlist(lapply(G.vals, d.sel.1, s=.3, Ve=1))
tt.2 <- unlist(lapply(G.vals, d.sel.2, s=.3, Ve=1, k=0.1))
plot(G.vals, tt.2, type="l")
points(G.vals, tt.1, type="l")

m.1 <- function(Z, s, k, V, Ve){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    wg <- d.sel.2(Z, s, k, Ve)
    result <- wg
    return(result)
}

m.2 <- function(Z, s, k, V, Ve, M3){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    wg <- d.sel.2(Z, s, k, Ve)
    result <- wg * (1 + (V/2)*(B + A^2))
    return(result)
}

m.3 <- function(Z, s, k, V, Ve, M3){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    wg <- d.sel.2(Z, s, k, Ve)
    result <- wg * (1 + (V/2)*(B + A^2) + (M3/6) * (3*A*B + A^3))
    return(result)
}

m.4 <- function(Z, s, k, V, Ve, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    wg <- d.sel.2(Z, s, k, Ve)
    result <- wg * (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
                        (M4/24)*(3*B^2 + 6*A^2*B + A^4))
    return(result)
}

m.1(1, .3, .02, 1)
m.2(1, .3, .02, 1)
m.3(1, .3, .02, 1, 0)
m.4(1, .3, .02, 1, 0, 3)

L1 <- function(Z, s, k, V, Ve, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- A + (V/2)*(3*A*B + A^3) + (M3/6)*(3*B^2 + 6*A^2*B + A^4) +
        (M4/24)*(15*B^2*A + 10*A^3*B + A^5)
    return(result/D)
}

L2 <- function(Z, s, k, V, Ve, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- B + A^2
    return(result/(2*D))
}

L3 <- function(Z, s, k, V, Ve, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- 3*A*B + A^3
    return(result/(6*D))
}

L4 <- function(Z, s, k, V, Ve, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*Ve)
    B <- (2*k) / (1 - 2*k*Ve)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- 3*B^2 + 6*A^2*B + A^4
    return(result/(24*D))
}

delta.z <- function(Z, s, k, V, Ve, M3, M4){
    l1 <- L1(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=M4)
    l2 <- L2(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=M4)
    l3 <- L3(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=M4)
    result <- V*l1 + M3*l2 + (M4 - 3*V)*l3
    return(result)
}

pdf("sel_resp.pdf", width=9, height=3)
par(mfrow = c(1,3))
Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .001; k <- .01
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)

plot((m4.set - 3*V)/(V^2), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess kurtosis", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.01")))
points((m4.set - 1.5*V)/((0.5*V)^2), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V)/(4*V^2), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))

Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .001; k <- .05
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)
plot((m4.set - 3*V)/(V^2), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess kurtosis", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.05")))
points((m4.set - 1.5*V)/((0.5*V)^2), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V)/(4*V^2), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))

Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .001; k <- .1
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)
plot((m4.set - 3*V)/(V^2), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess kurtosis", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.1")))
points((m4.set - 1.5*V)/((0.5*V)^2), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V)/(4*V^2), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))
dev.off()



pdf("sel_resp_mom4.pdf", width=9, height=3)
par(mfrow = c(1,3))
Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .3; k <- .01
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)

plot((m4.set - 3*V)/(V^2), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess fourth moment", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.01")))
points((m4.set - 1.5*V), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))

Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .3; k <- .05
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)
plot((m4.set - 3*V), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess fourth moment", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.05")))
points((m4.set - 1.5*V), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))

Z <- 1; V <- 1; Ve <- 1; M3 <- 0; M4 <- 4; s <- .3; k <- .1
d.test.1 <- mapply(delta.z, Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=V, Ve=Ve, M3=M3, M4=3*V)
d.test.2 <- mapply(delta.z, Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=0.5*V, Ve=0.5*Ve, M3=M3, M4=3*0.5*V)
d.test.3 <- mapply(delta.z, Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=m4.set) /
    delta.z(Z=Z, s=s, k=k, V=2*V, Ve=2*Ve, M3=M3, M4=3*2*V)
plot((m4.set - 3*V), d.test.1, type="l", 
     xlim=c(-1.5,4), ylim=c(0.8,1.4), log="y", xlab="excess fourth moment", 
     ylab="Response to selection relative to normal", lwd=2, lty=2, 
     main = expression(paste(epsilon, "= 0.1")))
points((m4.set - 1.5*V), d.test.2, type="l", lwd=2, lty=1)
points((m4.set - 6*V), d.test.3, type="l", lwd=2, lty=3)
legend("topleft", lwd=rep(2,3), lty=1:3, legend=paste("CV = ", c(0.7, 1, 1.4)))
dev.off()

