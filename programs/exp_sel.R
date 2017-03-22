d.sel.1 <- function(G, s, V){
    return(exp(s*G) * exp(s^2*V/2))
}

d.sel.2 <- function(G, s, k, V){
    p1 <- exp((2*G*(G*k+s)+s^2*V) / (2-4*k*V))
    return(sqrt(1-2*k*V)^(-1)*p1)
}

G.vals <- seq(.01, 3, by=.01)
tt.1 <- unlist(lapply(G.vals, d.sel.1, s=.1, V=1))
tt.2 <- unlist(lapply(G.vals, d.sel.2, s=.1, V=1, k=0.01))
plot(G.vals, tt.2, type="l")
points(G.vals, tt.1, type="l")

m.1 <- function(Z, s, k, V){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    wg <- d.sel.2(Z, s, k, V)
    result <- wg
    return(result)
}

m.2 <- function(Z, s, k, V, M3){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    wg <- d.sel.2(Z, s, k, V)
    result <- wg * (1 + (V/2)*(B + A^2))
    return(result)
}

m.3 <- function(Z, s, k, V, M3){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    wg <- d.sel.2(Z, s, k, V)
    result <- wg * (1 + (V/2)*(B + A^2) + (M3/6) * (3*A*B + A^3))
    return(result)
}

m.4 <- function(Z, s, k, V, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    wg <- d.sel.2(Z, s, k, V)
    result <- wg * (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
                        (M4/24)*(3*B^2 + 6*A^2*B + A^4))
    return(result)
}

m.1(1, .3, .02, 1)
m.2(1, .3, .02, 1)
m.3(1, .3, .02, 1, 0)
m.4(1, .3, .02, 1, 0, 3)


L1 <- function(Z, s, k, V, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- A + (V/2)*(3*A*B + A^3) + (M3/6)*(3*B^2 + 6*A^2*B + A^4) +
        (M4/24)*(15*B^2*A + 10*A^3*B + A^5)
    return(result/D)
}

L2 <- function(Z, s, k, V, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- B + A^2
    return(result/(2*D))
}

L3 <- function(Z, s, k, V, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- 3*A*B + A^3
    return(result/(6*D))
}

L4 <- function(Z, s, k, V, M3, M4){
    A <- (2*G*k + s) / (1 - 2*k*V)
    B <- (2*k) / (1 - 2*k*V)
    D <- (1 + (V/2)*(B + A^2) + (M3/6)*(3*A*B + A^3) + 
              (M4/24)*(3*B^2 + 6*A*B + A^4))
    result <- 3*B^2 + 6*A^2*B + A^4
    return(result/(24*D))
}

s <- .1; k <- .1

L1(1, s, k, 1, 0, 0)
l.1 <- L1(1, s, k, 1, 0, 3)

L2(1, s, k, 1, 0, 0)
l.2 <- L2(1, s, k, 1, 0, 3)

L3(1, s, k, 1, 0, 0)
l.3 <- L3(1, s, k, 1, 0, 4)

L4(1, s, k, 1, 0, 0)
l.4 <- L4(1, s, k, 1, 0, 3)

plot(c(l.1,l.2,l.3,l.4), log="y")
