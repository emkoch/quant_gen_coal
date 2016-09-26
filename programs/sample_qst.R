## FUNCTION TO SAMPLE A QST VALUE FROM ITS DISTRIBUTION
## n is a vector of sample sizes and tau is a matrix of expected coal times
sqst <- function(n, tau, m=1){
    k <- length(n)
    sigma <- c()
    s1 <- 0; s2 <- 0; s3 <- 0
    for (ii in 1:k) {
        sigma[ii] <- ( 2 / sum(n) ) * ( (n %*% tau[ii,]) - tau[ii,ii] )
        s1 <- ( 2 / sum(n) ) * ( (n %*% tau[ii,]) - tau[ii,ii] )
        sigma[ii] <- sigma[ii] - tau[ii,ii] * (n[ii] - 1)/n[ii]
        s2 <- tau[ii,ii] * (n[ii] - 1)/n[ii]
        tmp <- 0
        for (jj in 1:k) {
            tmp <- tmp + n[jj] * ( (n %*% tau[jj,]) - tau[jj,jj] )
        }
        s3 <- ( 1 / sum(n)^2 ) * tmp
        sigma[ii] <- sigma[ii] - ( 1 / sum(n)^2 ) * tmp
    }
    result <- list()
    result$V.between <- c(0)
    result$V.within <- c(0)
    result$Qst <- c(0)
    for (mm in 1:m) {
        if (mm %% 100 == 0){
            print( mm)
        }
        V.between <- 0
        for (ii in 1:k) {
            V.between <- V.between + sigma[ii]*rchisq(1,1)#*n[ii]
        }
        V.between <- V.between / (k-1)
        V.within <- 0
        for (ii in 1:k) {
            V.within <- V.within + (n[ii] - 1)/n[ii] * tau[ii,ii] * rchisq(1,df=n[ii]-1)
            ## for (jj in 1:n[ii]) {
            ##     V.within <- V.within + (n[ii] - 1)/n[ii] * tau[ii,ii] * rchisq(1,1)
            ## }
        }
        V.within <- V.within / ( sum(n) - k )
        result$V.between[mm] <- V.between
        result$V.within[mm] <- V.within
        
        result$Qst[mm] <- V.between / (V.within + V.between) 
    }
    return( result )
}

test.matrix <- function(within, between, n=10){
    result <- matrix(data=0, nrow=n, ncol=n)
    for (ii in 1:n) {
        for (jj in 1:n) {
            if (ii == jj) {
                result[ii,jj] <- within
            } else{
                result[ii,jj] <- between
            }
        }
    }
    return( result )
}

ring.matrix <- function(N, m, d){
    result <- matrix(data=0, nrow=d, ncol=d)
    for (ii in 1:d) {
        for (jj in 1:d) {
            b <- abs(ii - jj)
            result[ii,jj] <- 2*N*d + (d - b) * b/(2 * m)
        }
    }
    return(result)
}

doop <- c()
for (i in 1:5000) {
    doop[i] <- sqst( n=rep(20,20), test.matrix(within=1, between=1.3, n=20 ))
}    

my.mat <- ring.matrix(100, 0.05, 10)
my.mat <- test.matrix(1.0, 1.1, 10)
Nreps <- 100000
poop.10 <- sqst( n=rep(10,10), my.mat, m=Nreps)
poop.20 <- sqst( n=rep(20,10), my.mat, m=Nreps)
poop.50 <- sqst( n=rep(50,10), my.mat, m=Nreps)
poop.100 <- sqst( n=rep(100,10), my.mat, m=Nreps)
poop.200 <- sqst( n=rep(200,10), my.mat, m=Nreps)
poop.300 <- sqst( n=rep(300,10), my.mat, m=Nreps)
poop.1000 <- sqst( n=rep(1000,10), my.mat, m=Nreps)

colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")

pdf("Qst_test.pdf")
plot(density(poop.1000$Qst), main="", xlab="Qst", lty=4, col=colors[4], lwd=2.5, xlim=c(0,.4), cex.lab=1.6, cex.axis=1.4)
lines(density(poop.10$Qst), lty=1, col=colors[1], lwd=2.5)
#lines(density(poop.20$Qst))
lines(density(poop.50$Qst), lty=2, col=colors[2], lwd=2.5)
#lines(density(poop.100$Qst))
lines(density(poop.200$Qst), lty=3, col=colors[3], lwd=2.5)
legend("topright", legend=c(10,50,200,1000), cex=1.3, col=colors, lty=1:4, lwd=rep(2,4))
dev.off()
#lines(density(poop.300$Qst))
abline(v=(1.1-1)/1.1)
abline(v=mean(poop.1000$Qst))

plot(density(poop.1000$V.between))
lines(density(poop.10$V.between), lty=2, col=colors[1], lwd=2.5)
#lines(density(poop.20$V.between))
lines(density(poop.50$V.between), lty=3, col=colors[2], lwd=2.5)
#lines(density(poop.100$V.between))
lines(density(poop.200$V.between), lty=4, col=colors[3], lwd=2.5)


plot(density(poop.10$V.within))
lines(density(poop.1000$V.within), lty=2, col=colors[1], lwd=2.5)
#lines(density(poop.20$V.within))
lines(density(poop.50$V.within), lty=3, col=colors[2], lwd=2.5)
#lines(density(poop.100$V.within))
lines(density(poop.200$V.within), lty=4, col=colors[3], lwd=2.5)
