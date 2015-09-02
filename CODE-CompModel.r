## Example code to run compositional model in Allen et al. (submitted)
## check http://github.com/tjleininger/comp-data/downloads for latest version

require(MASS)
require(tmvtnorm)
require(coda)
require(bayesm)
require(msm)
require(fields)
require(matlab)

## INPUT
# Y <- (n x D) matrix containing n D-variate observations, the last column is used as a baseline (no zeros allowed)
# X <- (n x p) matrix or dataframe with desired covariates
# prox <- (n x n) neighborhood matrix; (i,j) entry = 1 indicates location i and j are neighbors, 0 otherwise
# burn <- desired number of burn-in MCMC iterations
# Length <- desired number of MCMC

n <- nrow(Y)  # number of observations
D <- ncol(Y)
d <- D - 1
pred.Y <- array(0, c(Length, n, D))
p <- ncol(X)  # number of covariates
kthin <- 1  # thinning parameter
niter <- burn + Length * kthin  # total number of MCMC iterations
zero.ind <- Y[, 1:d] == 0
n.zero <- sum(sum(zero.ind))
which.zero <- which(rowSums(zero.ind[, 1:d]) > 0)
if (n.zero != 1) full.zero <- rowSums(zero.ind[which.zero, ]) == d
if (n.zero == 1) full.zero <- sum(zero.ind[which.zero, ]) == d
neg.ind <- c(t(zero.ind[which.zero, ]))
diag(prox) <- 0
n.i <- rowSums(prox)
max.neigh <- max(n.i)
if (sum(n.i == 0) > 0) stop("Remove observations with no neighbors.")
sX <- cbind(1, scale(X, scale = T))  # scale X to improve model fitting
XpX <- t(sX) %*% sX
XpXiXp <- solve(XpX) %*% t(sX)


## Priors
lambda <-10^3 ## prior variance for beta_kl
m.V <- D+1  ## for V = cov mat for Z
W.V <- diag(1,D-1)  ## for V 
m.sig <- D+1      ## for Sigma = cov of phi_i
W.sig <- diag(1,D-1)   ## for Sigma

# Initialize parameters
curr.Z <- sqrt(Y[, 1:d]/Y[, D])
curr.Z[curr.Z == 0] <- -mean(curr.Z)
curr.beta <- matrix(0, p + 1, d)
betas <- array(0, c(Length, p + 1, d))
curr.V <- W.V/m.V
curr.Vi <- solve(curr.V)
Vs <- array(0, c(Length, d, d))
curr.phi <- matrix(0, n, d)
phis <- array(0, c(Length, n, d))
curr.sig <- W.sig/m.sig
curr.sigi <- solve(curr.sig)
sigs <- array(0, c(Length, d, d))
DwmW <- diag(n.i) - prox
Lts.phi <- array(0, c(max.neigh, d, d))

##################
### BEGIN MCMC ###
##################
cat("Starting MCMC\n", date(), "\n")
for (mcmcit in 2:niter) {
  
  ## update beta
  temp <- curr.Vi %x% XpX
  Lt <- chol(temp + diag(1/lambda, d * (p + 1)))
  b <- temp %*% c(XpXiXp %*% (curr.Z - curr.phi))
  curr.beta <- matrix(solve(Lt, solve(t(Lt), b)) + solve(Lt, rnorm(d * (p + 1))), p + 1, d)
  curr.xb <- sX %*% curr.beta
  
  ## update V
  temp <- curr.Z - curr.xb - curr.phi
  curr.V <- rwishart(m.V + n, solve(W.V + t(temp) %*% temp))$IW
  curr.Vi <- solve(curr.V)
  
  ## update Z iteratively (to truncate)
  if (n.zero > 0) {
    for (i in 1:length(which.zero)) {
      ii <- which.zero[i]
      if (full.zero[i]) { ## all zeros
        curr.Z[ii, ] <- rtmvnorm(1, c(curr.xb[ii, ] + curr.phi[ii, ]), curr.V, upper = rep(0, d), algorithm = "gibbs", 
                                 burn.in.samples = 100)
      } else { ## some zeros; need conditional dist.
        ind1 <- which(zero.ind[ii, ])
        ind2 <- which(!zero.ind[ii, ])
        umean <- c(curr.xb[ii, ] + curr.phi[ii, ])
        temp <- curr.V[ind1, ind2] %*% solve(curr.V[ind2, ind2])
        cvar <- curr.V[ind1, ind1] - temp %*% curr.V[ind2, ind1]
        cmean <- umean[ind1] - temp %*% (umean[ind2] - c(curr.Z[ii, ind2]))
        if (length(cmean) > 1) {
          curr.Z[ii, ind1] <- rtmvnorm(1, c(cmean), cvar, upper = rep(0, length(ind1)), algorithm = "gibbs", burn.in.samples = 100)
        } else {
          curr.Z[ii, ind1] <- rtnorm(1, cmean, sqrt(cvar), upper = 0)
        }
      }
    }
  }
  
  ## update phi
  for (temp in 1:max.neigh) {
    Lts.phi[temp, , ] <- chol(temp * curr.sigi + curr.Vi)
  }
  for (i in 1:n) {
    Lt <- Lts.phi[n.i[i], , ]
    temp <- curr.Vi %*% (curr.Z[i, ] - curr.xb[i, ]) + curr.sigi %*% t(prox[i, ] %*% curr.phi)
    curr.phi[i, ] <- solve(Lt, solve(t(Lt), temp)) + solve(Lt, rnorm(d))
  }
  curr.phi <- curr.phi - matrix(colMeans(curr.phi), n, d, byrow = T)
  
  ## update sigma
  curr.sig <- rwishart(m.sig + n - 1, solve(W.sig + t(curr.phi) %*% DwmW %*% curr.phi))$IW
  curr.sigi <- solve(curr.sig)
  
  if (mcmcit > burn) {
    if ((mcmcit - burn)%%kthin == 0) {
      ## store draws
      kp <- (mcmcit - burn)/kthin
      betas[kp, , ] <- curr.beta
      Vs[kp, , ] <- curr.V
      phis[kp, , ] <- curr.phi
      sigs[kp, , ] <- curr.sig
      
      ## get predictive draws
      pred <- curr.xb + curr.phi + matrix(rnorm(n * d), n, d) %*% chol(curr.V)
      pred[pred < 0] <- 0
      pred <- pred^2
      pred.Y[kp, , ] <- cbind(pred, 1)/(1 + rowSums(pred))
    }
  }
  if (mcmcit%%(niter/10) == 0) 
    cat(100 * mcmcit/niter, "% done \n", sep = "")
}
cat("End of MCMC \n", date(), "\n")


##########################################################
### Plot posterior means of each component E(y_k|data) ###
### while varying one covariate                        ### 
##########################################################

## Input
lseq <- 20 ## number of evaluations per covariate
Ylabels <- paste('Y',1:D,sep='')
Xlabels <- paste('X',1:p,sep='')

## Calculate posterior means
xrange<-matrix(0,2,p)
xrange[1,]<-apply(X,2,min)
xrange[2,]<-apply(X,2,max)
xrange.std<-matrix(0,2,p)
xrange.std[1,]<-apply(sX[,-1],2,min)
xrange.std[2,]<-apply(sX[,-1],2,max)
Eyk <- array(0,c(lseq,D,p))
pY <- matrix(0,Length,d)
cVs <- array(0,c(niter,d,d))
for(iter in 1:Length){
  cVs[iter,,]<-chol(Vs[iter,,])
}
for(j in 1:p){
  xseq<-seq(xrange.std[1,j],xrange.std[2,j],l=lseq)
  for(l in 1:lseq){
    x <- xseq[l]
    for(iter in 1:Length){
      pY[iter,]<-betas[iter,1,]+x*betas[iter,j+1,]+rnorm(d)%*%cVs[iter,,]
    }
    pY[pY<0] <- 0
    pY2 <- cbind(pY^2,1)
    pY2 <- pY2/rowSums(pY2)
    Eyk[l,,j]<-colMeans(pY2)
  }
}
par(mfrow=c(1,3),mar=c(3,4,2,2))
plot(1,1,col='white',pch='.',axes=F,xlab='',ylab='')
legend(1,1,Ylabels,col=1:D, xjust=0.5, yjust=0.5,cex=1.1,lwd=1)
for(j in 1:(dim(betas)[2]-1)){
  xseq<-seq(xrange[1,j],xrange[2,j],l=dim(Eyk)[1])
  matplot(xseq,Eyk[,,j],type='l', main=Xlabels[j],ylab='',col=1:D, xlab='',ylim=c(0,1),cex.main=1.5,cex=0.7,lwd=1,lty=1)
  mtext(expression(paste('E(', y[k],'|x,data)')),side=2,line=2,cex=.75)
}

