# Example compositional data to be used with CODE-CompModel.r
# Does fairly well given the small sample size

library(MASS)

set.seed(123)
n <- 25
burn <- 1000
Length <- 5000
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

## Fix parameters for simulating
V <- matrix(c(1, 0.25, 0.25, 1), 2, 2)
SIGMA <- matrix(c(1, -0.25, -0.25, 1), 2, 2)
BETA <- 1 * matrix(c(1, -1, 0.5, 0.5, 0, 0.5), 3, 2)
PHI <- matrix(0, n, 2)

## Generate phis
prox <- 1 * (as.matrix(dist(expand.grid(1:5, 1:5), diag = T, upper = T)) < 1.5)
diag(prox) <- 0
n.i <- rowSums(prox)
for (iter in 1:50) {
  for (i in 1:n) {
    PHI[i, ] <- mvrnorm(1, prox[i, ] %*% PHI/n.i[i], SIGMA/n.i[i])
  }
  PHI <- PHI - matrix(colMeans(PHI), n, 2, byrow = T)
}

## Generate Ys
Y <- cbind(1, scale(X)) %*% BETA + PHI + mvrnorm(n, c(0, 0), V)
Y[Y < 0] <- 0
Y <- Y^2
Y <- cbind(Y, 1)/(1 + rowSums(Y))

rm(n)
rm(n.i)