
unbiased_stICA <- function(X,k=10,alpha) {

  #
  #    [A B W] = unbiased_stICA(X, k,alpha)
  #
  # Compute factorization of X =A*B' using Jade algorithm
  #
  # Inputs:
  # - X (matrix p*n) is the matrix to factorise
  # - alpha (scalar between 0 and 1) is the trade-off between spatial ICA (alpha = 0) et temporal ICA
  # (alpha=1) (default alpha =0.5)
  # - k (scalar integer) is the number of components to estimate (default = 10)
  #
  # Ouputs:
  # - A (matrix p*k), B (matrix n*k) such that A*B' is an approximation of X
  # - W orthogonal matrix k*k minimizing the objective function :
  #       f(W) = sum_i ( alpha ||off(W*C_i( Dk^(1-alpha)*Vk)* W')||_F^2
  #                   + (1-alpha)||off( W* C_i(Dk^alpha*Uk)*W') ||_F^2  )
  #
  #       where - C_i(X) is fourth order cumulant-like matrices of X
  #             - X = U* D*V^T is the svd decomposition of X, U/D/V_k is
  #              the U/D/V matrix keeping only the k first components

  # Author: Emilie Renard, december 2013

  # References:
  #  E. Renard, A. E. Teschendorff, P.-A. Absil
  # 'Capturing confounding sources of variation in DNA methylation data by
  #  spatiotemporal independent component analysis', submitted to ESANN 2014
  #  (Bruges)
  #
  #  based on code from F. Theis, see
  #  F.J. Theis, P. Gruber, I. Keck, A. Meyer-Baese and E.W. Lang,
  #  'Spatiotemporal blind source separation using double-sided approximate joint diagonalization',
  #  EUSIPCO 2005 (Antalya), 2005.
  #
  #  Uses Cardoso's diagonalization algorithm based on iterative Given's
  #  rotations (matlab-file RJD), see
  #  J.-F. Cardoso and A. Souloumiac, 'Jacobi angles for simultaneous diagonalization',
  #  SIAM J. Mat. Anal. Appl., vol 17(1), pp. 161-164, 1995

  library(JADE)
  library(corpcor)

  jadeCummulantMatrices <- function(X) {
    # calcs the n(n+1)/2 cum matrices used in JADE, see A. Cichocki and S. Amari. Adaptive blind signal and image
    # processing. John Wiley & Sons, 2002. book, page 173
    # (pdf:205), C.1
    # does not need whitened data X (n x t)
    # Args:
    #   matrix X
    # Returns:
    #   M as n x n x (n(n+1)/2) array, each n x n matrix is one of the n(n+1)/2 cumulant matrices used in JADE

    # adapted code from Matlab implementation of F. Theis, see
    #  F.J. Theis, P. Gruber, I. Keck, A. Meyer-Baese and E.W. Lang,
    #  'Spatiotemporal blind source separation using double-sided approximate joint diagonalization',
    #  EUSIPCO 2005 (Antalya), 2005.


    n <- nrow(X)
    t <- ncol(X)

    M <- array(0,c(n,n,n*(n+1)/2))
    scale <- matrix(1,n,1)/t  # for convenience


    R <- cov(t(X)) # covariance

    k <- 1
    for (p in 1:n){
      #case q=p
      C <- ((scale %*% (X[p,]*X[p,]))*X) %*% t(X)
      E <- matrix(0,n,n)
      E[p,p] <- 1
      M[,,k] <- C - R %*% E %*% R - sum(diag(E %*% R)) * R - R %*% t(E) %*% R
      k <- k+1
      #case q<p
      if (p > 1) {
        for (q in 1:(p-1)){
          C <- ((scale %*% (X[p,]*X[q,]))*X) %*% t(X) * sqrt(2)
          E <- matrix(0,n,n)
          E[p,q] <- 1/sqrt(2)
          E[q,p] <- E[p,q]
          M[,,k] <- C - R %*% E %*% R - sum(diag(E %*% R)) * R - R %*% t(E) %*% R
          k <- k+1
        }
      }
    }

    return(M)

  }






  p <- nrow(X)
  n <- ncol(X)


  dimmin <- min(n,p)

  if (dimmin < k) {
    k <- dimmin
  }
  if (alpha <0 | alpha >1){
    stop("alpha not in [0 1]")
  }

  # Remove the spatiotemporal mean

  Xc <- X - matrix(rep(colMeans(X,dims=1),p),nrow = p,byrow=T);
  Xc <- Xc - matrix(rep(rowMeans(Xc,dims=1),n),nrow = p);

  # SVD of Xc and dimension reduction: keeping only the k first
  # components
  udv <- svd(Xc,k,k)
  D <- diag(udv$d[1:k]); if (k==1) {D <- udv$d[1]}
  U <- udv$u;
  V <- udv$v;


  # Estimation of the cumulant matrices
  nummat <- k*(k+1)/2;
  M <- array(0,c(k,k,2*nummat));
  Bt <- D^(1-alpha) %*% t(V)
  if (alpha == 1) { Bt <- t(V)}
  At <- D^(alpha) %*% t(U)
  if (alpha == 0) { At <- t(U)}
  M[,,1:nummat] <- jadeCummulantMatrices(Bt);
  M[,,(nummat+1):(2*nummat)] <- jadeCummulantMatrices(At)




  # normalization within the groups in order to allow for comparisons using
  # alpha
  M[,,1:nummat] <- alpha*M[,,1:nummat]/mean(sqrt(apply(M[,,1:nummat]*M[,,1:nummat],3,sum)));
  M[,,(nummat+1):(2*nummat)] <- (1-alpha)*M[,,(nummat+1):(2*nummat)]/mean(sqrt(apply(M[,,(nummat+1):(2*nummat)]*M[,,(nummat+1):(2*nummat)],3,sum)));



  # Joint diagonalization
  Worth <- rjd(M,eps = 1e-06, maxiter = 1000);
  Wo <-t (Worth$V);


  #     Computation of A and B

  A0 <- U %*% D^(alpha) %*% solve(Wo);
  B0 <- V%*% D^(1-alpha) %*% t(Wo);
  if (alpha == 1) { B0 <- V %*% t(Wo)}
  if (alpha == 0) { A0 <- U %*% solve(Wo)}

  # Add transformed means
  meanCol <- matrix(colMeans(X,dims=1),ncol =1); # spatial means
  meanRows <- matrix(rowMeans(X,dims=1),ncol = 1); # temporal means

  meanB <- pseudoinverse(A0) %*% (meanRows);
  meanA <- pseudoinverse(B0) %*% (meanCol);

  Bfin <- B0 + matrix(rep(meanB,n),nrow = n,byrow=T)
  Afin <- A0 + matrix(rep(meanA,p),nrow = p,byrow=T)



  return(list(A=Afin,B=Bfin,W=Wo))



}
