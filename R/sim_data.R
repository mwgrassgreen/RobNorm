#' @title sim.dat.fn
#' @description To simulation an expression matrix
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param row.frac Outlier fraction in the rows of the expression matrix.
#' @param col.frac Outlier fraction in the columns of the expression matrix.
#' @param mu.up The up-shifted mean of the outliers.
#' @param mu.down The down-shifted mean of the outliers.
#' @param n Number of rows (genes).
#' @param m Number of colunms (samples).
#' @param nu.fix Logic value indicating underying nu=0 (default: TRUE).
#' @return Simulated data and its parameters.
#' @import invgamma
#' @export


sim.dat.fn = function(row.frac, col.frac, mu.up, mu.down, n, m, nu.fix=TRUE) {
	mu.00 = rnorm(n, 0, 1)
	var.00 = rinvgamma(n, 5, 2)
    X.0 = matrix(rnorm(n*m, outer(mu.00, rep(1, m)), sqrt(outer(var.00, rep(1, m)))), n, m) # the null matrix
    
    if (nu.fix) {
    	nu.00 = rep(0,m)
    } else {
       nu.00 = rnorm(m, 0, 1)
	   nu.00[1:round(m*0.2)] = rnorm(round(m * 0.2), 1, 1)	
    }
    B = outer(rep(1, n), nu.00) # the sample effect matrix

	S = matrix(0, n, m) 
	if (row.frac*col.frac > 0) {
	   	   bk.nm = round(n*row.frac * m*col.frac)
		   a = rbinom(1, bk.nm, 0.8)
		   S[ 1:round(n*row.frac), 1:(m*col.frac)] = sample(c(rep(mu.up, a), rep(0, bk.nm-a)), bk.nm) # the shifted mean of the signal mx
		   a = rbinom(1, bk.nm, 0.8)
		   S[ (n-round(n*row.frac)+1):n, (m-m*col.frac+1):m] = sample(c(rep(mu.down, a), rep(0, bk.nm-a)), bk.nm) # the signal mx of shifted mean 		  
	}
	
	X = X.0 + B + S
	S.ind = matrix(0, nrow(X), ncol(X))
	S.ind[S != 0] = 1
	
	rownames(X) = paste("prt", 1:nrow(X), sep=".")
	colnames(X) = paste("s", 1:ncol(X), sep=".")
    return(list(dat=X, mu.00=mu.00, var.00=var.00, nu.00=nu.00, sig.ind=S.ind))	
}

