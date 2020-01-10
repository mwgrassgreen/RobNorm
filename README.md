# RobNorm
RobNorm is a R package to robustly normalize an expression matrix, analyzed in the paper "RobNorm: Model-Based Robust Normalization for High-Throughput Proteomics from Mass Spectrometry Platform". 

Please contact Meng Wang by email <mengw1@stanford.edu> for questions. 

## Installation
`library(devtools)`

`install_github("mwgrassgreen/RobNorm")`

## Usage
`library(RobNorm)`

`norm.result = RobNorm(X.0, gamma.0=0.5, tol=10^(-4), step=200)`

## Example
To simulate an expression matrix 

`sim.result = sim.dat.fn(row.frac=0.2, col.frac=0.2, mu.up=3, mu.down=-3, n=5000, m=200,  nu.fix=TRUE)`
 
`X.0 = sim.result$dat`

`norm.result = RobNorm(X.0, gamma.0=0.5)`

`X.0.norm = norm.result$norm.data`

To compare sample boxplots before and after normalization

`par(mfrow=c(2,1))`

`boxplot(X.0, main="Sample boxplots before normalization", ylab="expression", xlab="sample", cex.main=1.5, cex.lab=1.5)`

`boxplot(X.0.norm, main="Sample boxplots after normalization", ylab="expression", xlab="sample", cex.main=1.5, cex.lab=1.5)`

Since in the simulation nu.fix=TRUE meaning the underlying nu = 0, the boxplots before and after normalization.