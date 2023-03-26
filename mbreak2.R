

library(multiMNbreaks)


# Set your T by n matrix of dependent variables@                                  
y <- data1[,1:2]

# Set your T by q matrix of regressors  
z <- data1[,3:5]

# Specify the number of observations
Ts <- nrow(y)

# Specify the matrix S (here the first and second regressors are used in the 
# first equation and the first and third regressors are used in the second equation)@                      
S <- rbind(c(1,0,0,0),
           c(0,1,0,0),
           c(0,0,0,0),
           c(0,0,1,0),
           c(0,0,0,0),
           c(0,0,0,1))

# set the maximum number of breaks allowed. If use the WDmax test set m=5, since the critical values reported
# correspond to this case
m <- 2


# If you have restrictions, specify the matrix R (here only the constants (first
# regressor) is allowed to change in each equation). 
# if you do not have restrictions, set R=eye(cols(_S)*(m+1))@
# NOTE: This R is for when m=2
R <- rbind(c(1,0,0,0,0,0,0,0),
           c(0,1,0,0,0,0,0,0),
           c(0,0,1,0,0,0,0,0),
           c(0,0,0,1,0,0,0,0),
           c(0,0,0,0,1,0,0,0),
           c(0,1,0,0,0,0,0,0),
           c(0,0,0,0,0,1,0,0),
           c(0,0,0,1,0,0,0,0),
           c(0,0,0,0,0,0,1,0),
           c(0,1,0,0,0,0,0,0),
           c(0,0,0,0,0,0,0,1),
           c(0,0,0,1,0,0,0,0))



# set the trimming parameter that specifies the minimal length of a segment as a proportion of the sample size.
trm <- 0.2

# set brv=1 when allowing breaks in the covariance matrix of the errors; otherwise, set brv=0
brv <- 0

# set brbeta=1 when allowing breaks in regression coefficients; otherwise, set brbeta=0.                  
brbeta <- 1





# set vauto=1 when applying a correction for serial correlation in the errors
# (in which case Andrews' (1991) method is used to construct the robust 
# covariance matrix. 
vauto <- 0

# set to prewhit=1 if you want to apply an AR(1) pre-whitening when estimating the 
# robust covariance matrix (see Andrews and Monahan (1992).
prewhit <- 0

# Set hetq=1 if you want to allow the distribution of the regressors to change across
# regimes (this is used only when constructing confidence intervals for the estimates 
# of the break dates. For the construction of the test hetq=1 always.
hetq <- 0

# /*Call the main procedure*/
#res <- mainp(m,ncol(S),z,y,ncol(y),trm,Ts,brv,brbeta,vauto,hetq,prewhit)

n <- ncol(y)


h <- round(trm*Ts)  #trimming
# transform the data
transf_out <- transf(y,z,S)

maty <- transf_out$maty
matz <- transf_out$matz

i <- 1 # number of breaks

mdll <- est(maty,matz,n,i,Ts,trm,diag(ncol(matz)*(i+1)),brbeta, brv)















