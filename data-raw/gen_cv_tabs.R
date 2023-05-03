
seed <- 250493
workers <- 7
`%dopar%` <- foreach::`%dopar%`
#' -----------------------------------------------------------------------------
# Critical values for supF & supLR when q= [55, 61], k = [1, 9], & trm = 0.05
q_min <- 55
q_max <- 61
k_max <- 9
n <- 1000
B <- 10000
trm <- 0.05
supF_cv_05 <- matrix(0,4*((q_max-q_min)+1),k_max)
doParallel::registerDoParallel(workers)
for (k in 1:k_max){
  set.seed(seed)
  out_tmp <- foreach::foreach(wi = 1:workers ,.inorder = FALSE, .packages = "multiMNbreaks") %dopar% {
    cv_tmp <- bp_supF_critval(k, q_min+wi-1, trm, n, B)    
    sort(cv_tmp)[round(c(0.9,0.95,0.975,0.99)*B)]
  }
  for (xi in 1:((q_max-q_min)+1)){
    supF_cv_05[((xi-1)*4+1):(xi*4),k] <- out_tmp[[xi]]
  }
}
doParallel::stopImplicitCluster()

rownames(supF_cv_05) <- paste0(sort(rep(q_min:q_max,4)),rep(c(" - 90.0%"," - 95.5%"," - 97.5%"," - 99.0%"),q_max-q_min+1))
colnames(supF_cv_05) <- 1:k_max
supLR_cv_05 <- t(matrix(1:k_max,k_max,4*(q_max-q_min+1)))*matrix(sort(rep(q_min:q_max,4)),4*(q_max-q_min+1),k_max)*supF_cv_05

usethis::use_data(supF_cv_05, overwrite = TRUE)
usethis::use_data(supLR_cv_05, overwrite = TRUE)

#' -----------------------------------------------------------------------------
# Critical values for seqF when q= [55, 61], k = [1, 9], & trm = 0.05
q_min <- 55
q_max <- 61
k_max <- 9
n <- 1000
B <- 10000
trm <- 0.05
seqF_cv_05 <- matrix(0,4*((q_max-q_min)+1),k_max)
doParallel::registerDoParallel(workers)
for (k in 1:k_max){
  set.seed(seed)
  out_tmp <- foreach::foreach(wi = 1:workers ,.inorder = FALSE, .packages = "multiMNbreaks") %dopar% {
    cv_tmp <- bp_supSeq_critval(k,q_min+wi-1,trm,n,B)
    sort(cv_tmp)[round(c(0.9,0.95,0.975,0.99)*B)]
  }
  for (xi in 1:((q_max-q_min)+1)){
    seqF_cv_05[((xi-1)*4+1):(xi*4),k] <- out_tmp[[xi]]
  }
}
doParallel::stopImplicitCluster()

rownames(seqF_cv_05) <- paste0(sort(rep(q_min:q_max,4)),rep(c(" - 90.0%"," - 95.5%"," - 97.5%"," - 99.0%"),q_max-q_min+1))
colnames(seqF_cv_05) <- 1:k_max
seqF_cv_05
seqLR_cv_05 <- t(matrix(1:k_max,k_max,4*(q_max-q_min+1)))*matrix(sort(rep(q_min:q_max,4)),4*(q_max-q_min+1),k_max)*seqF_cv_05

usethis::use_data(seqF_cv_05, overwrite = TRUE)
usethis::use_data(seqLR_cv_05, overwrite = TRUE)

#' -----------------------------------------------------------------------------
# Critical values for supF & supLR when q= [55, 61], k = [1, 9], & trm = 0.10
q_min <- 55
q_max <- 61
k_max <- 9
n <- 1000
B <- 10000
trm <- 0.10

st <- proc.time()
supF_cv_10 <- matrix(0,4*((q_max-q_min)+1),k_max)
doParallel::registerDoParallel(workers)
for (k in 1:k_max){
  set.seed(seed)
  out_tmp <- foreach::foreach(wi = 1:workers ,.inorder = FALSE, .packages = "multiMNbreaks") %dopar% {
    cv_tmp <- bp_supF_critval(k, q_min+wi-1, trm, n, B)    
    sort(cv_tmp)[round(c(0.9,0.95,0.975,0.99)*B)]
  }
  for (xi in 1:((q_max-q_min)+1)){
    supF_cv_10[((xi-1)*4+1):(xi*4),k] <- out_tmp[[xi]]
  }
}
doParallel::stopImplicitCluster()
end <- proc.time() - st
end

rownames(supF_cv_10) <- paste0(sort(rep(q_min:q_max,4)),rep(c(" - 90.0%"," - 95.5%"," - 97.5%"," - 99.0%"),q_max-q_min+1))
colnames(supF_cv_10) <- 1:k_max
supLR_cv_10 <- t(matrix(1:k_max,k_max,4*(q_max-q_min+1)))*matrix(sort(rep(q_min:q_max,4)),4*(q_max-q_min+1),k_max)*supF_cv_10

usethis::use_data(supF_cv_10, overwrite = TRUE)
usethis::use_data(supLR_cv_10, overwrite = TRUE)

#' -----------------------------------------------------------------------------
# Critical values for seqF when q= [55, 61], k = [1, 9], & trm = 0.10
q_min <- 55
q_max <- 61
k_max <- 9
n <- 1000
B <- 10000
trm <- 0.10
seqF_cv_10 <- matrix(0,4*((q_max-q_min)+1),k_max)
doParallel::registerDoParallel(workers)
for (k in 1:k_max){
  set.seed(seed)
  out_tmp <- foreach::foreach(wi = 1:workers ,.inorder = FALSE, .packages = "multiMNbreaks") %dopar% {
    cv_tmp <- bp_supSeq_critval(k,q_min+wi-1,trm,n,B)
    sort(cv_tmp)[round(c(0.9,0.95,0.975,0.99)*B)]
  }
  for (xi in 1:((q_max-q_min)+1)){
    seqF_cv_10[((xi-1)*4+1):(xi*4),k] <- out_tmp[[xi]]
  }
}
doParallel::stopImplicitCluster()

rownames(seqF_cv_10) <- paste0(sort(rep(q_min:q_max,4)),rep(c(" - 90.0%"," - 95.5%"," - 97.5%"," - 99.0%"),q_max-q_min+1))
colnames(seqF_cv_10) <- 1:k_max
seqLR_cv_10 <- t(matrix(1:k_max,k_max,4*(q_max-q_min+1)))*matrix(sort(rep(q_min:q_max,4)),4*(q_max-q_min+1),k_max)*seqF_cv_10

usethis::use_data(seqF_cv_10, overwrite = TRUE)
usethis::use_data(seqLR_cv_10, overwrite = TRUE)








