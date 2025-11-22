library(mvtnorm)
library(MASS)
#generate TOX and ORR the same way in DODII using multinomial 
gen_pis <- function(p_ORR, p_TOX, phi) {
  stopifnot(phi > 0)
  pi <- matrix(NA, ncol = 1, nrow = 4)
  quad_a <- 1 - phi
  quad_b <- 1 - (1 - phi) * (p_ORR + p_TOX)
  quad_c <- -phi * p_ORR * p_TOX
  
  if (phi != 1) {
    pi[1,1] <- (-quad_b + sqrt(quad_b^2 - 4 * quad_a * quad_c)) / (2 * quad_a)
  } else {
    pi[1,1] <- p_ORR * p_TOX
  }
  
  pi[2,1] <- p_ORR - pi[1,1]  
  pi[3,1] <- p_TOX - pi[1,1] 
  pi[4,1] <- 1 - sum(pi[1:3,1])
  return(pi)
}

gen_phi <- function(p_ORR, p_TOX, rho_tox_eff){
  # Convert to thresholds
  theta_tox <- qnorm(1 - p_TOX)
  theta_eff <- qnorm(1 - p_ORR)
  
  # Compute bivariate normal CDF for threshold region
  Phi2 <- function(q1, q2, rho) {
    pmvnorm(lower = c(-Inf, -Inf), upper = c(q1, q2), corr = matrix(c(1, rho, rho, 1), 2))
  }
  
  # Compute joint probabilities
  pi00 <- Phi2(theta_tox, theta_eff, rho_tox_eff)
  pi10 <- pnorm(theta_tox) - pi00
  pi01 <- pnorm(theta_eff) - pi00
  pi11 <- 1 - pnorm(theta_tox) - pnorm(theta_eff) + pi00
  
  # Compute phi
  phi <- (pi11 * pi00) / (pi10 * pi01)
  phi <- as.numeric(phi)
}

sim_data= function(N1,N2,lambda,seed,
                   p_ORR,p_TOX,Sigma){ 
  #N1: sample size for each arm in stage 1
  #N2: sample size for each arm in stage 2
  #lambda_c: parameter of exponential distribution in control arm
  #HR: a vector of c(HR_low, HR_high)
  #H: null (type I error) or alternative (power)
  #Sigma: correlation matrix of [tox,eff,surv]
  #Sigma is the same across doses, but phi is different across doses in DODII
  
  N=N1+N2
  colname =c("Time","Group","Stage","orr","tox");
  arm0 = as.data.frame(matrix(0,N,length(colname))); colnames(arm0)=colname;
  arm0$Stage[1:N1]=1; arm0$Stage[(N1+1):N]=2;
  arm1=arm0; arm2=arm0 
  
  arm0$Group=0 # control
  arm1$Group=1 # low dose
  arm2$Group=2 # high dose
  
  for (dose in 1:3) {
    set.seed(dose+seed)
    Z <- mvrnorm(N, mu = rep(0, 3), Sigma = Sigma)
    
    # Convert to binary toxicity and efficacy using thresholds
    theta_tox <- qnorm(1 - p_TOX[dose])
    theta_eff <- qnorm(1 - p_ORR[dose])
    
    tox <- as.numeric(Z[, 1] > theta_tox)
    eff <- as.numeric(Z[, 2] > theta_eff)
    U <- pnorm(Z[, 3])
    surv <- -log(U) / lambda[dose]
    
    # Assign to the correct arm
    if (dose == 1) {
      arm0$tox <- tox
      arm0$orr <- eff
      arm0$Time <- surv
    } else if (dose == 2) {
      arm1$tox <- tox
      arm1$orr <- eff
      arm1$Time <- surv
    } else if (dose == 3) {
      arm2$tox <- tox
      arm2$orr <- eff
      arm2$Time <- surv
    }
    
  }
  
  dat= rbind(arm0,arm1,arm2)
  
  return(dat)
}

surv.test = function(dat, altern="two.sided"){ 
  ##generate one-sdied z statistic/p-value
  data=subset(dat, select=c("Time","ind","Group"))
  survf=logrank.test(data$Time,data$ind,data$Group,alternative = altern)
  pval = survf$test$p
  ## 10/19/2024 Haiyang Sheng: p-value for one-sided log-rank test is just half of the p-value for two-sided log-rank test, if the direction of HR is as expected.
  return(pval)
}

OS_Test = function(pvals, pcutS){
  #OS test for 2 superiority
  ##pvals: all pvals
  ##pcutF: boundary for futility
  ##pcutS: boundary for superiority
  decision=NA
  if (pvals[1]<pcutS[1]){#win at look 1
    decision=1
  }else if (pvals[2]<pcutS[2]){#win at look 2
    decision=1
  }else {
    decision=0
  }
  return(decision)
}

dunnett=function(p, m){
  ##p: p-value to be adjusted
  f =AdjustPvalues(c(p, rep(0.999, m-1)),proc = "DunnettAdj",par = parameters(n = Inf))[1]
  return(f)
}

wInvnorm= function(p,n){
  ##p: vector of pvalues
  ##vector of sample sizes ;same length as p
  N=sum(n)
  p1=p[1]; p2=p[2]; w1= sqrt(n[1]/N); w2= sqrt(n[2]/N)
  f = 1-pnorm(w1*qnorm(1-p1)+w2*qnorm(1-p2)) #
  return(f)
}

BOP2_gatekeeper <- function(num_ORR, num_TOX, L_i, stop_mat) {
  return(
    (num_ORR > stop_mat[L_i, 2])*(num_TOX < stop_mat[L_i, 3])
  )
}

pick_winner <- function(n_High_ORR, n_Low_ORR, n_High_TOX, n_Low_TOX, L_i, 
                        ORR_boundary_interim_tbl_list, TOX_boundary_interim_tbl_list){
  tox = TOX_boundary_interim_tbl_list[[L_i]] %>%
    filter(num_High_TOX == n_High_TOX, num_LOW_TOX == n_Low_TOX) %>%
    pull(TOX_decision)
  
  orr = ORR_boundary_interim_tbl_list[[L_i]] %>%
    filter(num_High_ORR == n_High_ORR, num_LOW_ORR == n_Low_ORR) %>%
    pull(ORR_decision)
  
  if (tox == "TOX_high_sig_not_worse" && orr == "ORR_high_sig_better") {
    return(2)  # High dose wins
  } else if (tox == "TOX_high_sig_worse" && orr == "ORR_high_sig_not_better") {
    return(1)  # Low dose wins
  } else {
    return(c(1, 2))  # No clear winner
  }
  
}