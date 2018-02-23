GeneQC <- function(dat){

  ran <- 2:3

  pen <- 10

  colnames(dat) <- c('Gene ID', 'D1', 'D2', 'D3', 'D4', 'AD1')

  dat$D4 <- (dat$D4)/2 #Modify D4

  dat$D1 <- dat$AD1 #Modify D1

  dat <- dat[, 1:5] #Remove original D1

  datz <- dat[dat$D2 == 0, ] #Remove genes with no mapping uncertainty
  dat <- dat[dat$D2 != 0, ]

  mod <- lm(dat$D4 ~ dat$D1*dat$D2*dat$D3)
  smod <- lm(dat$D4 ~ dat$D3:dat$D1 + dat$D3:dat$D2)

  # Determine significant coefficients
  sig <- function(md){
    colmn <- rep(0, dim(summary(md)$coefficients)[1])
    for(i in 1:dim(summary(md)$coefficients)[1]){
      if(summary(md)$coefficients[i,4] < 0.05){
        colmn[i] <- (summary(md)$coefficients[i,1])
      }
      else(colmn[i] <- (0))
    }
    return(colmn)
  }

  dat_co <- sig(mod)
  dat_sco <- sig(smod)

  #D-score
  dscore <- function(x, y){
    d <- rep(NA, dim(x)[1])
    for(i in 1:dim(x)[1]){
      d[i] <- y[1] + y[2]*x$D1[i] + y[3]*x$D2[i] + y[4]*x$D3[i] + y[5]*x$D1[i]*x$D2[i] + y[6]*x$D1[i]*x$D3[i] + y[7]*x$D2[i]*x$D3[i] + y[8]*x$D1[i]*x$D2[i]*x$D3[i]
    }
    return(d)
  }

  dscoresim <- function(x, y){
    sd <- rep(NA, dim(x)[1])
    for(i in 1:dim(x)[1]){
      sd[i] <- y[1] + y[2]*x$D1[i]*x$D3[i] + y[3]*x$D2[i]*x$D3[i]
    }
    return(sd)
  }

  d <- dscore(dat, dat_co)
  sd <- dscoresim(dat, dat_sco)

  # Mixture Model Fitting

  #######################################################################################################

  #######################################################################################################



  M <- d

  M2 <- M
  if (min(M)<0){
    M2 <- M + abs(min(M)) + 0.000001 # make sure all values are greater than 0
  }

  # Mixture Normal Distribution Fitting

  set.seed(123456)
  B <-list()
  for (KK in ran){
    B[[KK]] <- mixtools::normalmixEM(M,k=KK)
  }

  n <- length(M)   # number of observations
  Logl <-as.numeric()
  KK <- as.numeric()
  for (i in ran) {
    if (length(B[[i]])!=0){
      k <-i
      ll <-B[[i]]$loglik
      Logl <-c(ll,Logl)
      KK <-c(k,KK)
    }
  }

  BIC_norm <-pen*KK*log(n)-2*Logl
  DF_Norm <-data.frame(LOGLIK = Logl,K = KK,BIC = BIC_norm)
  OPK_norm <-DF_Norm[which(DF_Norm$BIC ==min(DF_Norm$BIC)),]$K
  BIC_norm_opt <-DF_Norm[which(DF_Norm$BIC ==min(DF_Norm$BIC)),]$BIC

  PARS <- rbind(B[[OPK_norm]]$mu,  B[[OPK_norm]]$sigma)
  row.names(PARS) <- c('mu', 'sigma')


  # Mixture Gamma Distribution Fitting
  # k-means initialization
  set.seed(1234567)
  x <-M2

  ## initially grouping data using kmeans  ##

  kmeans.init <- function(x, lambda = NULL, alpha = NULL, beta = NULL, k = 2){
    kfit <-kmeans(x,k)   # kmeans clustering
    df <-data.frame(x,K=kfit$cluster) # save values and cluster label together
    names(df) <-c("Dscore","K")

    if(k==1){
      x.bar=mean(x)
      x2.bar=mean(x^2)
    } else{
      x.part <-list()
      lambda <-as.numeric()
      for (i in 1:k){
        x.part[[i]] <-subset(df,K==i,select=Dscore:K)$Dscore
        lambda[i] <-length(x.part[[i]])/dim(df)[1]
      }

      x.bar=sapply(x.part,mean)
      x2.bar=sapply(lapply(x.part,"^",2),mean)
    }
    if(is.null(alpha)){
      alpha=x.bar^2/(x2.bar-x.bar^2)
    }

    if(is.null(beta)){
      beta=(x2.bar-x.bar^2)/x.bar
    }

    list(lambda=lambda, alpha=alpha, beta=beta, k=k)
  }

  # 1121 NOTE :Change epsilon from 1e-08 to 1e-06
  ##  update the fitting with mixture of Gamma distribution
  GammamixEM <- function (x, lambda = NULL, alpha = NULL, beta = NULL, k = 2,
                          epsilon = 1e-06, maxit = 1000, maxrestarts=10, verb = FALSE) {
    x <- as.vector(x)
    tmp <- kmeans.init(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k)    # initial grouping based on Kmeans
    lambda <- tmp$lambda
    alpha <- tmp$alpha
    beta <- tmp$beta
    theta <- c(alpha,beta)
    k <- tmp$k
    iter <- 0
    mr <- 0
    diff <- epsilon+1
    n <- length(x)
    dens <- NULL
    dens <- function(lambda, theta, k){  # this function obtains the pdf for mixture of Gamma
      temp<-NULL
      alpha=theta[1:k]
      beta=theta[(k+1):(2*k)]
      for(j in 1:k){
        temp=cbind(temp,dgamma(x,shape=alpha[j],scale=beta[j]))  }  #  pdf of Gamma
      temp=t(lambda*t(temp))
      temp
    }
    old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))    ## log likelihood function for mixture of Gamma
    ll <- old.obs.ll
    gamma.ll <- function(theta, z,lambda, k) -sum(z*log(dens(lambda,theta,k)))  # define an objective function for the following nlm()
    while(diff > epsilon && iter < maxit){
      dens1=dens(lambda,theta,k)
      z=dens1/apply(dens1,1,sum)  # z is the probability that x belong to ith compoment
      lambda.hat=apply(z,2,mean)
      out=try(suppressWarnings(nlm(gamma.ll,p=theta,lambda=lambda.hat,k=k,z=z)),   # nlm try to find optimal paramters for gamma.all
              silent=TRUE)
      if(class(out)=="try-error"){
        cat("Note: Choosing new starting values.", "\n")
        if(mr==maxrestarts) break(paste("Trying different number of components?","\n"))
        mr <- mr+1
        tmp <- kmeans.init(x = x, k=k)
        lambda <- tmp$lambda
        alpha <- tmp$alpha
        beta <- tmp$beta
        theta <- c(alpha,beta)
        k <- tmp$k
        iter <- 0
        diff <- epsilon+1
        old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
        ll <- old.obs.ll
      } else{
        theta.hat=out$estimate
        alpha.hat=theta.hat[1:k]
        beta.hat=theta.hat[(k+1):(2*k)]
        new.obs.ll <- sum(log(apply(dens(lambda.hat, theta.hat, k),1,sum)))
        diff <- new.obs.ll-old.obs.ll
        old.obs.ll <- new.obs.ll
        ll <- c(ll,old.obs.ll)
        lambda=lambda.hat
        theta=theta.hat
        alpha=alpha.hat
        beta=beta.hat
        iter=iter+1
        if (verb) {
          cat("iteration =", iter, " log-lik diff =", diff, " log-lik =",
              new.obs.ll, "\n")
        }
      }
    }
    if (iter == maxit) {
      cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    theta=rbind(alpha,beta)
    rownames(theta)=c("alpha","beta")
    colnames(theta)=c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll,
           posterior = z, all.loglik=ll, ft="gammamixEM")
    class(a) = "mixEM"
    a
  }

  # fit the data with different number of components
  A <-list()
  for (KK in ran){
    A[[KK]] <- GammamixEM(x,k=KK)
  }

  n <- length(M2)   # num of obesrvations
  DF_LL <- as.data.frame(do.call(rbind,lapply(A,'[[',4)))  # extract the log likelihood

  if(sum(DF_LL < 0) < 2){
    subll <- which(DF_LL >= 0)
  } else{
    subll <- which(DF_LL < 0)  # confine to negative loglikelihood
  }

  bic <- c(subll+1)*pen*log(n)-2*DF_LL[subll,]
  # bic2 <-c(2:9)*3*log(n)-2*DF_LL
  # OPTK2 <-which(bic==min(bic))+1
  OPTK <- subll[which(bic==min(bic))]+1  # choose the K than minimize bic
  THETAS <-A[[OPTK]]$gamma.pars
  lambdas <-A[[OPTK]]$lambda
  BIC_gamma_opt <-min(bic)

  if (BIC_norm_opt < BIC_gamma_opt){
    order(PARS[1,])
    oPARS <- PARS[,order(PARS[1,])]

    m <- oPARS[1,]
    sd <- oPARS[2,]
    x <- rep(NA, length(m)-1)
    for(i in 1:(length(m)-1)){
      x1 <- -(m[i+1]*sd[i]^2-m[i]*sd[i+1]^2)/(sd[i+1]^2-sd[i]^2)+sqrt((2*sd[i]^2*sd[i+1]^2*log(sd[i+1]/sd[i])-m[i]^2*sd[i+1]^2+m[i+1]^2*sd[i]^2)/(sd[i+1]^2-sd[i]^2)+((m[i+1]*sd[i]^2-m[i]*sd[i+1]^2)/(sd[i+1]^2-sd[i]^2))^2)
      x2 <- -(m[i+1]*sd[i]^2-m[i]*sd[i+1]^2)/(sd[i+1]^2-sd[i]^2) - sqrt((2*sd[i]^2*sd[i+1]^2*log(sd[i+1]/sd[i])-m[i]^2*sd[i+1]^2+m[i+1]^2*sd[i]^2)/(sd[i+1]^2-sd[i]^2)+((m[i+1]*sd[i]^2-m[i]*sd[i+1]^2)/(sd[i+1]^2-sd[i]^2))^2)
      if(x1 < m[i+1] & x1 > m[i]){
        x[i] <- x1
      } else{
        x[i] <- x2
      }
    }
    cutoffs <- x
  } else {
    gamma_mode1 <- function(x){
      if(x[1] > 1){
        mo <- (x[1] - 1)*x[2]
      }
      else{
        mo <- 0
      }
      mo
    }

    oTHETAS <- THETAS[,order(apply(FUN =  gamma_mode1,MARGIN =  2, X = THETAS))]

    cutoffs <- rep(NA, (dim(oTHETAS)[2]-1))
    for(i in 1:(dim(oTHETAS)[2]-1)){
      modes <- apply(FUN =  gamma_mode1 ,MARGIN =  2, X = oTHETAS[,i:(i+1)])
      gs <- seq(min(modes), max(modes), length.out = 1000)
      lv <- dgamma(gs, shape = oTHETAS[1,i], scale = oTHETAS[2,i]) > dgamma(gs, shape = oTHETAS[1,i+1], scale = oTHETAS[2,i+1])
      cutoffs[i] <- mean(gs[min(which(lv != lv[1]))], gs[min(which(lv != lv[1]))+1])
    }
  }

  if( (abs(cutoffs[2] - cutoffs[1])/(range(d)[2] - range(d)[1])) < 0.20 ){
    cutoffs[1] <- (min(d) + cutoffs[2])/2
  }

  gamma_mode1 <- function(x){
    if(x[1] > 1){
      mo <- (x[1] - 1)*x[2]
    }
    else{
      mo <- 0
    }
    mo
  }

  if (BIC_norm_opt < BIC_gamma_opt){
    parameters <- PARS[,order(PARS[1,])]
  } else{
    parameters <- THETAS[,order(apply(FUN =  gamma_mode1,MARGIN =  2, X = THETAS))]
  }

  list1 <- list(d, sd, THETAS[,order(apply(FUN =  gamma_mode1,MARGIN =  2, X = THETAS))], PARS[,order(PARS[1,])], cutoffs, parameters)

  categ <- rep(NA, length(list1[[1]]))
  for(i in 1:length(list1[[1]])){
    if(list1[[1]][i] < list1[[5]][1]){
      categ[i] <- 'Low'
    } else if(list1[[1]][i] >= list1[[5]][1] & list1[[1]][i] < list1[[5]][2]){
      categ[i] <- 'Medium'
    } else if(list1[[1]][i] >= list1[[5]][1]){
      categ[i] <- 'High'
    }
  }

  sign <- rep(NA, length(list1[[1]]))
  for(i in 1:length(list1[[1]])){
    if(categ[i] == 'Low'){
      if(all(list1[[6]] == list1[[3]])){
        sign[i] <- pgamma(list1[[1]][i], shape = list1[[6]][1,2], scale = list1[[6]][2,2])
      } else{
        sign[i] <- pnorm(list1[[1]][i], mean = list1[[6]][1,2], sd = list1[[6]][2,2])
      }
    } else if(categ[i] == 'High'){
      if(all(list1[[6]] == list1[[3]])){
        sign[i] <- 1-pgamma(list1[[1]][i], shape = list1[[6]][1,2], scale = list1[[6]][2,2])
      } else{
        sign[i] <- 1-pnorm(list1[[1]][i], mean = list1[[6]][1,2], sd = list1[[6]][2,2])
      }
    } else if(categ[i] == 'Medium'){
      if(all(list1[[6]] == list1[[3]])){
        sign[i] <- max(1-pgamma(list1[[1]][i], shape = list1[[6]][1,1], scale = list1[[6]][2,1]), pgamma(list1[[1]][i], shape = list1[[6]][1,3], scale = list1[[6]][2,3]))
      } else{
        sign[i] <- max(1-pnorm(list1[[1]][i], mean = list1[[6]][1,1], sd = list1[[6]][2,1]), pnorm(list1[[1]][i], mean = list1[[6]][1,3], sd = list1[[6]][2,3]))
      }
    }
  }

  list2 <- list(d, sd, THETAS[,order(apply(FUN =  gamma_mode1,MARGIN =  2, X = THETAS))], PARS[,order(PARS[1,])], cutoffs, parameters, categ, sign, dat, datz)

  .res <- rbind(list2[[9]], list2[[10]])
  .res <- cbind(.res, c(list2[[1]], rep(0, dim(list2[[10]])[1])))
  .res <- cbind(.res, c(list2[[7]], rep('None', dim(list2[[10]])[1])))
  .res <- cbind(.res, c(list2[[8]], rep(0, dim(list2[[10]])[1])))
  colnames(.res) <- c('Gene ID', 'D1', 'D2', 'D3', 'D4', 'D-score', 'Category', 'Significance')

  list3 <- list(d, sd, THETAS[,order(apply(FUN =  gamma_mode1,MARGIN =  2, X = THETAS))], PARS[,order(PARS[1,])], cutoffs, parameters, categ, sign, .res[,c(1:4,6:8)])
  list3[[9]]
}
