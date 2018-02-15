GeneQC_plot <- function(dat, thickness = 2){
  hist(dat[[1]], breaks = 50, prob = TRUE, main = deparse(substitute(dat)), xlab = 'D-score')
  if(all(dat[[4]] == dat[[6]])){
    for(i in 1:dim(dat[[4]])[2]){
      curve(dnorm(x, mean = dat[[4]][1,i], sd = dat[[4]][2,i]), add = TRUE)
    }
    abline(v = dat[[5]], col = 'red')
  } else if(all(dat[[3]] == dat[[6]])){
    for(i in 1:dim(dat[[4]])[2]){
      curve(dgamma(x, shape = dat[[3]][1,i], scale = dat[[3]][2,i]), add = TRUE)
    }
    abline(v = dat[[5]], col = 'red')
  }
}
