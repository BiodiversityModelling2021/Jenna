library(dplyr)
library(tidyr)
library(ggplot2)

#deterministic growth function
# a is a constant
# s is also a constant
# L is the light data, simulate from uniform


growth <- function(a, s, L){
  (a*L)/((a/s) + L)
}

true_beta <- .3
#values from paper
true_constant_a <-  260.6
true_constant_s <-  7.56
sample_size <- 5000


simulated_tree_data <- tibble(
  L = runif(n = sample_size, min = 0, max = 200),
  growth = rlnorm(n = sample_size,
                    meanlog = log(growth(true_constant_a, true_constant_s, L)), 
                    sdlog = true_beta))


data_figure <- simulated_tree_data %>% 
  ggplot(aes(x = L, y = growth)) + geom_point() + xlab("light availability")
data_figure

data_figure + 
  stat_function(fun = function(x) growth(
    L = x,
    a = true_constant_a,
    s = true_constant_s), lwd = 2, col = "green")

#likelihood function

#testing

loglike_tree <- function(constant_a, constant_s, L, simulated_tree_data, beta){
  alpha <- growth(a = constant_a,
                           s = constant_s,
                           L = L)
  log_like_datapoints <- dlnorm(x = simulated_tree_data, 
                                meanlog = log(alpha), 
                                sdlog = beta, 
                                log = TRUE)
  return(sum(log_like_datapoints))
}

loglike_tree(250, 5, simulated_tree_data$L, simulated_tree_data$growth, 0.2)


#function that takes known value for alpha and gives a beta
#invgamma conjugate to lognormal


sample_beta <- function(y, alpha){
  one_variance <- extraDistr::rinvgamma(1, 
                                        alpha = 50 + length(y)/2,
                                        beta = 4 + sum((log(y) - log(alpha))^2 ) / 2)
  one_sd <- sqrt(one_variance)
  return(one_sd)
}



sample_constantS_MetropolisH <- function(constantS_previous_val,
                                        constantA_known,
                                        beta_known,
                                        light_dat = simulated_tree_data$L,
                                        growth_dat = simulated_tree_data$growth,
                                        MH_tune){
  # 
  
  prev_log_prior_constantS <- calculate_constant_prior_logprob(constantS_previous_val)
  
  prev_log_like_data <- loglike_tree(
    constant_a = constantA_known,
    # the next line is for the parameter we are sampling! 
    constant_s = constantS_previous_val, 
    L = light_dat,
    beta = beta_known,
    simulated_tree_data = growth_dat)
  
  ## propose new!
  constantS_propose_val <- MH_propose_new_val(constantS_previous_val,
                                             sig_tune = MH_tune)
  
  
  propose_log_prior_constantS <- calculate_constant_prior_logprob(constantS_propose_val)
  
  propose_log_like_data <- loglike_tree(
    constant_a = constantA_known,
    # the next line is for the parameter we are sampling! 
    constant_s = constantS_propose_val, 
    L = light_dat,
    beta = beta_known,
    simulated_tree_data = growth_dat)
  
  ## Correction! the special step of MH
  correction <- calculate_correction(new_value = constantS_previous_val,
                                     old_value = constantS_propose_val,
                                     tune_parameter = MH_tune)
  
  
  p_accept <- MetropolisHastings_accept_probability(
    new_numerator_log = propose_log_prior_constantS + propose_log_like_data,
    old_numerator_log = prev_log_prior_constantS + prev_log_like_data,
    correction = correction)
  
  next_val <- ifelse(runif(1) < p_accept,
                     yes = constantS_propose_val, 
                     no = constantS_previous_val)
  
  return(next_val)
}



sample_constantA_MetropolisH <- function(constantS_known,
                                         constantA_previous_val,
                                         beta_known,
                                         light_dat = simulated_tree_data$L,
                                         growth_dat = simulated_tree_data$growth,
                                         MH_tune){
  # 
  
  prev_log_prior_constantA <- calculate_constant_prior_logprob(constantA_previous_val)
  
  prev_log_like_data <- loglike_tree(
    constant_s = constantS_known,
    # the next line is for the parameter we are sampling! 
    constant_a = constantA_previous_val, 
    L = light_dat,
    beta = beta_known,
    simulated_tree_data = growth_dat)
  
  ## propose new!
  constantA_propose_val <- MH_propose_new_val(constantA_previous_val,
                                             sig_tune = MH_tune)
  
  
  propose_log_prior_constantA <- calculate_constant_prior_logprob(constantA_propose_val)
  
  propose_log_like_data <- loglike_tree(
    constant_s = constantS_known,
    # the next line is for the parameter we are sampling! 
    constant_a = constantA_propose_val, 
    L = light_dat,
    beta = beta_known,
    simulated_tree_data = growth_dat)
  
  ## Correction! the special step of MH
  correction <- calculate_correction(new_value = constantA_previous_val,
                                     old_value = constantA_propose_val,
                                     tune_parameter = MH_tune)
  
  
  p_accept <- MetropolisHastings_accept_probability(
    new_numerator_log = propose_log_prior_constantA + propose_log_like_data,
    old_numerator_log = prev_log_prior_constantA + prev_log_like_data,
    correction = correction)
  
  next_val <- ifelse(runif(1) < p_accept,
                     yes = constantA_propose_val, 
                     no = constantA_previous_val)
  
  return(next_val)
}




# proposal function to make new parameter values.
# a needs to be positive or else the whole function becomes negative
# gamma is always positive, a could be zero, gamma has no zero

MH_propose_new_val <- function(old_value, sig_tune) {
  rgamma(1, shape = old_value^2/sig_tune^2, rate = old_value/sig_tune^2)
}

MH_propose_density <- function(jump_to, jump_from, sig_tune) {
  dgamma(jump_to, shape = jump_from^2/sig_tune^2, rate = jump_from/sig_tune^2)
}


# calculate the correction factor
calculate_correction <- function(new_value, old_value, tune_parameter){
  # calculate the correction -- asymmetrical proposal
  new_from_old <- MH_propose_density(jump_to = new_value,
                                     jump_from = old_value,
                                     sig_tune = tune_parameter)
  
  old_from_new <- MH_propose_density(jump_to = old_value,
                                     jump_from = new_value,
                                     sig_tune = tune_parameter)
  
  return(old_from_new / new_from_old)
}

MetropolisHastings_accept_probability <- function(new_numerator_log, 
                                                  old_numerator_log, 
                                                  correction) {
  # working on the LOG scale so we SUBTRACT values (the same as dividing these two things)
  difference_log_scale <- new_numerator_log - old_numerator_log
  # use EXP to undo the LOG -- convert back to probability
  r <- exp(difference_log_scale)
  R <- r * correction
  p_accept <- min(1, R)
  return(p_accept)
}


# we also need a function for the prior
calculate_constant_prior_logprob <- function(x){
  dgamma(x = x, shape = 3^2/2^2, rate = 3/2^2, log = TRUE)
}


sample_constantS_MetropolisH(5, 5, 0.3, MH_tune = 0.1)


sample_constantA_MetropolisH(5, 5, 0.3, MH_tune = 0.1)

## Gibbs part

n_iter <- 5000
gibbs_chain <- list(constantA = numeric(n_iter),
                    constantS = numeric(n_iter),
                    beta     = numeric(n_iter))

# start value
gibbs_chain$constantA[1] <- 500
gibbs_chain$constantS[1] <- 10
gibbs_chain$beta[1] <- 1


for (k in 2:n_iter) {
  
  # constant a
  #this takes the all three previous values
  gibbs_chain$constantA[k] <- sample_constantA_MetropolisH(
    constantA_previous_val = gibbs_chain$constantA[k-1],
    constantS_known = gibbs_chain$constantS[k-1],
    beta_known = gibbs_chain$beta[k - 1],
    MH_tune = 2)
  
  #constant s
  gibbs_chain$constantS[k] <- sample_constantS_MetropolisH(
    constantS_previous_val = gibbs_chain$constantS[k - 1],
    constantA_known = gibbs_chain$constantA[k],
    beta_known = gibbs_chain$beta[k - 1],
    MH_tune = 0.1)
  
  #beta
  gibbs_chain$beta[k] <- sample_beta(y = simulated_tree_data$growth,
                                     alpha = growth(
                                       a = gibbs_chain$constantA[k],
                                       s = gibbs_chain$constantS[k],
                                       L = simulated_tree_data$L
                                     ))
}


data <- do.call(cbind, gibbs_chain)
real_param <- tribble(~parameter, ~truth,
                      "beta", true_beta,
                      "constantS", true_constant_s,
                      "constantA", true_constant_a)

gibbs_chain %>% as_tibble() %>% 
  mutate(iter = 1:n_iter) %>% 
  filter(iter > 300) %>% 
  pivot_longer(-iter, names_to = "parameter", values_to = "value") %>% 
  ggplot(aes(x = iter, y = value)) + 
  geom_line()+
  facet_wrap(~parameter, scales = "free_y") + 
  geom_hline(aes(yintercept = truth),
             data = real_param, lwd = 3, col = "red")


post <- gibbs_chain %>% as_tibble() %>% 
  mutate(iter = 1:n_iter) %>% 
  filter(iter > 300)

postmedian <- lapply(post, median)

sequ <- seq(0, 200, 1)

#plot functions with parameters estimated from gibbs and true values over data

plot(simulated_tree_data$L, simulated_tree_data$growth, col = "gray", xlab = "Light availability", ylab = "Hemlock growth")
lines(growth(postmedian$constantA, postmedian$constantS, sequ), col = "black")
lines(growth(true_constant_a, true_constant_s, sequ), col = "red")
legend("topright", legend = c("Estimated parameters", "True parameters"), col = c("black", "red"), lty = 1, cex = 0.75)


