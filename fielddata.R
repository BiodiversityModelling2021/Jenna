#Thursday class
library(dplyr)
library(ggplot2)


data <- read.csv("fielddata.csv")
data <- data[1:20, 1:2]
data <- rename(data, holes = Number.of.holes.per.radish.leaf)
data$holes <- as.numeric(data$holes)
ggplot(data, aes(holes)) + geom_histogram(bins = 6)

#estimate lambda using MLE

#equation for poisson regression
pos <- function(x, a, b){exp(a+b*x)}

loglikpos <- function(df, lambda){
  ll <- sum(dpois(df, lambda, log = TRUE))
  sum(ll)
}

holes <- data$holes
seq.lamb <- seq(1, 100, by = 0.01)

list.lik <- list()
list.lamb <-list()

for (i in 1:length(seq.lamb)){
  a <- seq.lamb[i]
  lik <- loglikpos(holes, a)
  list.lik[[i]] <- lik
  list.lamb[[i]] <- a
}

lik <- do.call(rbind, list.lik)
lamb <- do.call(rbind, list.lamb)

results <-  cbind(lik, lamb)
results <- as.data.frame(results)
results <-  results %>% rename(ll = V1, lambda = V2)

max(results$ll)

max <- results[which.max(results$ll),]

#have lambda not sure how to predict
hist(holes)
curve(dpois(x, max[1,2]), xlim = c(0, 100), add = TRUE)
#lol
curve(dpois(x, max[1,2])*100, xlim = c(0, 100), add = TRUE)
#aha I think that is it?

m1 <- glm(holes ~ 1, family = poisson(link = identity))
summary(m1)
#the intercept from GLM = the MLE lambda value

#liklihood profile
plot(results$lambda, results$ll)

##simulated annealing


