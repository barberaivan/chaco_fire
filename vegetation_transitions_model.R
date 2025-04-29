# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(deeptime)
library(ggh4x)
library(terra)
library(mgcv)
theme_set(theme_classic())

# Functions ---------------------------------------------------------------

softmax <- function(x) exp(x) / sum(exp(x)) # inverse multinomial-logit link

normalize <- function(x) x / sum(x)

# to compute mean and 95 % CI from samples.
mean_ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  tmp <- c(mean(x), qs) %>% unname
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

# just the ci
ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  names(qs) <- paste(name, c("lower", "upper"), sep = "")
  return(qs)
}

r2bern <- function(p) {
  var_mu <- var(p)
  var_y <- mean(p * (1 - p))
  return(var_mu / (var_mu + var_y))
}

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.35),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# Data --------------------------------------------------------------------

d <- read.csv("data/R5_muestreoFinal_6meses.csv")
# head(d)

cover_levels <- c("Forest", "Shrubland", "Grassland", "Pasture", "Agriculture")
K <- length(cover_levels)
k <- K-1

d$cov_prev_factor <- plyr::revalue(as.character(d$cob_prev),
                                   c("1" = cover_levels[1],
                                     "2" = cover_levels[2],
                                     "4" = cover_levels[3],
                                     "5" = cover_levels[4],
                                     "6" = cover_levels[5]))
d$cov_prev_factor <- factor(d$cov_prev_factor,
                            levels = cover_levels)

d$cov_next_factor <- plyr::revalue(as.character(d$cob_next),
                                   c("1" = cover_levels[1],
                                     "2" = cover_levels[2],
                                     "4" = cover_levels[3],
                                     "5" = cover_levels[4],
                                     "6" = cover_levels[5]))
d$cov_next_factor <- factor(d$cov_next_factor,
                            levels = cover_levels)

d$cov_next_num <- as.numeric(d$cov_next_factor) - 1 # ref cat starts at zero in gam


fire_levels <- c("Burned", "Unburned")
d$fire <- plyr::revalue(as.character(d$fuego_focal),
                        c("0" = fire_levels[2],
                          "1" = fire_levels[1]))
d$fire <- factor(d$fire, levels = fire_levels)


# Remove NA in relevant variables
rows_use <- complete.cases(d[, c("fire", "cov_prev_factor", "cov_next_num",
                                 "fuego_prev", "fuego_focal", "fuego_next")])
d <- d[rows_use, ]

# remove burned in previous or next year to avoid vegetation misclassification
filter_fire <- d$fuego_next == 0 & d$fuego_prev == 0
d <- d[filter_fire, ]

# Transition model as a function of fire ----------------------------------

vegmod1 <- gam(
  list(cov_next_num
       ~ cov_prev_factor * fire,
       ~ cov_prev_factor * fire,
       ~ cov_prev_factor * fire,
       ~ cov_prev_factor * fire),
  family = multinom(K = k), data = d
)

anova(vegmod1)


# Predictions (simulations to compute CI)
pdata <- expand.grid(
  cov_prev_factor = cover_levels,
  fire = fire_levels
)

nsim <- 10000

# make linear predictor matrix:
lpmat <- predict(vegmod1, pdata, type = "lpmatrix")
lpmat <- lpmat[, 1:(ncol(lpmat) / k)] # all matrices are the same

np <- nrow(pdata)

# create an array for simulated predicted probabilites.
pred_arr <- array(NA, dim = c(np, K, nsim),   # it's an array with nsim tpm
                  dimnames = list(
                    case = 1:np,
                    cov_next_factor = cover_levels,
                    sim = 1:nsim
                  ))

# the model parameter estimates, under hypothetical infinite replicates of the
# study, are assumed to follow a multivariate normal distribution, with the
# MLE as their mean and a certain variance-covariance matrix.
# We simulate parameter vectors from this distribution:
set.seed(123)
coef_samples <- rmvn(nsim, coef(vegmod1),
                     V = vcov(vegmod1, unconditional = F, freq = T)) %>% t

# each linear predictor has the same number of coefficients, so
coef_ids <- rep(1:k, each = nrow(coef_samples) / k)

for(i in 1:nsim) {
  # i = 1
  # order the coefficients vector in matrix form. It enters by column, as the
  # default for matrix()
  # i = 1 # just to test the loop
  coef_temp <- coef_samples[, i]
  # turn potential Inf into something reasonable (exp(700) is finite, but a bit
  # above return Inf and breaks the loop)
  coef_temp[coef_temp > 700] <- 700

  coef_mat <- matrix(coef_temp, ncol = k)
  # add reference linear predictor
  linpred <- cbind(rep(0, np), lpmat %*% coef_mat)
  pred_arr[, , i] <- apply(linpred, 1, softmax) %>% t
}

# compute ci and longanize array:
pred_long_ci <- as.data.frame.table(apply(pred_arr, 1:2, ci))

# get the MLE in a similar way, but before match the names.
prob_hat <- predict(vegmod1, pdata, type = "response")
dimnames(prob_hat) <- dimnames(pred_arr)[1:2]
pred_long_mle <- as.data.frame.table(prob_hat)
pred_long_mle <- cbind(data.frame("Var1" = rep("p_mle", np)),
                       pred_long_mle)

# merge with ci data_veg:
pred_long <- rbind(pred_long_ci, pred_long_mle)
# a bit wider for ggplot:
pred <- pivot_wider(pred_long, names_from = "Var1", values_from = "Freq")

# merge with predictor variables
pdata$case <- 1:np
pred$case <- as.numeric(as.character(pred$case))
pred_veg <- left_join(pred, pdata, by = "case")

pred_veg$cov_next_factor <- factor(pred_veg$cov_next_factor,
                                   levels = cover_levels)
pred_veg$cov_prev_factor <- factor(pred_veg$cov_prev_factor,
                                   levels = cover_levels)

pred_veg$column_name <- "Next land cover"
pred_veg$row_name <- "Previous land cover"


# Transition probability quotients (burned / unburned)

prob_hat <- predict(vegmod1, pdata, type = "response")

prob_q <- prob_hat[pdata$fire == "Burned", ] / prob_hat[pdata$fire == "Unburned", ]
dimnames(prob_q) <- list(
  cov_prev_factor = cover_levels,
  cov_next_factor = cover_levels
)

# compute p-value. First, get quotient distribution
prob_q_arr <- pred_arr[pdata$fire == "Burned", , ] / pred_arr[pdata$fire == "Unburned", , ]
dimnames(prob_q_arr) <- dimnames(prob_q)

pfun <- function(x) {
  cdf <- ecdf(x)
  pp <- cdf(1)
  if(pp < 0.5) return(pp / 2) else return((1 - pp) / 2)
}

pval_mat <- apply(prob_q_arr, 1:2, pfun)

# Longanize to plot
probq_long <- as.data.frame.table(prob_q)
names(probq_long)[3] <- "q"
pval_long <- as.data.frame.table(pval_mat)
names(pval_long)[3] <- "pval"

qp <- cbind(data.frame(y = 0.3, x = 1),
            probq_long, pval = pval_long$pval)

qp$qtext <- format(round(qp$q, 3), nsmall = 3) |> as.character()
qp$ptext <- format(round(qp$pval, 3), nsmall = 3) |> as.character()
qp$ptext[qp$ptext == "0.000"] <- "<0.001"
qp$ptext <- paste("(", qp$ptext, ")", sep = "")
qp$text <- paste(qp$qtext, qp$ptext, sep = "\n")

qp$column_name <- "Next land cover"
qp$row_name <- "Previous land cover"


# Plot everything

ggplot(pred_veg, aes(cov_next_factor, p_mle, ymin = p_lower, ymax = p_upper,
                     fill = fire)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.6,
           position = "dodge") +
  geom_linerange(position = position_dodge(width = 0.9)) +

  geom_text(data = qp, mapping = aes(x = x, y = y, label = text),
            inherit.aes = F, size = 5.5/.pt) +

  facet_nested(rows = vars(row_name, cov_prev_factor),
               cols = vars(column_name, cov_next_factor),
               scales = "free_x",
               axes = "x",
               nest_line = element_line(colour = "gray50",
                                        linewidth = 0.2)) +
  ylab("Transition probability (%)") +
  xlab("Fire") +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005),
                     labels = scales::label_percent(suffix = "")) +
  scale_fill_viridis(discrete = T, option = "D", end = 0.5,
                     name = "Fire", direction = 1) +
  nice_theme() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.y = unit(4, "mm"),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("figures/03) vegetation-transitions-fire.png",
       width = 14, height = 15, units = "cm")

# Assess N
table(d$cov_prev_factor, d$fire)
sum(d$fuego_focal) / nrow(d)