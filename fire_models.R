# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(deeptime)
library(ggh4x)
library(terra)
library(mgcv)
library(logitnorm)
library(doMC)
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

# vectorized logit_norm mean
logit_norm_vect <- function(mu, sigma) {
  if(length(sigma) == 1) sigma <- rep(sigma, length(mu))
  pmean <- numeric(length(mu))

  for(i in 1:length(mu)) {
    pmean[i] <- momentsLogitnorm(mu[i], sigma[i])["mean"]
  }
  return(pmean)
}

# compute the mean of the logit-normal distribution in parallel.
# Intended to be used with a matrix of logit-means and a point estimate
# for the standard deviation.
logit_norm_mean <- function(mu, sigma, cores = parallel::detectCores()-1) {

  registerDoMC(cores = cores)

  arg_list <- lapply(1:ncol(mu), function(j) {
    r <- cbind(mu[, j], sigma)
    colnames(r) <- c("mu", "sigma")
    return(r)
  })

  # compute moments in parallel
  means_list <- foreach(cc = arg_list) %dopar% {
    logit_norm_vect(cc[, "mu"], cc[, "sigma"])
  }

  return(do.call("cbind", means_list))
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

d$fire <- d$fuego_focal

# year as factor and as continuous variable
year_mid <- round(median(unique(d$anio)))
d$year_num <- d$anio - year_mid
d$year_factor <- factor(as.character(d$anio),
                        levels = unique(d$anio))

d$spei <- d$speiHumedo

# Remove NA in relevant variables
rows_use <- complete.cases(d[, c("fire", "cov_prev_factor", "year_num",
                                 "year_factor", "speiHumedo")])
d <- d[rows_use, ]

# remove burned in previous year to avoid vegetation misclassification
filter_fire <- d$fuego_prev == 0
d <- d[filter_fire, ]


# Create auxiliary variables to remove reference level from
# continuous predictors

veg_mat <- model.matrix(~ cov_prev_factor - 1, d)
year_mat <- veg_mat * d$year_num
spei_mat <- veg_mat * d$spei

colnames(year_mat) <- paste("year", cover_levels, sep = "_")
colnames(spei_mat) <- paste("spei", cover_levels, sep = "_")

dfull <- cbind(d, as.data.frame(year_mat), as.data.frame(spei_mat))

# Burn probability models -------------------------------------------------

# year-only model
m1 <- bam(fire ~ year_num + s(year_factor, bs = "re"),
          family = "binomial", method = "REML", data = d)


# full model, with two parameterizations:
# the first is used to extract coefficients by land cover type, and the
# second, to compute the anova.

m5 <- bam(fire ~ cov_prev_factor - 1 +
                 s(year_factor, bs = "re") +
                 year_Forest + year_Shrubland + year_Grassland +
                 year_Pasture + year_Agriculture +
                 spei_Forest + spei_Shrubland + spei_Grassland +
                 spei_Pasture + spei_Agriculture,
          family = "binomial", method = "REML", data = dfull)

# reparameterize to facilitate predictions and anova
m5b <- bam(fire ~ year_num * cov_prev_factor +
                 spei * cov_prev_factor +
                 s(year_factor, bs = "re"),
          family = "binomial", method = "REML", data = dfull)

# Extract coeff table
sm5 <- summary(m5)
sm5$p.table
sm1 <- summary(m1)
sm1$p.table

anova(m5b)
anova(m1)


# Year effects predictions ------------------------------------------------

# Full model
pdata <- expand.grid(
  cov_prev_factor = cover_levels,
  spei = 0,
  year_num = seq(min(d$year_num), max(d$year_num), length.out = 150),
  year_factor = "2011"
)

pdata$cov_prev_factor <- factor(pdata$cov_prev_factor, cover_levels)
pdata$year <- pdata$year_num + year_mid

b_full <- coef(m5b)
coef_use <- grep("year_factor", names(b_full), invert = T)
b <- b_full[coef_use]

lpmat <- predict(m5b, pdata, type = "lpmatrix")[, coef_use]
V <- vcov(m5b, unconditional = T)[coef_use, coef_use]
sigma_year <- gam.vcomp(m5b)[1, 1]

set.seed(6756)
b_sim <- rmvn(10000, b, V) |> t()

lp_mle <- (lpmat %*% b) |> c()
lp_sim <- lpmat %*% b_sim

# psim <- logit_norm_mean(lp_sim, sigma_year)
# saveRDS(psim, "exports/fire_predictions_psim_year.rds")
psim <- readRDS("exports/fire_predictions_psim_year.rds")
pdata$p_mle <- logit_norm_vect(lp_mle, sigma_year)
pdata$p_lower <- apply(psim, 1, quantile, probs = 0.025)
pdata$p_upper <- apply(psim, 1, quantile, probs = 0.975)

# The same with the marginal model
pdata2 <- expand.grid(
  year_num = seq(min(d$year_num), max(d$year_num), length.out = 150),
  year_factor = "2011"
)
pdata2$year <- pdata2$year_num + year_mid

b_full <- coef(m1)
coef_use <- grep("year_factor", names(b_full), invert = T)
b <- b_full[coef_use]

lpmat <- predict(m1, pdata2, type = "lpmatrix")[, coef_use]
V <- vcov(m1, unconditional = T)[coef_use, coef_use]
sigma_year <- gam.vcomp(m1)[1, 1]

set.seed(6756)
b_sim <- rmvn(10000, b, V) |> t()

lp_mle <- (lpmat %*% b) |> c()
lp_sim <- lpmat %*% b_sim

# psim2 <- logit_norm_mean(lp_sim, sigma_year)
# saveRDS(psim2, "exports/fire_predictions_psim_year_marginal.rds")
psim2 <- readRDS("exports/fire_predictions_psim_year_marginal.rds")
pdata2$p_mle <- logit_norm_vect(lp_mle, sigma_year)
pdata2$p_lower <- apply(psim2, 1, quantile, probs = 0.025)
pdata2$p_upper <- apply(psim2, 1, quantile, probs = 0.975)

pdata2$cov_prev_factor <- "All land cover types"

pboth <- rbind(
  pdata[, c("cov_prev_factor", "year_num", "year", "p_mle", "p_lower", "p_upper")],
  pdata2[, c("cov_prev_factor", "year_num", "year", "p_mle", "p_lower", "p_upper")]
)
pboth$cov_prev_factor <- factor(pboth$cov_prev_factor,
                                levels = c(cover_levels, "All land cover types"))

# Aggregate data
data_agg_year <- aggregate(fire ~ cov_prev_factor + year_num, d, mean)
data_agg_year2 <- aggregate(fire ~ year_num, d, mean)
data_agg_year2$cov_prev_factor <- "All land cover types"

data_agg_year <- rbind(
  data_agg_year, data_agg_year2[, colnames(data_agg_year)]
)
data_agg_year$year <- data_agg_year$year_num + year_mid
data_agg_year$cov_prev_factor <- factor(data_agg_year$cov_prev_factor,
                                        levels = c(cover_levels,
                                                   "All land cover types"))

# Map colors (not used)
colores <- c(
  "#10620d", # forest
  "#9cca88", # shrubland
  "#dec313", # grassland
  "#4492ce", # pasture
  "#bd4343", # agriculture
  "black"    # all
)

colores <- c(
  viridis(5, option = "C", end = 0.7), "black"
)

# Plot
fire_year <-
ggplot(pboth, aes(year, p_mle, ymin = p_lower, ymax = p_upper,
                  color = cov_prev_factor, fill = cov_prev_factor)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(aes(year, fire, color = cov_prev_factor), data = data_agg_year,
             inherit.aes = F, alpha = 0.6) +
  # scale_color_viridis(discrete = TRUE, option = "B", end = 0.8,
  #                     name = "Land cover") +
  # scale_fill_viridis(discrete = TRUE, option = "B", end = 0.8,
  #                    name = "Land cover") +
  scale_fill_manual(values = colores) +
  scale_color_manual(values = colores) +

  facet_wrap(vars(cov_prev_factor), nrow = 2,
             axes = "all") +
             # axes = "margins") +
  ylab("Annual burn probability (%)") +
  xlab("Year") +
  scale_y_continuous(labels = scales::label_percent(suffix = ""),
                     expand = c(0, 0), limits = c(-0.005, 0.18)) +
  nice_theme() +
  theme(legend.position = "none",
        panel.spacing.x = unit(2, "mm"),
        panel.spacing.y = unit(4, "mm"),
        axis.text = element_text(size = 8)) +
  labs(tag = "A")
fire_year

# ggsave("figures/05a) fire year.png",
#        width = 16, height = 12, units = "cm")

# SPEI effects predictions ------------------------------------------------

pdata <- expand.grid(
  cov_prev_factor = cover_levels,
  spei = seq(min(d$spei), max(d$spei), length.out = 150),
  year_num = 0,
  year_factor = "2011"
)

pdata$cov_prev_factor <- factor(pdata$cov_prev_factor, cover_levels)

b_full <- coef(m5b)
coef_use <- grep("year_factor", names(b_full), invert = T)
b <- b_full[coef_use]

lpmat <- predict(m5b, pdata, type = "lpmatrix")[, coef_use]
V <- vcov(m5b, unconditional = T)[coef_use, coef_use]
sigma_year <- gam.vcomp(m5b)[1, 1]

set.seed(6756)
b_sim <- rmvn(10000, b, V) |> t()

lp_mle <- (lpmat %*% b) |> c()
lp_sim <- lpmat %*% b_sim

# psim <- logit_norm_mean(lp_sim, sigma_year)
# saveRDS(psim, "exports/fire_predictions_psim_spei.rds")
psim <- readRDS("exports/fire_predictions_psim_spei.rds")
pdata$p_mle <- logit_norm_vect(lp_mle, sigma_year)
pdata$p_lower <- apply(psim, 1, quantile, probs = 0.025)
pdata$p_upper <- apply(psim, 1, quantile, probs = 0.975)


# Aggregate data
data_agg_spei <- do.call("rbind", lapply(cover_levels, function(cc) {
  dsub <- d[d$cov_prev_factor == cc, ]
  dsub <- dsub[order(dsub$spei), ]
  gg <- rep(letters[1:15], each = ceiling(nrow(dsub) / 15))
  dsub$group <- gg[1:nrow(dsub)]

  dagg <- aggregate(cbind(fire, spei) ~ group, dsub, mean)
  dagg$cov_prev_factor <- cc
  return(dagg)
}))

data_agg_spei$cov_prev_factor <- factor(data_agg_spei$cov_prev_factor,
                                        levels = cover_levels)

# Plot
fire_spei <-
ggplot(pdata, aes(spei, p_mle, ymin = p_lower, ymax = p_upper,
                  color = cov_prev_factor, fill = cov_prev_factor)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(aes(spei, fire, color = cov_prev_factor), data = data_agg_spei,
             inherit.aes = F, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = "C", end = 0.7,
                      name = "Land cover") +
  scale_fill_viridis(discrete = TRUE, option = "C", end = 0.7,
                     name = "Land cover") +
  facet_wrap(vars(cov_prev_factor), nrow = 2,
             axes = "all") +
  ylab("Annual burn probability (%)") +
  xlab("Standardized Precipitation Evapotranspiration Index (SPEI)") +
  scale_y_continuous(labels = scales::label_percent(suffix = ""),
                     expand = c(0.001, 0.001), limits = c(0, 0.15))+#max(pdata$p_upper))) +
  nice_theme() +
  theme(legend.position = "none",
        panel.spacing.x = unit(2, "mm"),
        panel.spacing.y = unit(4, "mm"),
        axis.text = element_text(size = 8)) +
  labs(tag = "B")
fire_spei

# ggsave("figures/05b) fire spei.png",
#        width = 16, height = 12, units = "cm")



# Merge plots -------------------------------------------------------------

pp <- deeptime::ggarrange2(
  fire_spei, fire_year, nrow = 2,
  labels = c("A", "B"),
  label.args = list(gp = grid::gpar(face = "plain", cex = 1.1))
)

ggsave("figures/04) fire-year-spei.png", plot = pp,
       width = 16, height = 20, units = "cm")