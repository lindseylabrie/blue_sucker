# load packages
library(tidyverse)
library(rethinking)
library(brms)
library(janitor)
library(readr)

blue_sucker_2021_data <- read_csv("~/Desktop/blue_sucker/data/blue_sucker_2021_data.csv")%>% 
  clean_names()%>% 
  filter(!is.na(lab_sex)) %>% 
  filter(!is.na(gonad_weight_g)) %>% 
  select(fish_id, length_mm, weight_g, lab_sex, gonad_weight_g, ovary_1_g, ovary_2_g, 
         o1_sample_1, o1_sample_2, o1_sample_3, o2_sample_1, o2_sample_2, o2_sample_3) %>% #selecting which columns to pay attention to
  mutate(o1_average_count = as.integer((o1_sample_1+o1_sample_2+o1_sample_3)/3),
         o2_average_count = as.integer((o2_sample_1+o2_sample_2+o2_sample_3)/3)) %>% # I know there's a better way. Getting the mean egg counts.
  mutate(ovary_1_total_eggs = as.integer(ovary_1_g*o1_average_count),
         ovary_2_total_eggs = as.integer(ovary_2_g*o2_average_count),
         combined_egg_total = as.integer(ovary_1_total_eggs+ovary_2_total_eggs)) %>%  # multiplying weight of ovary times average egg count, then adding them together to get total eggs in each fish
  mutate(egg_total_simplified = combined_egg_total/10000,
         length_c = length_mm - mean(length_mm),
         length_s = length_c/sd(length_mm),
         weight_s = (weight_g - mean(weight_g))/sd(weight_g))
  #simplifying the total eggs so the computer doesn't have to make such large calculations
  #centering the length, standardizing the length, standardizing the weight

d <- blue_sucker_2021_data #simplifying what it's called.
  
d %>% 
ggplot(aes(x=length_mm, y=weight_g, color=lab_sex)) +
  geom_point() +
  geom_smooth()
# a sexy graph of lengths, weights, and sexes


# log(mu) = alpha + b_w*x_w + b_l*x_l + b_wl*(x_w*x_l) (underscore denotes subscript)

# need priors for alpha, beta_w, beta_l, beta_wl

N = 1000  # number of simulations (change as needed)

# simulate priors
priors <- tibble(a = rnorm(N, 0.1, 0.5),
                 bw = rnorm(N, 0, 0.1),
                 bl = rnorm(N, 0, 0.1),
                 bwl = rnorm(N, 0, 0.5),
                 shape = rexp(N, 0.1), # if a poisson, no shape here
                 sim = 1:N)

priors

# data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
# x <- data_to_plot$length - mean(data_to_plot$length)
length_s <- d$length_s
weight_s <- d$weight_s

# combine and simulate
prior_and_x <- priors %>% expand_grid(length_s = length_s) %>%    # combine priors and x's
  left_join(d %>% distinct(length_s, weight_s)) %>% 
  distinct(sim, weight_s, length_s, a, bw, bl, bwl, shape) %>% 
  mutate(mu = exp(a + bw*weight_s + bl*length_s+ bwl*length_s*weight_s),                              # simulate regressions
         y = rgamma(nrow(.), scale = mu/shape, shape = shape))            # simulate data (e.g., y_rep)

library(scales)

# plot
prior_and_x %>% 
  ggplot(aes(x = weight_s, y = mu, group = sim)) + 
  geom_line() +
  # geom_point(aes(y = y)) +
  labs(y = "sim")+
  scale_y_log10(labels = comma) +
  NULL

# What is fecundity?
## Fecundity is the number of eggs a female fish will lay in a spawning season.

# What I'm trying to do in this scenario: use length and weight as predictors for egg count (fecundity)

# What else could be done? 

## Use only length as a predictor for egg count.
## Use only weight as a predictor for egg count.
## See if total length would be a good predictor of male gonad weight.
## See if total weight would be a good predictor of male gonad weight.
## See if length would be a predictor of GSI (Gonadosomatic Index)?: GSI = (gonad weight/wet weight)*100
## See if total weight be a predictor of GSI? 
###GSI could be used for both males and females.

get_prior(egg_total_simplified ~ length_s + weight_s + length_s*weight_s, 
           data = d,
           family = Gamma(link="log"))


bsr_brm <- brm(combined_egg_total ~ length_s + weight_s + length_s*weight_s, 
                data = d,
                family = negbinomial(link="log"),
                prior = c(prior(normal(11, 1), class = "Intercept"),
                          prior(normal(0, 2), class = "b", coef="length_s"),
                          prior(normal(0, 2), class = "b", coef="weight_s"),
                          prior(normal(0,0.5), class ="b", coef="length_s:weight_s")), 
                cores = 1, chains = 1, iter = 1000,
                sample_prior = "yes")

bsr_brm
plot(conditional_effects(bsr_brm), points = T)

################### more graphs #########################

####### Use only weight as a predictor for egg count #######

N = 1000  # number of simulations (change as needed)

# simulate priors
priors1 <- tibble(a = rnorm(N, 11.5, 0.5),
                 bw = rnorm(N, 0, 0.5),
                 phi = rexp(N, 0.001),
                 sim = 1:N)

prior_and_x_weight <- priors1 %>% expand_grid(weight_s = weight_s) %>%    # combine priors and x's
  # left_join(d %>% distinct(length_s, weight_s)) %>% 
  distinct(sim, weight_s, a, bw, phi) %>% 
  mutate(mu = exp(a + bw*weight_s),                              # simulate regressions
         y = rnbinom(nrow(.), size = phi, mu = mu))

# plot
prior_and_x_weight %>% 
  ggplot(aes(x = weight_s, y = mu, group = sim)) + 
  geom_line(alpha=0.1) +
  geom_point(aes(y = y)) +
  labs(y = "sim")+
  scale_y_log10(labels = comma) +
  geom_hline(yintercept = c(250000,
                            50000,
                            100000), color = "red") +
  NULL

get_prior(combined_egg_total ~ weight_s, 
          data = d,
          family = negbinomial(link="log"))

bsr_brm_weight <- brm(combined_egg_total ~ weight_s, 
               data = d,
               family = negbinomial(link="log"),
               prior = c(prior(normal(11.5, 0.5), class = "Intercept"),
                         prior(normal(0, 0.5), class = "b", coef="weight_s"),
                         prior(exponential(0.001), class = "shape")), 
               cores = 1, chains = 1, iter = 1000,
               file = "models/bsr_brm_weight.rds",
               file_refit = "on_change",
               sample_prior = "yes")

plot(conditional_effects(bsr_brm_weight), points = T)


####### Use only length as a predictor for egg count ########

N = 1000  # number of simulations (change as needed)

# simulate priors
priors2 <- tibble(a = rnorm(N, 11.5, 0.3),
                  bl = rnorm(N, 0, 0.4),
                  phi = rexp(N, 0.001),
                  sim = 1:N)

prior_and_x_length <- priors2 %>% expand_grid(length_s = length_s) %>%    # combine priors and x's
  # left_join(d %>% distinct(length_s, weight_s)) %>% 
  distinct(sim, length_s, a, bl, phi) %>% 
  mutate(mu = exp(a + bl*length_s),                              # simulate regressions
         y = rnbinom(nrow(.), size = phi, mu = mu))

# plot
prior_and_x_length %>% 
  ggplot(aes(x = length_s, y = mu, group = sim)) + 
  geom_line(alpha=0.1) +
  geom_point(aes(y = y)) +
  labs(y = "sim")+
  scale_y_log10(labels = comma) +
  geom_hline(yintercept = c(250000,
                            50000,
                            100000), color = "red") +
  NULL

bsr_brm_length <- brm(combined_egg_total ~ length_s, 
                      data = d,
                      family = negbinomial(link="log"),
                      prior = c(prior(normal(11.5, 0.3), class = "Intercept"),
                                prior(normal(0, 0.4), class = "b", coef="length_s"),
                                prior(exponential(0.001), class = "shape")), 
                      cores = 1, chains = 1, iter = 1000,
                      file = "models/bsr_brm_length.rds",
                      file_refit = "on_change",
                      sample_prior = "yes")

plot(conditional_effects(bsr_brm_length), points = T)


### Danielle's model ###

bsr_brm_length_danielle <- brm(combined_egg_total ~ length_s, 
                      data = d,
                      family = Gamma(link="log"),
                      prior = c(prior(normal(12, 0.2), class = "Intercept"),
                                prior(normal(0, 0.1), class = "b", coef="length_s"),
                                prior(exponential(0.0005), class = "shape")), 
                      cores = 1, chains = 1, iter = 1000,
                      file = "models/bsr_brm_length_danielle.rds",
                      file_refit = "on_change",
                      sample_prior = "yes")

plot(conditional_effects(bsr_brm_length_danielle), points = T)


### Lauren's model ###

bsr_brm3 <- brm(egg_total_simplified ~ length_s + weight_s + length_s*weight_s, 
                 data = d,
                 family = Gamma(link="log"),
                 prior = c(prior(normal(10, 5), class = "Intercept"),
                           #prior(normal(0, 0.1), class = "b", coef="length_s"),
                           #prior(normal(0, 0.1), class = "b", coef="length_s:weight_s"),
                           prior(normal(20, 15), class = "b", coef="weight_s"),#is this what i do when i want to use different priors for different betas?
                           prior(exponential(0.001), class = "shape")), 
                 cores = 4, chains = 2, iter = 1000,
                 sample_prior = "yes")
plot(gaus_brm3)
gaus_brm3
plot(conditional_effects(gaus_brm3), points = T)





# length and weight are highly correlated, which is why the length model looks so wonky.

posteriors_bsr_brm <- as_draws_df(bsr_brm) %>% as_tibble() %>% expand_grid(length_s = length_s) %>%    # combine priors and x's
  left_join(d %>% distinct(length_s, weight_s)) %>% 
  mutate(mu = exp(b_Intercept + b_weight_s*0 + b_length_s*length_s + 
                    `b_length_s:weight_s`*length_s*0),
         mu_prior = exp(prior_Intercept + prior_b_weight_s*0 + prior_b_length_s*length_s + 
                          `prior_b_length_s:weight_s`*length_s*0)
         )

ggplot(posteriors_bsr_brm %>% filter(.draw <=1000), aes(x=length_s, y=mu_prior)) + 
  geom_point() +
  geom_line(aes(group  = .draw), alpha = 0.1)



ggplot(d, aes(x = length_mm, y = weight_g)) + 
  geom_point()

