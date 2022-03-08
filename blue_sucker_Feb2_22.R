# load packages
library(tidyverse)
library(rethinking)
library(brms)
library(janitor)
library(readr)
library(tidybayes)

blue_sucker_2021_data <- read_csv("data/blue_sucker_2021_data.csv")%>% 
  clean_names()%>% 
  filter(!is.na(lab_sex)) %>% 
  filter(!is.na(gonad_weight_g)) %>% 
  select(fish_id, length_mm, weight_g, lab_sex, gonad_weight_g, ovary_1_g, ovary_2_g, 
         o1_sample_1, o1_sample_2, o1_sample_3, o2_sample_1, o2_sample_2, o2_sample_3) %>% #selecting which columns to pay attention to
  mutate(o1_average_count = as.integer((o1_sample_1+o1_sample_2+o1_sample_3)/3),
         o2_average_count = as.integer((o2_sample_1+o2_sample_2+o2_sample_3)/3), # I know there's a better way. Getting the mean egg counts.
         ovary_1_total_eggs = as.integer(ovary_1_g*o1_average_count),
         ovary_2_total_eggs = as.integer(ovary_2_g*o2_average_count),
         combined_egg_total = as.integer(ovary_1_total_eggs+ovary_2_total_eggs),  # multiplying weight of ovary times average egg count, then adding them together to get total eggs in each fish
         egg_total_simplified = combined_egg_total/10000,
         length_c = length_mm - mean(length_mm),
         length_s = length_c/sd(length_mm),
         weight_s = (weight_g - mean(weight_g))/sd(weight_g),
         gsi = (gonad_weight_g/weight_g)*100,
         length_s_squared = length_s*length_s) 
  #GSI = (gonad weight/wet weight)*100
  #simplifying the total eggs so the computer doesn't have to make such large calculations
  #centering the length, standardizing the length, standardizing the weight

d <- blue_sucker_2021_data #simplifying what it's called.
  
LengthsWeights <- d %>% 
ggplot(aes(x=length_mm, y=weight_g, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="Lengths and Weights by Sex, 2021",
       x="Length (mm)",
       y= "Wet Weight (g)")

ggsave(LengthsWeights, file = "plots/LengthsWeighs.png", dpi=750,  width = 5, height = 3,
       units = "in")


WeightGSI <- d %>% 
  ggplot(aes(x=weight_g, y=gsi, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="GSI by Sex and Weight (g), 2021",
       x="Wet Weight (g)",
       y= "GSI")

ggsave(WeightGSI, file = "plots/WeightGSI.png", dpi=750,  width = 5, height = 3,
       units = "in")

LengthGSI <- d %>% 
  ggplot(aes(x=length_mm, y=gsi, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="GSI by Sex and Length (mm), 2021",
       x="Length (mm)",
       y= "GSI")

ggsave(LengthGSI, file = "plots/LengthGSI.png", dpi=750,  width = 5, height = 3,
       units = "in")

StandardLengthGSI <- d %>% 
  ggplot(aes(x=length_s, y=gsi, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="GSI by Sex and Length (standardized)",
       x="Length (standardized)",
       y= "GSI")
#same as the regular length plot

ggsave(StandardLengthGSI, file = "plots/StandardLengthGSI.png", dpi=750,  width = 5, height = 3,
       units = "in")

StandardLengthEggTotal <- d %>% 
  ggplot(aes(x=length_s, y=combined_egg_total)) +
  geom_point() +
  geom_smooth()+
  labs(title="Lengths and estimated total egg counts",
       x="Length (standardized)",
       y= "Estimated egg total")

ggsave(StandardLengthEggTotal, file = "plots/StandardLengthEggTotal.png", dpi=750,  width = 5, height = 3,
       units = "in")


# length as predictor of gsi
get_prior(gsi ~ length_s + length_s*lab_sex + (1|fish_id), 
          data = d,
          family = negbinomial(link="log"))

length_gaus <- brm(gsi ~ length_s + length_s*lab_sex + (1|fish_id), 
                           data = d,
                           family = gaussian(),
                           cores = 1, chains = 4, iter = 5000,
                           sample_prior = "yes",
                           file="models/4_chain_length_gaus.rds",
                           file_refit = "on_change")


plot(conditional_effects(length_gaus), points = T)

pp_check(length_gaus)
pp_check(length_gaus, type="stat")
pp_check(length_gaus, type="stat_grouped", group="lab_sex")

bayes_R2(length_gaus)

#weight as a predictor of gsi
weight_gaus <- brm(gsi ~ weight_s + weight_s*lab_sex + (1|fish_id), 
                   data = d,
                   family = gaussian(),
                   cores = 1, chains = 4, iter = 5000,
                   sample_prior = "yes",
                   file="models/4_chains_weight_gaus.rds",
                   file_refit = "on_change")

plot(conditional_effects(weight_gaus), points = T)

pp_check(weight_gaus)
pp_check(weight_gaus, type="stat")
pp_check(weight_gaus, type="stat_grouped", group="lab_sex")

bayes_R2(weight_gaus)


# What is fecundity?
## Fecundity is the number of eggs a female fish will lay in a spawning season.

# What I'm trying to do in this scenario: use length and weight as predictors for egg count (fecundity)

# What else could be done? 

## Use only length as a predictor for egg count.
## Use only weight as a predictor for egg count.
## See if total length would be a good predictor of male gonad weight.
## See if length weight be a good predictor of male gonad weight.
## See if length be a predictor of GSI (Gonadosomatic Index)?: GSI = (gonad weight/wet weight)*100
## Could total weight be a predictor of GSI? 
###GSI could be used for both males and females.

# Using length as a predictor

get_prior(egg_total_simplified ~ length_s + I(length_s^2) + (1|fish_id), 
           data = d,
           family = negbinomial(link="log"))


length_bsr_negbinom <- brm(combined_egg_total ~ length_s + I(length_s^2) + (1|fish_id), 
                           data = d,
                           family = negbinomial(link="log"),
                           prior = c(prior(normal(13, 1), class = "Intercept"),
                                      prior(normal(-2, 3), class = "b", coef="length_s")),
                           cores = 1, chains = 4, iter = 1000,
                           sample_prior = "yes",
                           file="models/length_bsr_negbinom.rds",
                           file_refit = "on_change")

plot(conditional_effects(length_bsr_negbinom, re_formula = NULL), points = T)

length_bsr_negbinom
plot(conditional_effects(length_bsr_negbinom), points = T)

pp_check(length_bsr_negbinom)
pp_check(length_bsr_negbinom, type = "hist")

saveRDS(length_bsr_negbinom, "models/length_bsr_negbinom.rds")

as_draws_df(length_bsr_negbinom)

### Trying with weight_s

weight_bsr_negbinom <- brm(combined_egg_total ~ weight_s + I(weight_s^2) + (1|fish_id), 
                           data = d,
                           family = negbinomial(link="log"),
                           prior = c(prior(normal(13, 5), class = "Intercept"),
                                     prior(normal(0, 1), class = "b", coef="weight_s")),
                           cores = 1, chains = 4, iter = 1000,
                           sample_prior = "yes",
                           file="models/weight_bsr_negbinom.rds",
                           file_refit = "on_change")
weight_bsr_negbinom
plot(conditional_effects(weight_bsr_negbinom, re_formula=NULL), points = T)

pp_check(weight_bsr_negbinom)

cond_effect_weight <- conditional_effects(weight_bsr_negbinom)
cond_effect_weight$weight_s

cond_effect_weight$weight_s %>% 
  ggplot(aes(x=weight_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = weight_bsr_negbinom$data, aes(x=weight_s, y=combined_egg_total))+
  theme_default()

cond_data_weight <- weight_bsr_negbinom$data %>% distinct(weight_s, combined_egg_total)

posts_weight <- add_epred_draws(weight_bsr_negbinom, newdata= weight_bsr_negbinom$data %>% 
                                  distinct(weight_s, fish_id), re_formula = NULL)

d_weight <- d %>% distinct(weight_g, weight_s)

PosteriorWeight <- posts_weight %>%
  group_by(weight_s) %>% 
  left_join(d_weight) %>% 
  median_qi(.epred) %>% 
  ggplot(aes(x = weight_s, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = combined_egg_total)) +
  labs(title= "Blue Sucker Fecundity Prediction",
       x="Weight (standardized)",
       y="Predicted total egg count")

# try figuring out how to fix it so it's not so bumpy

ggsave(PosteriorAll, file = "plots/PosteriorAll.png", dpi = 750, width = 7, height = 5,
       units = "in")








# GSI gaus?
gsi_length_negbinom <- brm(gsi ~ length_s + length_s*lab_sex + I(length_s^2)+ (1|fish_id), 
                   data = d,
                   family = gaussian(),
                   cores = 1, chains = 1, iter = 5000,
                   sample_prior = "yes")
# ,
#                    file="models/gsi_length_gaus.rds",
#                    file_refit = "on_change")

plot(conditional_effects(gsi_length_gaus), points = T)





















# #simplified egg count
# eggs_gaus <- brm(egg_total_simplified ~ length_s + (1|fish_id), 
#                    data = d,
#                    family = gaussian(),
#                    cores = 1, chains = 1, iter = 1000,
#                    sample_prior = "yes",
#                    file="models/eggs_gaus.rds",
#                    file_refit = "on_change")
# 
# plot(conditional_effects(eggs_gaus), points = T)
# 
# ## JUST using LENGTH as a predictor for egg count
# 
# # log(mu) = alpha + b_l*x_l (underscore denotes subscript)
# 
# # need priors for alpha, beta_l
# 
# N = 1000  # number of simulations (change as needed)
# 
# # simulate priors
# priors <- tibble(a = rnorm(N, 0, 0.5),
#                  bl = rnorm(N, 0, 0.1),
#                  shape = rexp(N, 0.1), # if a poisson, no shape here
#                  sim = 1:N)
# 
# priors
# 
# # data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
# # x <- data_to_plot$length - mean(data_to_plot$length)
# length_mm <- d$length_mm
# 
# # combine and simulate
# prior_and_x <- priors %>% expand_grid(length_s = length_s) %>%    # combine priors and x's
#   left_join(d %>% distinct(length_s)) %>% 
#   distinct(sim, length_s, a, bl, shape) %>% 
#   mutate(mu = exp(a  + bl*length_s),                              # simulate regressions
#          y = rgamma(nrow(.), scale = mu/shape, shape = shape))            # simulate data (e.g., y_rep)
# 
# library(scales)
# 
# # plot
# prior_and_x %>% 
#   ggplot(aes(x = length_s, y = mu, group = sim)) + 
#   geom_line() +
#   # geom_point(aes(y = y)) +
#   labs(y = "sim")+
#   scale_y_log10(labels = comma) +
#   NULL



# # length and weight are highly correlated, which is why the length model looks so wonky.
# 
# posteriors_bsr_brm <- as_draws_df(bsr_brm) %>% as_tibble() %>% expand_grid(length_s = length_s) %>%    # combine priors and x's
#   left_join(d %>% distinct(length_s, weight_s)) %>% 
#   mutate(mu = exp(b_Intercept + b_weight_s*0 + b_length_s*length_s + 
#                     `b_length_s:weight_s`*length_s*0),
#          mu_prior = exp(prior_Intercept + prior_b_weight_s*0 + prior_b_length_s*length_s + 
#                           `prior_b_length_s:weight_s`*length_s*0)
#          )
# 
# ggplot(posteriors_bsr_brm %>% filter(.draw <=1000), aes(x=length_s, y=mu_prior)) + 
#   geom_point() +
#   geom_line(aes(group  = .draw), alpha = 0.1)
# 
# 
# 
# ggplot(d, aes(x = length_mm, y = weight_g)) + 
#   geom_point()
# 
# 
# 
# ############ Working in class ###################
# 
# # Simulating Weight and Length #
# 
# N = 1000  # number of simulations (change as needed)
# 
# # simulate priors
# priors <- tibble(a = rnorm(N, 11, 1),
#                  bw = rnorm(N, 0, 2),
#                  bl = rnorm(N, 0, 2),
#                  shape = rexp(N, 0.1), # if a poisson, no shape here
#                  sim = 1:N)
# 
# priors
# 
# # data (only the x values, since we're simulating y and mu and pretending we don't have them yet)
# # x <- data_to_plot$length - mean(data_to_plot$length)
# length_mm <- d$length_mm
# weight_g <- d$weight_g
# 
# # combine and simulate
# prior_and_x <- priors %>% expand_grid(length_mm = length_mm) %>%    # combine priors and x's
#   left_join(d %>% distinct(length_mm, weight_g)) %>% 
#   distinct(sim, weight_g, length_mm, a, bw, bl, shape) %>% 
#   mutate(mu = exp(a + bw*weight_g + bl*length_mm),                              # simulate regressions
#          y = rgamma(nrow(.), scale = mu/shape, shape = shape))            # simulate data (e.g., y_rep)
# 
# library(scales)
# 
# # plot
# prior_and_x %>% 
#   ggplot(aes(x = weight_g, y = mu, group = sim)) + 
#   geom_line() +
#   # geom_point(aes(y = y)) +
#   labs(y = "sim")+
#   scale_y_log10(labels = comma) 
# 
# get_prior(combined_egg_total ~ length_mm + weight_g, 
#           data = d,
#           family = Gamma(link="log"))
# 
# bsr_brm <- brm(combined_egg_total ~ length_mm + weight_g, 
#                data = d,
#                family = negbinomial(link="log"),
#                prior = c(prior(normal(13, 1), class = "Intercept"),
#                          prior(normal(0, 2), class = "b", coef="length_mm"),
#                          prior(normal(0, 2), class = "b", coef="weight_g")),
#                cores = 1, chains = 1, iter = 1000,
#                sample_prior = "yes",
#                file="models/bsr_brm_test.rds",
#                refit = "on_change")
# 
# bsr_brm
# plot(conditional_effects(bsr_brm), points = T)
# 
# saveRDS(bsr_brm, "models/bsr_brm_test.rds")
# dir.create("models")
# 
# 
# bsr_brm <- brm(combined_egg_total ~ length_s + length_s*length_s, 
#                data = d,
#                family = negbinomial(link="log"),
#                prior = c(prior(normal(13, 1), class = "Intercept"),
#                          prior(normal(0, 3), class = "b", coef="length_s")),
#                cores = 1, chains = 1, iter = 1000,
#                sample_prior = "yes",
#                file="models/bsr_brm_test.rds",
#                file_refit = "on_change")
# bsr_brm
# plot(conditional_effects(bsr_brm), points = T)
