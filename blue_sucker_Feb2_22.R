# load packages
library(Matrix)
library(dbplyr)
library(tidyverse)
library(brms)
library(dplyr)
library(janitor)
library(readr)
library(tidybayes)
library(ggplot2)


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

###### Graphing, visualizing the RAW data #######

d <- blue_sucker_2021_data #simplifying what it's called.

# Lengths and weights for both sexes: 

LengthsWeights <- d %>% 
ggplot(aes(x=length_mm, y=weight_g, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="Lengths and Weights by Sex, 2021",
       x="Length (mm)",
       y= "Wet Weight (g)")

ggsave(LengthsWeights, file = "plots/LengthsWeighs.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Gonad weight and length for both sexes: 

GonadWeightsLengths <- d %>% 
  ggplot(aes(x=length_mm, y=gonad_weight_g, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="Lengths and Gonad Weights by Sex, 2021",
       x="Length (mm)",
       y= "Gonad Weight (g)")

ggsave(GonadLengthsWeights, file = "plots/GonadLengthsWeights.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Gonad weight and wet weight for both sexes:
GonadWeightsWet <- d %>% 
  ggplot(aes(x=weight_g, y=gonad_weight_g, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="Wet Weight and Gonad Weights by Sex, 2021",
       x="Wet weight (g)",
       y= "Gonad Weight (g)")

ggsave(GonadWeightsWet, file = "plots/GonadWeightsWet.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Weight vs total egg count

WeightEggCount <- d %>% 
  ggplot(aes(x=weight_s, y=combined_egg_total)) +
  geom_point() +
  geom_smooth()+
  labs(title="combined egg total by Weight (g), 2021",
       x="Wet Weight (g)",
       y= "combined egg count")

ggsave(WeightEggCount, file = "plots/WeightEggCount.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Weight and GSI for both sexes

WeightGSI <- d %>% 
  ggplot(aes(x=weight_g, y=gsi, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="GSI by Sex and Weight (g), 2021",
       x="Wet Weight (g)",
       y= "GSI")

ggsave(WeightGSI, file = "plots/WeightGSI.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Length and GSI for both sexes

LengthGSI <- d %>% 
  ggplot(aes(x=length_mm, y=gsi, color=lab_sex)) +
  geom_point() +
  geom_smooth()+
  labs(title="GSI by Sex and Length (mm), 2021",
       x="Length (mm)",
       y= "GSI")

ggsave(LengthGSI, file = "plots/LengthGSI.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Standardized length and GSI for both sexes

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

# Standard length and estimated egg total based on rounded approximations of egg counts

StandardLengthEggTotal <- d %>% 
  ggplot(aes(x=length_s, y=combined_egg_total)) +
  geom_point() +
  geom_smooth()+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  labs(title="Lengths and estimated total egg counts",
       x="Length (standardized)",
       y= "Estimated egg total")

ggsave(StandardLengthEggTotal, file = "plots/StandardLengthEggTotal.png", dpi=750,  width = 5, height = 3,
       units = "in")

d %>% 
  ggplot(aes(x=length_mm, y=combined_egg_total)) +
  geom_point() +
  geom_smooth()+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  labs(title="Lengths and estimated total egg counts",
       x="Length (standardized)",
       y= "Estimated egg total")

# Egg total vs gonad weight of both ovaries combined

GonadWeightEgg <- d %>% 
  ggplot(aes(x=gonad_weight_g, y=combined_egg_total)) +
  geom_point() +
  geom_smooth()+
  labs(title="Combined ovary weight and estimated total egg counts",
       x="Total Gonad Weight (g)",
       y= "Estimated Egg Total")

ggsave(GonadWeightEgg, file = "plots/GonadWeightEgg.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Weight vs egg count per gram in ovary 1

O1WeightEgg <- d %>% 
  ggplot(aes(x=ovary_1_g, y=o1_average_count)) +
  geom_point() +
  geom_smooth()+
  labs(title="Ovary one",
       x="Ovary weight (g)",
       y= "Avg egg count per gram")+
  theme_linedraw()

ggsave(O1WeightEgg, file = "plots/O1WeightEgg.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Weight vs egg count per gram in ovary 2

O2WeightEgg <- d %>% 
  ggplot(aes(x=ovary_2_g, y=o2_average_count)) +
  geom_point() +
  geom_smooth()+
  labs(title="Ovary two",
       x="Ovary weight (g)",
       y= "Avg egg count per gram")+
  theme_linedraw()

ggsave(O2WeightEgg, file = "plots/O2WeightEgg.png", dpi=750,  width = 5, height = 3,
       units = "in")

# Checking to see if the ovaries have a bias of some sort in terms of measuring eggs
# need to fix this probably

O1 = data.frame(x = d$ovary_1_g,y=d$o1_average_count)
O2 = data.frame(x = d$ovary_2_g,y=d$o1_average_count)

Ov1Ov2comp <- ggplot(O1,aes(x,y)) +
  geom_point(alpha=0.6, color='darkolivegreen4') +
  geom_point(alpha=0.6, data=O2,color='dodgerblue3')+
  geom_smooth(alpha=0,color='darkolivegreen4')+
  geom_smooth(data=O2,alpha=0, color='dodgerblue3')+
  xlab("Ovary weight") +
  scale_y_continuous("Average egg count per gram", limits = c(150,350))+
  labs(title="Comparison of ovaries one (green) and two (blue)")+
  theme_linedraw()

# HOW to make the dots line up based on individuals (i.e. do the individuals have 
# differences in ovary 1 vs ovary 2)?

ggsave(Ov1Ov2comp, file = "plots/Ov1Ov2comp.png", dpi=750,  width = 5, height = 3,
       units = "in")

ovary_data <- d %>% 
  select(fish_id, ovary_1_g,ovary_2_g,o2_average_count,o1_average_count) %>% 
  mutate(ovary_diff = ovary_2_g - ovary_1_g,
         egg_diff = o2_average_count - o1_average_count)

mean(ovary_data$egg_diff, na.rm = TRUE)
# 0.5263158
mean(ovary_data$ovary_diff, na.rm = TRUE)
# -1.905789

# Showing the difference between ovaries

OvaryDiff <- ovary_data %>% 
  ggplot(aes(x=ovary_diff, y=egg_diff))+
  geom_hline(linetype="twodash", yintercept = 2.263158, size=0.7, color="cyan4")+
  geom_vline(linetype="twodash", xintercept = -1.905789, size=0.7, color="red3")+
  geom_point(shape=19)+
  labs(x="Difference in weight between ovaries one and two (g)**",
       y="Difference in egg count*",
       caption = "*Mean difference in egg count denoted with blue dashed line
**Mean difference in ovary weight denoted with red dashed line")+
  theme_linedraw()+
  theme(plot.caption = element_text(hjust = 0))

ggsave(OvaryDiff, file = "plots/OvaryDiff.png", dpi=750,  width = 5, height = 3,
       units = "in")

##### some prelimiary tests, can pretty much ignore these #####

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


###### What is fecundity? ######

# Fecundity is the number of eggs a female fish will lay in a spawning season.

###### 6 DIFFERENT MODELS ######
# 1) TOTAL LENGTH as predictor of TOTAL EGG COUNT
# 2) WET WEIGHT as predictor of TOTAL EGG COUNT
# 3) TOTAL LENGTH as predictor of GSI
# 4) WET WEIGHT as predictor of GSI
# 5) TOTAL LENGTH as predictor of GONAD WEIGHT
# 6) WET WEIGHT as predictor of GONAD WEIGHT

######## TOTAL LENGTH as predictor of TOTAL EGG COUNT ###########

# getting priors
get_prior(egg_total_simplified ~ length_s + I(length_s^2), 
           data = d,
           family = negbinomial(link="log"))

#simulating priors
priors = tibble(length_beta = rnorm(100, 0.3, 0.15),
Ilength_beta2 = rnorm(100,-0.15 ,0.05),
                Intercept = rnorm(100, 11.5, 0.25),
                iter = 1:100)

prior_sims = priors %>%
  expand_grid(d %>% distinct(length_s)) %>%
  mutate(count_sims = Intercept + length_beta*length_s + Ilength_beta2*(length_s^2))

ggplot() +
  geom_line(data=prior_sims, aes(x = length_s, y = count_sims, group = iter))+
  geom_point(data=d,aes(x=length_s, y=log(combined_egg_total),color="red"))

# making the model
# add + (1|year) when we incorporate the new data and a prior for sigma: prior(exponential(1), class="sigma")

length_bsr_negbinom <- brm(combined_egg_total ~ length_s + I(length_s^2), 
                           data = d,
                           family = negbinomial(link="log"),
                           prior = c(prior(normal(11.5, 0.25), class = "Intercept"),
                                     prior(normal(0.2, 0.15), class = "b", coef="length_s"),
                                     prior(normal(-0.15,0.05), class = "b", coef="Ilength_sE2"),
                                     prior(exponential(0.1), class = "shape")),
                           cores = 1, chains = 1, iter = 1000,
                           sample_prior = "yes",
                           file="models/length_bsr_negbinom.rds",
                           file_refit = "on_change")

# conditional effects, taking all individuals into account
plot(conditional_effects(length_bsr_negbinom, re_formula = NULL), points = T)
# conditional effects, showing the mean difference
plot(conditional_effects(length_bsr_negbinom), points = T)

summary(length_bsr_negbinom)

pp_check(length_bsr_negbinom)
pp_check(length_bsr_negbinom, type = "hist")

saveRDS(length_bsr_negbinom, "models/length_bsr_negbinom.rds")

# conditional effects, manual plotting
as_draws_df(length_bsr_negbinom)

cond_effect_length <- conditional_effects(length_bsr_negbinom)
cond_effect_length$length_s

cond_effect_length$lenth_s %>% 
  ggplot(aes(x=length_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = length_bsr_negbinom$data, aes(x=length_s, y=combined_egg_total))+
  theme_default()
# it keeps saying "estimate__ not found"

cond_data_length <- length_bsr_negbinom$data %>% distinct(length_s, combined_egg_total)

posts_length <- add_epred_draws(length_bsr_negbinom, newdata= length_bsr_negbinom$data %>% 
                                  distinct(length_s) , re_formula = NA)


posts_length_all <- add_predicted_draws(length_bsr_negbinom, newdata= length_bsr_negbinom$data %>% 
                                          distinct(length_s) , re_formula = NA)

d_length <- d %>% distinct(length_mm, length_s)

PosteriorLength <- posts_length_all %>%
  group_by(length_s) %>% 
  left_join(d_length) %>% 
  median_qi(.prediction) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .prediction)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = combined_egg_total)) +
  labs(title= "Blue Sucker Fecundity Prediction",
       subtitle="Large grey bar incorporates the variation in individuals",
       x="Length (mm)",
       y="Predicted total egg count")

ggsave(PosteriorLength, file = "plots/PosteriorLength.png", dpi = 750, width = 7, height = 5,
       units = "in")
# This model incorporates all individuals, and not JUST the mean.


PosteriorLengthMean <- posts_length %>%
  group_by(length_s) %>% 
  left_join(d_length) %>% 
  median_qi(.epred) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = combined_egg_total)) +
  labs(title= "Blue Sucker Mean Fecundity Prediction",
       subtitle="Grey bar incorporates only the variation in the mean egg count",
       x="Length (mm)",
       y="Predicted total egg count")

ggsave(PosteriorLengthMean, file = "plots/PosteriorLengthMean.png", dpi = 750, width = 7, height = 5,
       units = "in")



######### WET WEIGHT as predictor of TOTAL EGG COUNT ########## 

#prior simulation
priors = tibble(weight_beta = rnorm(100, 0.25, 0.08),
                # Iweight_beta2 = rnorm(100,-0.15 ,0.05),
                Intercept = rnorm(100, 11.25, 0.25),
                iter = 1:100)

prior_sims = priors %>%
  expand_grid(d %>% distinct(weight_s)) %>%
  mutate(count_sims = Intercept + weight_beta*weight_s)

ggplot() +
  geom_line(data=prior_sims, aes(x = weight_s, y = count_sims, group = iter))+
  geom_point(data=d,aes(x=weight_s, y=log(combined_egg_total),color="red"))

weight_bsr_negbinom <- brm(combined_egg_total ~ weight_s, 
                           data = d,
                           family = negbinomial(link="log"),
                           prior = c(prior(normal(11.25, 0.25), class = "Intercept"),
                                     prior(normal(0.25, 0.08), class = "b", coef="weight_s"),
                                     prior(exponential(0.1), class="shape")),
                           cores = 1, chains = 4, iter = 1000,
                           # sample_prior = "yes",
                           file="models/weight_bsr_negbinom.rds",
                           file_refit = "on_change")
summary(weight_bsr_negbinom)

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

posts_weight <- add_epred_draws(weight_bsr_negbinom, newdata = weight_bsr_negbinom$data %>% 
                                  distinct(weight_s) , re_formula = NA)


posts_weight_all <- add_predicted_draws(weight_bsr_negbinom, newdata= weight_bsr_negbinom$data %>% 
                  distinct(weight_s) , re_formula = NA)

d_weight <- d %>% distinct(weight_g, weight_s)

PosteriorWeight <- posts_weight_all %>%
  group_by(weight_s) %>% 
  left_join(d_weight) %>% 
  median_qi(.prediction) %>% 
  mutate(weight_g = (weight_s*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x = weight_g, y = .prediction)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = combined_egg_total)) +
  labs(title= "Blue Sucker Fecundity Prediction",
       subtitle="Large grey bar incorporates the variation in individuals",
       x="Weight (g)",
       y="Predicted total egg count")

ggsave(PosteriorWeight, file = "plots/PosteriorWeight.png", dpi = 750, width = 7, height = 5,
       units = "in")

# This model incorporates all individuals, and not JUST the mean. We would not be surprised
# to see any range of egg counts for an individual weighing (for example) 3000 g, to be between
# ~75,000 and ~130,000 eggs. For the MEAN variation:

PosteriorWeightMean <- posts_weight %>%
  group_by(weight_s) %>% 
  left_join(d_weight) %>% 
  median_qi(.epred) %>% 
  mutate(weight_g = (weight_s*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x = weight_g, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = combined_egg_total)) +
  labs(title= "Blue Sucker Mean Fecundity Prediction",
       subtitle="Grey bar incorporates only the variation in the MEAN egg count",
       x="Weight (g)",
       y="Predicted total egg count")

ggsave(PosteriorWeightMean, file = "plots/PosteriorWeightMean.png", dpi = 750, width = 7, height = 5,
       units = "in")




######## TOTAL LENGTH as predictor of GSI ###########
#prior simulation

priors = tibble(beta = rnorm(100, 0.5, 0.2),
                beta2 = rnorm(100,-0.5 ,0.2),
                Intercept = rnorm(100, 13, 2),
                beta3 = rnorm(100, -8, 1),
                beta4 = rnorm(100,-0.4,0.5),
                beta5 = rnorm(100,0.4, 0.2),
                iter = 1:100)

prior_sims = priors %>%
  expand_grid(d %>% distinct(length_s, lab_sex)) %>%
  expand_grid(sex=c(0,1)) %>%
  mutate(gsi_sims = Intercept + beta*length_s + beta2*(length_s^2) + beta3*sex+
           beta4*length_s*sex + beta5*(length_s^2)*sex)

ggplot() +
  
  geom_line(data=prior_sims, aes(x = length_s, y = gsi_sims, group = interaction(sex,iter),
                                 
                                 color=as.factor(sex)))  #+

#geom_point(data=d,aes(x=length_s, y=gsi,shape=lab_sex))



get_prior(gsi ~ (length_s + I(length_s^2))*lab_sex,
          data = d,
          family = gaussian())



gsi_length <- brm(gsi ~ (length_s + I(length_s^2)) * lab_sex,
                  data = d,
                  family = gaussian(),
                  prior = c(prior(normal(13,2),class="Intercept"),
                            prior(normal(-0.5, 0.25), coef ="Ilength_sE2"),
                            prior(normal(-8, 1), coef = "lab_sexM"),
                            prior(normal(0.5, 0.25), coef = "length_s"),
                            prior(normal(0.4, 0.2), coef = "Ilength_sE2:lab_sexM"),
                            prior(normal(-0.4, 0.2), coef = "length_s:lab_sexM"),
                            prior(exponential(0.1), class="sigma")),
                  cores = 4, chains = 1, iter = 1000)
                  #sample_prior = "only")

summary(gsi_length)

plot(conditional_effects(gsi_length, re_formula=NULL), points = T)

gsi_length
pp_check(gsi_length)

cond_effect_gsi_l <- conditional_effects(gsi_length)
cond_effect_gsi_l$length_s

cond_effect_gsi_l$length_s %>% 
  ggplot(aes(x=length_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = gsi_length$data, aes(x=length_s, y=gsi))+
  theme_default()

cond_data_gsi_l <- gsi_length$data %>% distinct(length_s, gsi, fish_id)

posts_gsi_l <- add_epred_draws(gsi_length, newdata= gsi_length$data %>% 
                                 distinct(length_s, fish_id, lab_sex) , re_formula = NA)


posts_gsi_all <- add_predicted_draws(gsi_length, newdata=gsi_length$data %>% 
                                       distinct(length_s,fish_id,lab_sex) , re_formula = NA)

d_lengthgsi <- d %>% distinct(length_mm, length_s)

PosteriorGSIlength <- posts_gsi_all %>%
  group_by(length_s, lab_sex) %>% 
  left_join(d_lengthgsi) %>% 
  median_qi(.prediction) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x =length_mm, y = .prediction, fill = lab_sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker GSI Prediction",
       subtitle="Blue and pink bars incorporate the variation in individuals",
       x="Length (mm)",
       y="Predicted GSI")

ggsave(PosteriorGSIlength, file = "plots/PosteriorGSIlength.png", dpi = 750, width = 7, height = 5, units = "in")
# This model incorporates all individuals, and not JUST the mean. I LOVE the way this one looks.

PosteriorGSIlengthMean <- posts_gsi_l %>%
  group_by(length_s) %>% 
  left_join(d_lengthgsi) %>% 
  median_qi(.epred) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker Mean GSI Prediction",
       subtitle="Grey bar incorporates only the variation in the mean GSI",
       x="length (mm)",
       y="Predicted GSI")

#prior simulation
priors = tibble(beta = rnorm(100, 0.5, 0.2),
                beta2 = rnorm(100,-0.5 ,0.2),
                Intercept = rnorm(100, 13, 2),
                beta3 = rnorm(100, -8, 1),
                beta4 = rnorm(100,-0.4 ,0.5),
                beta5 = rnorm(100,-0.4 ,0.2),
                iter = 1:100)

prior_sims = priors %>%
  expand_grid(d %>% distinct(length_s)) %>%
  expand_grid(sex=c(0,1)) %>% 
  # mutate(prior_sexb = case_when(sex == 0 ~ 0,
  #                              TRUE ~ prior_sex)) %>% 
  mutate(gsi_sims = Intercept + beta*length_s +beta2*(length_s^2) + beta3*sex)

ggplot() +
  geom_line(data=prior_sims, aes(x = length_s, y = gsi_sims, group = interaction(sex,iter), 
                                 color=as.factor(sex)))+
  geom_point(data=d,aes(x=length_s, y=gsi,shape=lab_sex))


get_prior(gsi ~ length_s + (I(length_s^2))*lab_sex,
          data = d,
          family = gaussian())

gsi_length <- brm(gsi ~ (length_s + I(length_s^2) + (1 + length_s + I(length_s^2)|lab_sex)), 
                           data = d,
                           family = gaussian(),
                  piror = c(prior())
                           cores = 1, chains = 1, iter = 1000,
                           sample_prior = "yes")

plot(conditional_effects(gsi_length, re_formula=NULL), points = T)

gsi_length
pp_check(gsi_length)

cond_effect_gsi_l <- conditional_effects(gsi_length)
cond_effect_gsi_l$length_s

cond_effect_gsi_l$length_s %>% 
  ggplot(aes(x=length_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = gsi_length$data, aes(x=length_s, y=gsi))+
  theme_default()

cond_data_gsi_l <- gsi_length$data %>% distinct(length_s, gsi, fish_id)

posts_gsi_l <- add_epred_draws(gsi_length, newdata= gsi_length$data %>% 
                                  distinct(length_s, fish_id, lab_sex) , re_formula = NA)


posts_gsi_all <- add_predicted_draws(gsi_length, newdata=gsi_length$data %>% 
                                          distinct(length_s,fish_id,lab_sex) , re_formula = NA)

d_lengthgsi <- d %>% distinct(length_mm, length_s)

PosteriorGSIlength <- posts_gsi_all %>%
  group_by(length_s, lab_sex) %>% 
  left_join(d_lengthgsi) %>% 
  median_qi(.prediction) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x =length_mm, y = .prediction, fill = lab_sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker GSI Prediction",
       subtitle="Blue and pink bars incorporate the variation in individuals",
       x="Length (mm)",
       y="Predicted GSI")

ggsave(PosteriorGSIlength, file = "plots/PosteriorGSIlength.png", dpi = 750, width = 7, height = 5, units = "in")
# This model incorporates all individuals, and not JUST the mean. I LOVE the way this one looks.

PosteriorGSIlengthMean <- posts_gsi_l %>%
  group_by(length_s) %>% 
  left_join(d_lengthgsi) %>% 
  median_qi(.epred) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker Mean GSI Prediction",
       subtitle="Grey bar incorporates only the variation in the mean GSI",
       x="length (mm)",
       y="Predicted GSI")
# now it's saying that "lab_sex" cant be found

# ggsave(PosteriorGSIlengthMean, file = "plots/PosteriorGSIlengthMean.png", dpi = 750, width = 7, height = 5, units = "in")


######## WET WEIGHT as predictor of GSI ###########

get_prior(gsi ~ weight_s*lab_sex + I(weight_s^2) + (1|fish_id),
          data = d,
          family = gaussian())

gsi_weight <- brm(gsi ~ weight_s*lab_sex + I(weight_s^2) + (1|fish_id), 
                  data = d,
                  family = gaussian(),
                  cores = 4, chains = 4, iter = 7500,
                  sample_prior = "yes")

plot(conditional_effects(gsi_weight, re_formula=NULL), points = T)

gsi_weight
pp_check(gsi_weight)

cond_effect_gsi_w <- conditional_effects(gsi_weight)
cond_effect_gsi_w$weight_s

cond_effect_gsi_w$weight_s %>% 
  ggplot(aes(x=weight_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = gsi_weight$data, aes(x=weight_s, y=gsi))+
  theme_default()

cond_data_gsi_w <- gsi_weight$data %>% distinct(weight_s, gsi, fish_id)

posts_gsi_w <- add_epred_draws(gsi_weight, newdata= gsi_weight$data %>% 
                                 distinct(weight_s, fish_id, lab_sex) , re_formula = NA)


posts_gsi_allw <- add_predicted_draws(gsi_weight, newdata=gsi_weight$data %>% 
                                       distinct(weight_s,fish_id,lab_sex) , re_formula = NA)

d_weightgsi <- d %>% distinct(weight_g, weight_s)

PosteriorGSIweight<- posts_gsi_allw %>%
  group_by(weight_s, lab_sex) %>% 
  left_join(d_weightgsi) %>% 
  median_qi(.prediction) %>% 
  mutate(weight_g = (weight_s*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x =weight_g, y = .prediction, fill = lab_sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker GSI Prediction",
       subtitle="Blue and pink bars incorporate the variation in individuals",
       x="Weight (g)",
       y="Predicted GSI")

ggsave(PosteriorGSIweight, file = "plots/PosteriorGSIweight.png", dpi = 750, width = 7, height = 5, units = "in")
# This model incorporates all individuals, and not JUST the mean. I LOVE the way this one looks.

PosteriorGSIweightMean <- posts_gsi_w %>%
  group_by(weight_s) %>% 
  left_join(d_weightgsi) %>% 
  median_qi(.epred) %>% 
  mutate(weight_g = (weight_s*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x = weight_g, y = .epred), fill=lab_sex) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gsi)) +
  labs(title= "Blue Sucker Mean GSI Prediction",
       subtitle="Grey bar incorporates only the variation in the mean GSI",
       x="length (mm)",
       y="Predicted GSI")
# now it's saying that "lab_sex" cant be found

# ggsave(PosteriorGSIweightMean, file = "plots/PosteriorGSIweightMean.png", dpi = 750, width = 7, height = 5, units = "in")



######## TOTAL LENGTH as predictor of GONAD WEIGHT ###########

get_prior(gonad_weight_g ~ length_s*lab_sex + length_s + I(length_s^2) + (1|fish_id), 
          data = d,
          family = gaussian())


length_gonad_weight <- brm(gonad_weight_g ~ length_s*lab_sex + length_s + I(length_s^2) + (1|fish_id), 
                           data = d,
                           family = gaussian(),
                           cores = 4, chains = 4, iter = 7500,
                           sample_prior = "yes")
                           # file="models/length_gonad_weight.rds",
                           # file_refit = "on_change")

plot(conditional_effects(length_gonad_weight, re_formula = NULL), points = T)

length_gonad_weight
plot(conditional_effects(length_gonad_weight), points = T)

pp_check(length_gonad_weight)
pp_check(length_gonad_weight, type = "hist")

saveRDS(length_gonad_weight, "models/length_gonad_weight.rds")

as_draws_df(length_gonad_weight)

cond_effect_gonadwlength <- conditional_effects(length_gonad_weight)
cond_effect_gonadwlength$length_s

cond_effect_gonadwlength$length_s %>% 
  ggplot(aes(x=length_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = length_gonad_weight$data, aes(x=length_s, y=gonad_weight_g))+
  theme_default()

cond_data_gonadwlength <- length_gonad_weight$data %>% distinct(length_s, gonad_weight_g)

posts_gonadwlength <- add_epred_draws(length_gonad_weight, newdata= length_gonad_weight$data %>% 
                                  distinct(length_s, lab_sex) , re_formula = NA)


posts_length_gonadw <- add_predicted_draws(length_gonad_weight, newdata= length_gonad_weight$data %>% 
                                          distinct(length_s, lab_sex) , re_formula = NA)

d_lengthgonadw <- d %>% distinct(length_mm, length_s,gonad_weight_g)

PosteriorLengthGonadW <- posts_length_gonadw %>%
  group_by(length_s) %>% 
  left_join(d_lengthgonadw) %>% 
  median_qi(.prediction) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .prediction)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gonad_weight_g)) +
  labs(title= "Blue Sucker Gonad Weights and Lengths",
       subtitle="Large grey bar incorporates the variation in individuals",
       x="Length (mm)",
       y="Predicted gonad weight")

ggsave(PosteriorLengthGonadW, file = "plots/PosteriorLengthGonadW.png", dpi = 750, width = 7, height = 5,
       units = "in")
# This model incorporates all individuals, and not JUST the mean.
# need to get it to recognize the different sexes.

PosteriorLengthGonadWMean <- posts_gonadwlength %>%
  group_by(length_s) %>% 
  left_join(d_length) %>% 
  median_qi(.epred) %>% 
  mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>% 
  ggplot(aes(x = length_mm, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gonad_weight_g)) +
  labs(title= "Blue Sucker Gonad Weights and Lengths",
       subtitle="Grey bar incorporates only the variation in the mean gonad weight",
       x="Length (mm)",
       y="Predicted total gonad weight")

ggsave(PosteriorLengthGonadWMean, file = "plots/PosteriorLengthGonadWMean.png", dpi = 750, width = 7, height = 5,
       units = "in")


######## WET WEIGHT as predictor of GONAD WEIGHT ###########

get_prior(gonad_weight_g ~ weight_s*lab_sex + weight_s + (1|fish_id), 
          data = d,
          family = gaussian())


weight_gonad_weight <- brm(gonad_weight_g ~  weight_s*lab_sex + weight_s + (1|fish_id), 
                           data = d,
                           family = gaussian(),
                           cores = 4, chains = 4, iter = 7500,
                           sample_prior = "yes")
# file="models/weight_gonad_weight.rds",
# file_refit = "on_change")

plot(conditional_effects(weight_gonad_weight, re_formula = NULL), points = T)

weight_gonad_weight
plot(conditional_effects(weight_gonad_weight), points = T)

pp_check(weight_gonad_weight)
pp_check(weight_gonad_weight, type = "hist")

saveRDS(weight_gonad_weight, "models/weight_gonad_weight.rds")

as_draws_df(weight_gonad_weight)

cond_effect_weightgonadw <- conditional_effects(weight_gonad_weight)
cond_effect_weightgonadw$weight_s

cond_effect_weightgonadw$weight_s %>% 
  ggplot(aes(x=weight_s)) +
  geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
  geom_point(data = weight_gonad_weight$data, aes(x=weight_s, y=gonad_weight_g))+
  theme_default()

cond_data_weightgonadw <- weight_gonad_weight$data %>% distinct(weight_s, gonad_weight_g, weight_g)

posts_weightgonadw <- add_epred_draws(weight_gonad_weight, newdata= weight_gonad_weight$data %>% 
                                        distinct(weight_s, lab_sex) , re_formula = NA)


posts_weight_gonadw <- add_predicted_draws(weight_gonad_weight, newdata= weight_gonad_weight$data %>% 
                                             distinct(weight_s, lab_sex) , re_formula = NA)

d_weightgonadw <- d %>% distinct(weight_g, weight_s,gonad_weight_g)

PosteriorWeightGonadW <- posts_weight_gonadw %>%
  group_by(weight_s) %>% 
  left_join(d_weightgonadw) %>% 
  median_qi(.prediction) %>% 
  mutate(length_mm = (weight_s*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x = weight_g, y = .prediction)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gonad_weight_g)) +
  labs(title= "Blue Sucker Gonad Weights and Lengths",
       subtitle="Large grey bar incorporates the variation in individuals",
       x="Length (mm)",
       y="Predicted gonad weight")
# now it says weight_g can't be found


ggsave(PosteriorLengthGonadW, file = "plots/PosteriorWeightGonadW.png", dpi = 750, width = 7, height = 5,
       units = "in")
# This model incorporates all individuals, and not JUST the mean.
# need to get it to recognize the different sexes.

PosteriorWeightGonadWMean <- posts_gonadwlength %>%
  group_by(weight_s) %>% 
  left_join(d_weightgonadw) %>% 
  median_qi(.epred) %>% 
  mutate(weight_g = (weight_g*sd(d$weight_g)) + mean(d$weight_g)) %>% 
  ggplot(aes(x = weight_g, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.2) +
  geom_point(data = d, 
             aes(y = gonad_weight_g)) +
  labs(title= "Blue Sucker Gonad Weights and Wet Weight",
       subtitle="Grey bar incorporates only the variation in the mean gonad weight",
       x="Wet weight (g)",
       y="Predicted total gonad weight")

ggsave(PosteriorLengthGonadWMean, file = "plots/PosteriorLengthGonadWMean.png", dpi = 750, width = 7, height = 5,
       units = "in")




















###### Random other stuff ######

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
# Working in class
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
