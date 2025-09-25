library(dplyr)
library(ggplot2)
library(tidybayes)
library(tidytext)

source("extract_OH_draws.R")#no stickiness

res_no_kappa <-  extract_OH_draws(fit_no_kappa,dataList, width = 0.95)
res_global_kappa <-  extract_OH_draws(fit_global_kappa,dataList, width = 0.95)
res_id_kappa <-  extract_OH_draws(fit_id_kappa,dataList, width = 0.95)
res_fixed_kappa <-  extract_OH_draws(fit_fixed_kappa,dataList, width = 0.95)

res <- res_fixed_kappa
# -------------- plots --------------
#Global parameter plots
global_long <- res$draws$scalars %>%
  select(-.chain, -.iteration, -.draw) %>%
  tidyr::pivot_longer(everything(), names_to = "param", values_to = "value")
str(global_long)
ggplot(global_long, aes(value)) +
  geom_density() +
  facet_wrap(~ param, scales = "free") +
  theme_minimal()

# Lambda (raw)
ggplot(res$draws$re_id_raw, aes(x = lambda_raw_id, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.6) +
  labs(x = "lambda_raw_id", y = "density", title = "λ (raw) by treatment") +
  theme_minimal()

# Phi (res_id_kappa# Phi (raw)
ggplot(res$draws$re_id_raw, aes(x = phi_raw_id, fill = treatment, colour = treatment)) +
  geom_density(alpha = 0.6) +
  labs(x = "phi_raw_id", y = "density", title = "φ (raw) by treatment") +
  theme_minimal()

# Kappa (raw) – skip if you fixed κ (column will be NA)
if ("kappa_raw_id" %in% names(res$draws$re_id_raw) && any(is.finite(res$draws$re_id_raw$kappa_raw_id))) {
  ggplot(res$draws$re_id_raw, aes(x = kappa_raw_id, fill = treatment, colour = treatment)) +
    geom_density(alpha = 0.2) +
    labs(x = "kappa_raw_id", y = "density", title = "κ (raw) by treatment") +
    theme_minimal()
}

#Alpha (pre test) per ID
ggplot(res$summaries$sum_alpha, aes(x = reorder(id, rev(id)), y = mean, ymin = qlo, ymax = qhi, colour = treatment)) +
  geom_pointrange(linewidth = 0.5) +
  coord_flip() +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_x_reordered() +
  labs(x = "Individual",
       y = expression(alpha~"— posterior mean ± 80%"),
       title = expression(~alpha~" by ID")) +
  theme_minimal()


#Beta (pre test) per-ID and food type
ggplot(res$summaries$sum_beta,
       aes(x = reorder(id, rev(id)), y = mean, ymin = qlo, ymax = qhi, colour = treatment)) +
  geom_pointrange(position = position_dodge(width = 1)) +
  facet_wrap(~ food) +
  coord_flip() +
  labs(x = "Individual", y = "β (mean ± 80%)",
       title = expression(~beta~" by ID and food")) +
  theme_minimal()

#Attraction (pre test) per ID
attr_mean <- res$summaries$sum_ALL_F %>%
  dplyr::group_by(food) %>%
  dplyr::summarise(y = mean(mean), .groups = "drop")

ggplot(res$summaries$sum_ALL_F,
       aes(x = reorder(id, rev(id)), y = mean, ymin = qlo, ymax = qhi, colour = treatment)) +
  geom_pointrange(position = position_dodge(width = 1)) +
  facet_wrap(~ food, nrow = 1) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_hline(data = attr_mean, aes(yintercept = y),
             linetype = 2, linewidth = 0.6, colour = "grey40") +
  labs(x = "Individual", y = "Attraction (0–1)",
       title = "Pre-test attraction by ID and food") +
  theme_minimal()

#Prob_S by id and test trial
ggplot(res$summaries$sum_probS, aes(trial, mean, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.03, colour = NA) +
  labs(x = "Trial", y = "P(switch)") +
  ylim(0,1)+
  xlim(0,30)+
  guides(fill = "none") +
  theme_minimal()

#Calculate cumulative switch prob
sprob <- res$summaries$sum_probS

sprob <- sprob %>%
  arrange(treatment, id, trial) %>%             # ensure order
  group_by(treatment, id) %>%                   # reset per bird
  mutate(
    p_switch = mean,                            # rename for clarity
    p_stay   = 1 - p_switch,
    surv     = cumprod(p_stay),                 # P(still on B11 by trial t)
    cumprob_switch = 1 - surv,                  # P(ever switched by t)
    cumprob_stay   = surv
  ) %>%
  ungroup()


#Cumulative Prob SWITCH by id and test trial
ggplot(sprob, aes(trial, cumprob_switch, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  labs(x = "Trial", y = "P(switch)") +
  ylim(0,1)+
  xlim(0,30)+
  guides(fill = "none") +
  theme_minimal()

#Cumulative Prob STAY by id and test trial
ggplot(sprob, aes(trial, cumprob_stay, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  labs(x = "Trial", y = "P(switch)") +
  ylim(0,1)+
  xlim(0,30)+
  guides(fill = "none") +
  theme_minimal()

# alpha_t by id across trials
ggplot(res$summaries$sum_a_t, aes(trial, mean, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.03, colour = NA) +
  labs(x = "Trial", y = expression(alpha[t])) +
  theme_minimal()

# beta_t by id across trials
ggplot(res$summaries$sum_b_t, aes(trial, mean, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.05, colour = NA) +
  facet_wrap(~ food, nrow = 1) +
  labs(x = "Trial", y = expression(beta[t])) +
  ylim(0, 1) +
  theme_minimal()

# AF_t (food attractoion) by id across trials
ggplot(res$summaries$sum_AF_t, aes(trial, mean, color = treatment, fill = treatment, group = id)) +
  geom_line() +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.03, colour = NA) +
  facet_wrap(~ food, nrow = 1) +
  labs(x = "Trial", y = "Attraction (A_F)") +
  ylim(0, 1) +
  theme_minimal()

