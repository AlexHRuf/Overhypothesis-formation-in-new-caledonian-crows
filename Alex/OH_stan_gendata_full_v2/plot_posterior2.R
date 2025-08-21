library(tidyverse)
library(tidybayes)

# --- Handy labels from your generator / data list ---
treatments <- tibble(
  id = seq_len(dataList$N_id),
  treatment = dataList$treatment
)
tid_tbl <- tibble(id = seq_len(dataList$N_id), T_id = dataList$T_id)

# =============== α (a_s) ===============
alpha_long <- fit %>%
  gather_draws(a_s[id, food, trial]) %>%
  left_join(treatments, by = "id") %>%
  left_join(tid_tbl,    by = "id") %>%
  filter(trial <= T_id) %>%               # drop padded trials
  select(-T_id) %>%
  rename(alpha = .value)

# (A) Treatment-mean + 90% CI ribbon
alpha_plot <- alpha_long %>%
  group_by(treatment, food, trial) %>%
  mean_qi(alpha, .width = 0.9) %>%
  ggplot(aes(trial, alpha, color = treatment, fill = treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.18, color = NA) +
  facet_wrap(~ food, labeller = labeller(food = c(`1`="High-value", `2`="Low-value"))) +
  labs(title = "Posterior α trajectories (a_s)", y = "α", x = "Trial") +
  theme_minimal()

# =============== β (a_b) ===============
beta_long_raw <- fit %>%
  gather_draws(a_b[id, food, trial]) %>%
  left_join(treatments, by = "id") %>%
  left_join(tid_tbl,    by = "id") %>%
  filter(trial <= T_id) %>%
  select(-T_id) %>%
  rename(beta = .value)

# IMPORTANT: filter out padded rows by enforcing the simplex sum==1 within (draw,id,trial)
beta_long <- beta_long_raw %>%
  group_by(.draw, id, trial) %>%
  mutate(sum_beta = sum(beta)) %>%
  ungroup() %>%
  filter(abs(sum_beta - 1) < 1e-8) %>%   # remove padded/invalid rows
  select(-sum_beta)

beta_plot <- beta_long %>%
  group_by(treatment, food, trial) %>%
  mean_qi(beta, .width = 0.9) %>%
  ggplot(aes(trial, beta, color = treatment, fill = treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.18, color = NA) +
  facet_wrap(~ food, labeller = labeller(food = c(`1`="High-value", `2`="Low-value"))) +
  labs(title = "Posterior β trajectories (a_b)", y = "β", x = "Trial") +
  theme_minimal()

# =============== P(Switch) ===============
# Map observation index s -> (id, trial) using the original data
index_df <- tibble(
  s     = seq_len(dataList$N_test),
  id    = dataList$id,
  trial = dataList$trial
) %>% left_join(treatments, by = "id")

switch_long <- fit %>%
  gather_draws(Prob_S[s, k]) %>%
  filter(k == 2L) %>%                 # k=2 is switch
  left_join(index_df, by = "s") %>%
  rename(P_switch = .value)

switch_plot <- switch_long %>%
  group_by(treatment, trial) %>%
  mean_qi(P_switch, .width = 0.9) %>%
  ggplot(aes(trial, P_switch, color = treatment, fill = treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.18, color = NA) +
  labs(title = "Posterior P(switch)", y = "P(Switch)", x = "Trial") +
  theme_minimal()

# -------- Show the plots you want --------
alpha_plot
beta_plot
switch_plot


library(tidyverse); library(tidybayes)

trial_eff <- fit %>% gather_draws(trial_effect)
trial_eff %>% mean_qi(.value, .width = 0.9)
ggplot(trial_eff, aes(.value)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Posterior of trial_effect", x = "trial_effect")
