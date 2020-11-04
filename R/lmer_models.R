library(papaja)
library(scales)
library(broom)
library(broom.mixed)
library(here)
library(tidyverse)
library(nlme)
library(lme4)
library(kableExtra)
library(lmerTest)
library(merTools)
library(ggsci)
library(huxtable)
library(ggeffects)
library(mgcv)
options(scipen=999)
options(huxtable.long_minus = TRUE)
theme_set(theme_apa(base_size = 14) + theme(legend.position = "bottom"))

path <- here("mod_data", "NPS_data.rds")
df <- read_rds(path)

# make a dummy columns of exploratory phenotype 
df$viieminutine_lqik <- as.numeric(as.character(df$viieminutine_lqik))
df$exp_phenotype <- as.ordered(df$HE_LE)

# make first and second order orthogonal polynomials of time bins
t <- poly(unique(df$viieminutine_lqik), 4)
 # create orthogonal polynomial time variables in data frame
df[, paste0("ot", 1:4)] <- t[df$viieminutine_lqik, 1:4]

df_short <- df %>% dplyr::filter(kestvuskategooria %in% 
                                   c("<12 ms", "12-99 ms", "100-299 ms")) %>% 
  dplyr::mutate(kestvuskategooria = fct_drop(kestvuskategooria))

df_long <- df %>% dplyr::filter(kestvuskategooria %in% 
                                  c("300-999 ms", "≥1000 ms")) %>% 
  dplyr::mutate(kestvuskategooria = fct_drop(kestvuskategooria))

df_summary <- df %>% group_by(HE_LE, rott, viieminutine_lqik,
                              duration_cat) %>% 
  summarise(peak_frequency = mean(average_peak_frequency), 
            peak_amplitude = mean(average_peak_amplitude), 
            duration = mean(duration_ms), vocalisations = n())
################################################################################
# models of duration
# short calls
# using lme with autocorrelated residuals
# best model with polynomial decomposition of time and order 1 autocorrelation 
# of residuals 
m_dur_short_1 <- lme(duration_ms ~ (ot1 + ot2 + ot3)*HE_LE,
                     random = list(rott =~ ot1 + ot2),
                     data = df_short, method = "REML",
                     correlation = corAR1(form = ~ 1|rott))
summary(m_dur_short_1)
plot(m_dur_short_1, rott~resid(.))
qqnorm(m_dur_short_1, ~resid(.))
qqnorm(m_dur_short_1, ~ranef(.))
path <- here("models", "m_dur_short_lme.rds")
saveRDS(m_dur_short_1, path)

g_m_dur_short_1 <- ggplot(df_short,
             aes(viieminutine_lqik, duration_ms, shape = HE_LE, group = HE_LE,
                 colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_dur_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Duration (ms)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model" )
g_m_dur_short_1

path <- here("graphs", "g_dur_short_lme.rds")
saveRDS(g_m_dur_short_1, path)
################################################################################
# fitting a smooth function of time and random intercept and slope
# this model fits better
fmla_dur_short_2 <- duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr', 
                                            k = 6, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_dur_short_2 <- bam(fmla_dur_short_2, data = df_short, nthreads = 2, 
                     family = gaussian(link = "log"), method = "REML")
summary(m_dur_short_2)
AIC(m_dur_short_2)
gam.check(m_dur_short_2)
plot(m_dur_short_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_dur_short_gam.rds")
saveRDS(m_dur_short_2, path)

g_m_dur_short_2 <- ggplot(df_short,
                         aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                             group = HE_LE,
                             colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_dur_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Duration (ms)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_dur_short_2

path <- here("graphs", "g_dur_short_gam.rds")
saveRDS(g_m_dur_short_2, path)
################################################################################
# long calls
m_dur_long_1 <- lme(duration_ms ~ (ot1 + ot2 + ot3)*HE_LE,
                    random = list(rott =~ ot1 + ot2 + ot3),
                    data = df_long, method = "REML",
                    correlation = corAR1(form = ~ 1|rott),
                    control = lmeControl(opt = "optim", optimMethod = "SANN"))
summary(m_dur_long_1)
path <- here("models", "m_dur_long_lme.rds")
saveRDS(m_dur_long_1, path)

g_m_dur_long_1 <- ggplot(df_long,
                        aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                            group = HE_LE,
                            colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_dur_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Duration (ms)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model" )
g_m_dur_long_1

path <- here("graphs", "g_dur_long_lme.rds")
saveRDS(g_m_dur_long_1, path)
################################################################################
# fitting a smooth function of time and random intercept and slope
# this model fits better
fmla_dur_long_2 <- duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr', 
                                                    k = 5, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_dur_long_2 <- bam(fmla_dur_long_2, data = df_long, nthreads = 2, 
                     family = Gamma(), method = "REML")
summary(m_dur_long_2)
AIC(m_dur_long_2)
gam.check(m_dur_long_2)
plot(m_dur_long_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_dur_long_gam.rds")
saveRDS(m_dur_long_2, path)

g_m_dur_long_2 <- ggplot(df_long,
                          aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_dur_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Duration (ms)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_dur_long_2

path <- here("graphs", "g_dur_long_gam.rds")
saveRDS(g_m_dur_long_2, path)
################################################################################
################################################################################
# models of frequency
# short calls
m_freq_short_1 <- lme(average_peak_frequency ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2 + ot3),
                      data = df_short, method = "REML",
                      correlation = corAR1(form = ~ 1|rott))
summary(m_freq_short_1)
path <- here("models", "m_freq_short_lme.rds")
saveRDS(m_freq_short_1, path)

g_m_freq_short_1 <- ggplot(df_short,
                          aes(viieminutine_lqik, average_peak_frequency, 
                              shape = HE_LE, group = HE_LE,
                              colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_freq_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak frequency (Hz)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model")
g_m_freq_short_1

path <- here("graphs", "g_freq_short_lme.rds")
saveRDS(g_m_freq_short_1, path)
################################################################################
fmla_freq_short_2 <- average_peak_frequency ~ exp_phenotype + 
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_freq_short_2 <- bam(fmla_freq_short_2, data = df_short, nthreads = 2, 
                     family = Gamma(), method = "REML")
summary(m_freq_short_2)
AIC(m_freq_short_2)
gam.check(m_freq_short_2)
plot(m_freq_short_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_freq_short_gam.rds")
saveRDS(m_freq_short_2, path)

g_m_freq_short_2 <- ggplot(df_short,
                          aes(viieminutine_lqik, average_peak_frequency, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_freq_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak frequency (Hz)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_freq_short_2

path <- here("graphs", "g_freq_short_gam.rds")
saveRDS(g_m_freq_short_2, path)
################################################################################
# long calls
m_freq_long_1 <- lme(average_peak_frequency ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2 + ot3),
                      data = df_long, method = "ML",
                      correlation = corAR1(form = ~ 1|rott))
summary(m_freq_long_1)
path <- here("models", "m_freq_long_lme.rds")
saveRDS(m_freq_long_1, path)

g_m_freq_long_1 <- ggplot(df_long,
                           aes(viieminutine_lqik, average_peak_frequency, 
                               shape = HE_LE, group = HE_LE,
                               colour = HE_LE)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(aes(y = fitted(m_freq_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak frequency (Hz)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model")
g_m_freq_long_1

path <- here("graphs", "g_freq_long_lme.rds")
saveRDS(g_m_freq_long_1, path)
################################################################################
fmla_freq_long_2 <- average_peak_frequency ~ exp_phenotype + 
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_freq_long_2 <- bam(fmla_freq_long_2, data = df_long, nthreads = 2, 
                      family = Gamma(), method = "REML")
summary(m_freq_long_2)
AIC(m_freq_long_2)
gam.check(m_freq_long_2)
plot(m_freq_long_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_freq_long_gam.rds")
saveRDS(m_freq_long_2, path)

g_m_freq_long_2 <- ggplot(df_long,
                           aes(viieminutine_lqik, average_peak_frequency, shape = HE_LE, 
                               group = HE_LE,
                               colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_freq_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak frequency (Hz)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_freq_long_2

path <- here("graphs", "g_freq_long_gam.rds")
saveRDS(g_m_freq_long_2, path)
################################################################################
################################################################################
# models of amplitude
# short calls
m_ampl_short_1 <- lme(average_peak_amplitude ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2),
                      data = df_short, method = "REML",
                      correlation = corAR1(form = ~ 1|rott))
summary(m_ampl_short_1)
path <- here("models", "m_ampl_short_lme.rds")
saveRDS(m_ampl_short_1, path)

g_m_ampl_short_1 <- ggplot(df_short,
                           aes(viieminutine_lqik, average_peak_amplitude, 
                               shape = HE_LE, group = HE_LE,
                               colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_ampl_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak amplitude (dB)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model")
g_m_ampl_short_1

path <- here("graphs", "g_ampl_short_lme.rds")
saveRDS(g_m_ampl_short_1, path)
################################################################################
fmla_ampl_short_2 <- average_peak_amplitude ~ exp_phenotype + 
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_ampl_short_2 <- bam(fmla_ampl_short_2, data = df_short, nthreads = 2, 
                      family = gaussian(), method = "REML")
summary(m_ampl_short_2)
AIC(m_ampl_short_2)
gam.check(m_ampl_short_2)
plot(m_ampl_short_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_ampl_short_gam.rds")
saveRDS(m_ampl_short_2, path)

g_m_ampl_short_2 <- ggplot(df_short,
                           aes(viieminutine_lqik, average_peak_amplitude, shape = HE_LE, 
                               group = HE_LE,
                               colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_ampl_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak amplitude (dB)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_ampl_short_2

path <- here("graphs", "g_ampl_short_gam.rds")
saveRDS(g_m_ampl_short_2, path)
################################################################################
# long calls
m_ampl_long_1 <- lme(average_peak_amplitude ~ (ot1 + ot2 + ot3)*HE_LE,
                     random = list(rott =~ ot1 + ot2 + ot3),
                     data = df_long, method = "ML",
                     correlation = corAR1(form = ~ 1|rott))
summary(m_ampl_long_1)
path <- here("models", "m_ampl_long_lme.rds")
saveRDS(m_ampl_long_1, path)

g_m_ampl_long_1 <- ggplot(df_long,
                          aes(viieminutine_lqik, average_peak_amplitude, 
                              shape = HE_LE, group = HE_LE,
                              colour = HE_LE)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(aes(y = fitted(m_ampl_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak amplitude (dB)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Third order polynomial model")
g_m_ampl_long_1

path <- here("graphs", "g_ampl_long_lme.rds")
saveRDS(g_m_ampl_long_1, path)
################################################################################
fmla_ampl_long_2 <- average_peak_amplitude ~ exp_phenotype + 
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) + 
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_ampl_long_2 <- bam(fmla_ampl_long_2, data = df_long, nthreads = 2, 
                     family = gaussian(), method = "REML")
summary(m_ampl_long_2)
AIC(m_ampl_long_2)
gam.check(m_ampl_long_2)
plot(m_ampl_long_2, shade = TRUE, pages = 1, scale = 0)
path <- here("models", "m_ampl_long_gam.rds")
saveRDS(m_ampl_long_2, path)

g_m_ampl_long_2 <- ggplot(df_long,
                          aes(viieminutine_lqik, average_peak_amplitude, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(aes(y = fitted(m_ampl_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line") + scale_color_jco() +
  labs(x = "Time bin (5 min)",
       y = "Average peak amplitude (dB)",
       shape = "Phenotype", linetype = "Phenotype", colour = "Phenotype",
       title = "Cubic regression spline model")
g_m_ampl_long_2

path <- here("graphs", "g_ampl_long_gam.rds")
saveRDS(g_m_ampl_long_2, path)
################################################################################