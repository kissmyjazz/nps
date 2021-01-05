library(papaja)
library(scales)
library(broom)
library(broom.mixed)
library(here)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(merTools)
library(ggsci)
library(mgcv)
library(mgcViz)
library(extrafont)
library(cowplot)
library(ggplotify)
library(Cairo)
options(scipen=999)
options(huxtable.long_minus = TRUE)
theme_set(theme_apa(base_size = 14, base_family = "Times New Roman") + 
  theme(axis.title.y = element_text(size = 18), 
        axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", face = "plain"),
        legend.position = "none", 
        plot.margin = ggplot2::margin(b = 0.5, l = 0.5, r = 0.4, t = 0.7, unit = "cm"),
        axis.line = element_line(size = 0.8),
        legend.key.size = unit(0.1, "cm"),
        legend.background=element_blank()))

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
# ggplot constants
pd <- position_dodge(width = 0.15)
# x_scl <- scale_x_continuous(breaks = 1:6, 
#                    labels = c("0-5 min", "5-10 min", "10-15 min",
#                               "15-20 min", "20-25 min", "25-30 min"))
x_scl <- scale_x_continuous(breaks = 1:6, 
                            labels = c("0-5 min", "", "10-15 min", "",
                                       "20-25 min", ""))
c_scl <- scale_color_manual(values = c("LE" = "#008080", "HE" = "#ff8000"))
mean_stats <- stat_summary(fun = mean, geom = "point", size = 2.5, 
             position = pd)
cl_stats <- stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2,
               position = pd)
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
# summary(m_dur_short_1)
# plot(m_dur_short_1, rott~resid(.))
# qqnorm(m_dur_short_1, ~resid(.))
# qqnorm(m_dur_short_1, ~ranef(.))
# path <- here("models", "m_dur_short_lme.rds")
# saveRDS(m_dur_short_1, path)

g_m_dur_short_1 <- ggplot(df_short,
             aes(viieminutine_lqik, duration_ms, shape = HE_LE, group = HE_LE,
                 colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_dur_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl + 
  labs(x = NULL,
       y = "Duration (ms)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(legend.position = c(0.80, 0.90), axis.text.x = element_blank(), 
        axis.title.y = element_blank())

g_m_dur_short_1 

# path <- here("graphs", "g_dur_short_lme.rds")
# saveRDS(g_m_dur_short_1, path)
################################################################################
# fitting a smooth function of time and random intercept and slope
# this model fits better
fmla_dur_short_2 <- duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr',
                                            k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_dur_short_2 <- bam(fmla_dur_short_2, data = df_short, nthreads = 2,
                     family = gaussian(link = "log"), method = "REML")
# summary(m_dur_short_2)
# AIC(m_dur_short_2)
# gam.check(m_dur_short_2)
# plot(m_dur_short_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_dur_short_gam.rds")
# saveRDS(m_dur_short_2, path)

g_m_dur_short_2 <- ggplot(df_short,
                         aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                             group = HE_LE,
                             colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_dur_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Duration (ms)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.title.y = element_blank())

g_m_dur_short_2

# path <- here("graphs", "g_dur_short_gam.rds")
# saveRDS(g_m_dur_short_2, path)
################################################################################
# long calls
m_dur_long_1 <- lme(duration_ms ~ (ot1 + ot2 + ot3)*HE_LE,
                    random = list(rott =~ ot1 + ot2 + ot3),
                    data = df_long, method = "REML",
                    correlation = corAR1(form = ~ 1|rott),
                    control = lmeControl(opt = "optim", optimMethod = "SANN"))
# summary(m_dur_long_1)
# path <- here("models", "m_dur_long_lme.rds")
# saveRDS(m_dur_long_1, path)

g_m_dur_long_1 <- ggplot(df_long,
                        aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                            group = HE_LE,
                            colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_dur_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl + 
  x_scl +
  labs(x = NULL,
       y = "Duration (ms)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.text.x = element_blank())
g_m_dur_long_1

# path <- here("graphs", "g_dur_long_lme.rds")
# saveRDS(g_m_dur_long_1, path)
################################################################################
# fitting a smooth function of time and random intercept and slope
# this model fits better
fmla_dur_long_2 <- duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr',
                                                    k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_dur_long_2 <- bam(fmla_dur_long_2, data = df_long, nthreads = 2,
                     family = Gamma(), method = "REML")
# summary(m_dur_long_2)
# AIC(m_dur_long_2)
# gam.check(m_dur_long_2)
# plot(m_dur_long_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_dur_long_gam.rds")
# saveRDS(m_dur_long_2, path)

g_m_dur_long_2 <- ggplot(df_long,
                          aes(viieminutine_lqik, duration_ms, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_dur_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + c_scl +
  x_scl +
  labs(x = NULL,
       y = "Duration (ms)",
       shape = NULL, linetype = NULL, colour = NULL)
g_m_dur_long_2

# path <- here("graphs", "g_dur_long_gam.rds")
# saveRDS(g_m_dur_long_2, path)
################################################################################
################################################################################
# models of frequency
# short calls
m_freq_short_1 <- lme(average_peak_frequency ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2 + ot3),
                      data = df_short, method = "REML",
                      correlation = corAR1(form = ~ 1|rott))
# summary(m_freq_short_1)
# path <- here("models", "m_freq_short_lme.rds")
# saveRDS(m_freq_short_1, path)

g_m_freq_short_1 <- ggplot(df_short,
                          aes(viieminutine_lqik, average_peak_frequency, 
                              shape = HE_LE, group = HE_LE,
                              colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_freq_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\nfrequency (Hz)",
       shape = NULL, linetype = NULL, colour = NULL) +
  theme(legend.position = c(0.85, 0.22), axis.text.x = element_blank(), 
        axis.title.y = element_blank())
g_m_freq_short_1

# path <- here("graphs", "g_freq_short_lme.rds")
# saveRDS(g_m_freq_short_1, path)
################################################################################
fmla_freq_short_2 <- average_peak_frequency ~ exp_phenotype +
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_freq_short_2 <- bam(fmla_freq_short_2, data = df_short, nthreads = 2,
                     family = Gamma(), method = "REML")
# summary(m_freq_short_2)
# AIC(m_freq_short_2)
# gam.check(m_freq_short_2)
# plot(m_freq_short_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_freq_short_gam.rds")
# saveRDS(m_freq_short_2, path)

g_m_freq_short_2 <- ggplot(df_short,
                          aes(viieminutine_lqik, average_peak_frequency, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_freq_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\nfrequency (Hz)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.title.y = element_blank())
g_m_freq_short_2

# path <- here("graphs", "g_freq_short_gam.rds")
# saveRDS(g_m_freq_short_2, path)
################################################################################
# long calls
m_freq_long_1 <- lme(average_peak_frequency ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2 + ot3),
                      data = df_long, method = "ML",
                      correlation = corAR1(form = ~ 1|rott))
# summary(m_freq_long_1)
# path <- here("models", "m_freq_long_lme.rds")
# saveRDS(m_freq_long_1, path)

g_m_freq_long_1 <- ggplot(df_long,
                           aes(viieminutine_lqik, average_peak_frequency, 
                               shape = HE_LE, group = HE_LE,
                               colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_freq_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\nfrequency (Hz)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.text.x = element_blank())
g_m_freq_long_1

# path <- here("graphs", "g_freq_long_lme.rds")
# saveRDS(g_m_freq_long_1, path)
################################################################################
fmla_freq_long_2 <- average_peak_frequency ~ exp_phenotype +
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_freq_long_2 <- bam(fmla_freq_long_2, data = df_long, nthreads = 2,
                      family = Gamma(), method = "REML")
# summary(m_freq_long_2)
# AIC(m_freq_long_2)
# gam.check(m_freq_long_2)
# plot(m_freq_long_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_freq_long_gam.rds")
# saveRDS(m_freq_long_2, path)

g_m_freq_long_2 <- ggplot(df_long,
                           aes(viieminutine_lqik, average_peak_frequency, shape = HE_LE, 
                               group = HE_LE,
                               colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_freq_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\nfrequency (Hz)",
       shape = NULL, linetype = NULL, colour = NULL)
g_m_freq_long_2

# path <- here("graphs", "g_freq_long_gam.rds")
# saveRDS(g_m_freq_long_2, path)
################################################################################
################################################################################
# models of amplitude
# short calls
m_ampl_short_1 <- lme(average_peak_amplitude ~ (ot1 + ot2 + ot3)*HE_LE,
                      random = list(rott =~ ot1 + ot2),
                      data = df_short, method = "REML",
                      correlation = corAR1(form = ~ 1|rott))
summary(m_ampl_short_1)
# path <- here("models", "m_ampl_short_lme.rds")
# saveRDS(m_ampl_short_1, path)

g_m_ampl_short_1 <- ggplot(df_short,
                           aes(viieminutine_lqik, average_peak_amplitude, 
                               shape = HE_LE, group = HE_LE,
                               colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_ampl_short_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\namplitude (dB)",
       shape = NULL, linetype = NULL, colour = NULL) +
  theme(legend.position = c(0.80, 0.90), axis.text.x = element_blank(),
        axis.title.y = element_blank())

g_m_ampl_short_1

# path <- here("graphs", "g_ampl_short_lme.rds")
# saveRDS(g_m_ampl_short_1, path)
################################################################################
fmla_ampl_short_2 <- average_peak_amplitude ~ exp_phenotype +
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_ampl_short_2 <- bam(fmla_ampl_short_2, data = df_short, nthreads = 2,
                      family = gaussian(), method = "REML")
# summary(m_ampl_short_2)
# AIC(m_ampl_short_2)
# gam.check(m_ampl_short_2)
# plot(m_ampl_short_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_ampl_short_gam.rds")
# saveRDS(m_ampl_short_2, path)

g_m_ampl_short_2 <- ggplot(df_short,
                           aes(viieminutine_lqik, average_peak_amplitude, shape = HE_LE, 
                               group = HE_LE,
                               colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_ampl_short_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\namplitude (dB)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.title.y = element_blank())
g_m_ampl_short_2

# path <- here("graphs", "g_ampl_short_gam.rds")
# saveRDS(g_m_ampl_short_2, path)
################################################################################
# long calls
m_ampl_long_1 <- lme(average_peak_amplitude ~ (ot1 + ot2 + ot3)*HE_LE,
                     random = list(rott =~ ot1 + ot2 + ot3),
                     data = df_long, method = "ML",
                     correlation = corAR1(form = ~ 1|rott))
# summary(m_ampl_long_1)
# path <- here("models", "m_ampl_long_lme.rds")
# saveRDS(m_ampl_long_1, path)

g_m_ampl_long_1 <- ggplot(df_long,
                          aes(viieminutine_lqik, average_peak_amplitude, 
                              shape = HE_LE, group = HE_LE,
                              colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_ampl_long_1),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\namplitude (dB)",
       shape = NULL, linetype = NULL, colour = NULL) + 
  theme(axis.text.x = element_blank())
g_m_ampl_long_1

# path <- here("graphs", "g_ampl_long_lme.rds")
# saveRDS(g_m_ampl_long_1, path)
################################################################################
fmla_ampl_long_2 <- average_peak_amplitude ~ exp_phenotype +
  s(viieminutine_lqik, bs = 'cr', k = 6, by = exp_phenotype) +
  s(viieminutine_lqik, bs = 'cr', k = 6) +
  s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're')
m_ampl_long_2 <- bam(fmla_ampl_long_2, data = df_long, nthreads = 2,
                     family = gaussian(), method = "REML")
# summary(m_ampl_long_2)
# AIC(m_ampl_long_2)
# gam.check(m_ampl_long_2)
# plot(m_ampl_long_2, shade = TRUE, pages = 1, scale = 0)
# path <- here("models", "m_ampl_long_gam.rds")
# saveRDS(m_ampl_long_2, path)

g_m_ampl_long_2 <- ggplot(df_long,
                          aes(viieminutine_lqik, average_peak_amplitude, shape = HE_LE, 
                              group = HE_LE,
                              colour = HE_LE)) +
  mean_stats +
  cl_stats +
  stat_summary(aes(y = fitted(m_ampl_long_2),
                   linetype = HE_LE),
               fun = mean, geom = "line", position = pd) + 
  c_scl +
  x_scl +
  labs(x = NULL,
       y = "Average peak\namplitude (dB)",
       shape = NULL, linetype = NULL, colour = NULL)
g_m_ampl_long_2

# path <- here("graphs", "g_ampl_long_gam.rds")
# saveRDS(g_m_ampl_long_2, path)
################################################################################
################################################################################
# plot of the spline shape
m_dur_short_3 <- gam(duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr',
                                                     k = 6, by = exp_phenotype) +
                    s(viieminutine_lqik, bs = 'cr', k = 6) +
                    s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're'),
                    family = gaussian(link = "log"), method = "REML",
                    data = df_short)
b_dur_short_3 <- getViz(m_dur_short_3)
g_dur_short_3 <- plot(sm(b_dur_short_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  coord_cartesian(y = c(-0.3, 0.3)) + 
  x_scl + 
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_dur_short_3 <- g_dur_short_3 %>% gridPrint()

m_dur_long_3 <- gam(duration_ms ~ exp_phenotype + s(viieminutine_lqik, bs = 'cr',
                                                    k = 5, by = exp_phenotype) +
                      s(viieminutine_lqik, bs = 'cr', k = 6) +
                      s(rott, bs = 're') + s(rott, viieminutine_lqik, bs = 're'),
                    data = df_long,
                    family = Gamma(), method = "REML")
b_dur_long_3 <- getViz(m_dur_long_3)
g_dur_long_3 <- plot(sm(b_dur_long_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  x_scl + 
  scale_y_continuous(trans = "reverse") + 
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_dur_long_3 <- g_dur_long_3 %>% gridPrint()

b_freq_short_3 <- getViz(m_freq_short_2)
g_freq_short_3 <- plot(sm(b_freq_short_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  x_scl + 
  scale_y_continuous(trans = "reverse") +
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_freq_short_3 <- g_freq_short_3 %>% gridPrint()

b_freq_long_3 <- getViz(m_freq_long_2)
g_freq_long_3 <- plot(sm(b_freq_long_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  x_scl + 
  scale_y_continuous(trans = "reverse") + 
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_freq_long_3 <- g_freq_long_3 %>% gridPrint()

b_ampl_short_3 <- getViz(m_ampl_short_2)
g_ampl_short_3 <- plot(sm(b_ampl_short_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  x_scl + 
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_ampl_short_3 <- g_ampl_short_3 %>% gridPrint()

b_ampl_long_3 <- getViz(m_ampl_long_2)
g_ampl_long_3 <- plot(sm(b_ampl_long_3, 1)) + l_fitLine(colour = "red", lwd = 0.8) + 
  l_ciLine(level = 0.95, colour = "blue", linetype = 3, lwd = 0.8) + 
  x_scl + 
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        rect = element_rect(fill = "transparent"))
g_ampl_long_3 <- g_ampl_long_3 %>% gridPrint()
################################################################################
# cowplot grids
grid_dur <- plot_grid(g_m_dur_long_1, g_m_dur_short_1, g_m_dur_long_2, g_m_dur_short_2, 
          ncol = 2, labels = "AUTO", scale = 1,
          label_fontfamily = "Times New Roman", label_size = 16, rel_widths = c(1.11, 1))
grid_freq <- plot_grid(g_m_freq_long_1, g_m_freq_short_1, g_m_freq_long_2, g_m_freq_short_2, 
          ncol = 2, labels = "AUTO", 
          label_fontfamily = "Times New Roman", label_size = 16, rel_widths = c(1.17, 1))
grid_ampl <- plot_grid(g_m_ampl_long_1, g_m_ampl_short_1, g_m_ampl_long_2, g_m_ampl_short_2, 
             ncol = 2, labels = "AUTO", 
             label_fontfamily = "Times New Roman", label_size = 16, rel_widths = c(1.17, 1))

# draw insets
grid_dur_inset <- ggdraw(grid_dur) + draw_plot(g_dur_short_3, x = 0.81, y = 0.34, 
                                               width = 0.15, height = 0.15) +
  draw_plot(g_dur_long_3, x = 0.36, y = 0.34, 
            width = 0.15, height = 0.15)

grid_freq_inset <- ggdraw(grid_freq) + draw_plot(g_freq_short_3, x = 0.81, y = 0.12, 
                                               width = 0.15, height = 0.15) +
  draw_plot(g_freq_long_3, x = 0.36, y = 0.12, 
            width = 0.15, height = 0.15)

grid_ampl_inset <- ggdraw(grid_ampl) + draw_plot(g_ampl_short_3, x = 0.81, y = 0.33, 
                                                width = 0.15, height = 0.15) +
  draw_plot(g_ampl_long_3, x = 0.36, y = 0.33, 
            width = 0.15, height = 0.15)
################################################################################
# # save files
# path <- here("Figures", "grid_dur.pdf")
# save_plot(path, grid_dur, base_width = 7.08661, base_height = 4.72441)
# 
# path <- here("Figures", "grid_freq.pdf")
# save_plot(path, grid_freq, base_width = 7.08661, base_height = 4.72441)
# 
# path <- here("Figures", "grid_ampl.pdf")
# save_plot(path, grid_ampl, base_width = 7.08661, base_height = 4.72441)
# 
# path <- here("Figures", "grid_dur_splines.pdf")
# save_plot(path, grid_dur_inset, base_width = 7.08661, base_height = 4.72441)
# 
# path <- here("Figures", "grid_freq_splines.pdf")
# save_plot(path, grid_freq_inset, base_width = 7.08661, base_height = 4.72441)
# 
# path <- here("Figures", "grid_ampl_splines.pdf")
# save_plot(path, grid_ampl_inset, base_width = 7.08661, base_height = 4.72441)
################################################################################
# save files as jpegs
path <- here("Figures", "grid_dur.tiff")
save_plot(path, grid_dur, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)

path <- here("Figures", "grid_freq.tiff")
save_plot(path, grid_freq, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)

path <- here("Figures", "grid_ampl.tiff")
save_plot(path, grid_ampl, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)

path <- here("Figures", "grid_dur_splines.tiff")
save_plot(path, grid_dur_inset, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)

path <- here("Figures", "grid_freq_splines.tiff")
save_plot(path, grid_freq_inset, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)

path <- here("Figures", "grid_ampl_splines.tiff")
save_plot(path, grid_ampl_inset, base_width = 7.08661, base_height = 4.72441, 
          dpi = 600)