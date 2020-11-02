library(papaja)
library(scales)
library(broom)
library(broom.mixed)
library(here)
library(tidyverse)
library(lme4)
library(kableExtra)
library(lmerTest)
library(merTools)
library(ggsci)
options(options(scipen=999))
theme_set(theme_apa(base_size = 14) + theme(legend.position = "bottom"))

path <- here("raw_data", "NPS_madalad_Denisile.xlsx")
df <- readxl::read_xlsx(path, sheet = 1, col_names = TRUE)
df <- df %>% mutate_at(1:6, factor) %>% mutate(kestvuskategooria = 
             factor(kestvuskategooria, levels = 1:5, labels = c("<12 ms",
             "12-99 ms", "100-299 ms", "300-999 ms", "≥1000 ms")),
             keskmise_frequency_kategooria = factor(keskmise_frequency_kategooria,
             levels = 1:5, labels = c("<25000 Hz", "25000-29950 Hz",
             "30000-34950 Hz", "35000-39950 Hz", "40000-44950 Hz")),
             duration_cat = factor(ifelse(duration_ms < 300, "short", "long"))) %>% 
  dplyr::arrange(rott, viieminutine_lqik, peak_no_lqigus) 
# path <- here("mod_data", "NPS_data.rds")
# saveRDS(df, path)

df_summary_coarse <- df %>% group_by(HE_LE, rott, viieminutine_lqik,
                                     duration_cat) %>% 
  summarise(peak_frequency = mean(average_peak_frequency), 
            peak_amplitude = mean(average_peak_amplitude), 
            duration = mean(duration_ms), vocalisations = n())

df_short <- df %>% dplyr::filter(kestvuskategooria %in% 
                                   c("<12 ms", "12-99 ms", "100-299 ms")) %>% 
  dplyr::mutate(kestvuskategooria = fct_drop(kestvuskategooria))

df_long <- df %>% dplyr::filter(kestvuskategooria %in% 
                                  c("300-999 ms", "≥1000 ms")) %>% 
  dplyr::mutate(kestvuskategooria = fct_drop(kestvuskategooria))
################################################################################
g_short_avg_peak_freq_vs_amp <- df_short %>% ggplot(aes(x = average_peak_frequency,
                          y = average_peak_amplitude, colour = kestvuskategooria)) + 
  geom_point(alpha = 0.4) + scale_color_jco() + facet_grid(cols = vars(HE_LE))
g_short_avg_peak_freq_vs_amp
################################################################################
g_short_avg_peak_freq_vs_dur <- df_short %>% ggplot(aes(x = average_peak_frequency,
                                                        y = duration_ms, 
                                                        colour = kestvuskategooria)) + 
  geom_point(alpha = 0.4) + scale_color_jco() + facet_grid(cols = vars(HE_LE))
g_short_avg_peak_freq_vs_dur
################################################################################
g_avg_peak_freq_vs_dur <- df %>% ggplot(aes(y = average_peak_frequency,
                                                        x = duration_ms,
                                                        colour = duration_cat)) + 
  geom_point(alpha = 1, size = 0.3, stroke = 0) + 
  scale_color_jco() + facet_grid(cols = vars(HE_LE)) +
  labs(colour = "duration category", x = "duration (ms)", 
       y = "mean peak frequency (Hz)") +
  guides(colour = guide_legend(override.aes = list(size = 4)))
g_avg_peak_freq_vs_dur
path <- here("Figures", "g_avg_peak_freq_vs_dur.pdf")
ggsave(path, g_avg_peak_freq_vs_dur, width=6, height=4.5, dpi=300)
################################################################################
g_avg_peak_amp_vs_dur <- df %>% ggplot(aes(y = average_peak_amplitude,
                                            x = duration_ms,
                                            colour = duration_cat)) + 
  geom_point(alpha = 1, size = 0.3, stroke = 0) + 
  scale_color_jco() + facet_grid(cols = vars(HE_LE)) +
  labs(colour = "duration category", x = "duration (ms)", 
       y = "mean peak amplitude") +
  guides(colour = guide_legend(override.aes = list(size = 4)))
g_avg_peak_amp_vs_dur
path <- here("Figures", "g_avg_peak_amp_vs_dur.pdf")
ggsave(path, g_avg_peak_amp_vs_dur, width=6, height=4.5, dpi=300)
################################################################################

g_short_avg_peak_amp_vs_dur <- df_short %>% ggplot(aes(x = average_peak_amplitude,
                                                        y = duration_ms, 
                                                       colour = kestvuskategooria)) + 
  geom_point(alpha = 0.4) + scale_color_jco() + facet_grid(cols = vars(HE_LE))
g_short_avg_peak_amp_vs_dur
################################################################################
g_vocal_by_dur <- df_summary_coarse %>% 
  ggplot(aes(x = viieminutine_lqik, y = vocalisations, color = rott, 
             linetype = HE_LE)) + 
  geom_line(aes(group = rott)) + facet_grid(cols = vars(duration_cat))
g_vocal_by_dur
################################################################################
g_avg_peak_freq_vs_dur_id_le <- df %>% dplyr::filter(HE_LE == "LE") %>% 
  ggplot(aes(x = duration_ms, y = average_peak_frequency, colour = rott)) + 
  geom_point(alpha = 0.7, size = 0.2, stroke = 0) + scale_color_ucscgb() +
  labs(colour = "rat", x = "duration (ms)", 
       y = "mean peak frequency (Hz)", title = "LE") +
  guides(colour = guide_legend(override.aes = list(size = 4)))
g_avg_peak_freq_vs_dur_id_le
path <- here("Figures", "g_avg_peak_freq_vs_dur_id_le.pdf")
ggsave(path, g_avg_peak_freq_vs_dur_id_le, width=4, height=4, dpi=300)
################################################################################
g_avg_peak_freq_vs_dur_id_he <- df %>% dplyr::filter(HE_LE == "HE") %>% 
  ggplot(aes(x = duration_ms, y = average_peak_frequency, colour = rott)) + 
  geom_point(alpha = 0.7, size = 0.2, stroke = 0) + scale_color_ucscgb() +
  labs(colour = "rat", x = "duration (ms)", 
       y = "mean peak frequency (Hz)", title = "HE") +
  guides(colour = guide_legend(override.aes = list(size = 4)))
g_avg_peak_freq_vs_dur_id_he
path <- here("Figures", "g_avg_peak_freq_vs_dur_id_he.pdf")
ggsave(path, g_avg_peak_freq_vs_dur_id_he, width=4, height=4, dpi=300)
################################################################################
g_avg_peak_freq_vs_dur_id_he_faceted <- df %>% dplyr::filter(HE_LE == "HE") %>% 
  ggplot(aes(x = duration_ms, y = average_peak_frequency)) + 
  geom_point(alpha = 1, size = 0.2, stroke = 0) +
  geom_vline(xintercept = 300, colour = "coral", linetype = 3) +
  labs(x = "duration (ms)", 
       y = "mean peak frequency (Hz)", title = "HE") + 
  facet_wrap(vars(rott), ncol = 3)
g_avg_peak_freq_vs_dur_id_he_faceted 
path <- here("Figures", "g_avg_peak_freq_vs_dur_id_he_faceted.pdf")
ggsave(path, g_avg_peak_freq_vs_dur_id_he_faceted, width=9, height=9, dpi=300)
################################################################################
g_avg_peak_freq_vs_dur_id_le_faceted <- df %>% dplyr::filter(HE_LE == "LE") %>% 
  ggplot(aes(x = duration_ms, y = average_peak_frequency)) + 
  geom_point(alpha = 1, size = 0.2, stroke = 0) + 
  geom_vline(xintercept = 300, colour = "coral", linetype = 3) +
  labs(x = "duration (ms)", 
       y = "mean peak frequency (Hz)", title = "LE") + 
  facet_wrap(vars(rott), ncol = 3)
g_avg_peak_freq_vs_dur_id_le_faceted 
path <- here("Figures", "g_avg_peak_freq_vs_dur_id_le_faceted.pdf")
ggsave(path, g_avg_peak_freq_vs_dur_id_le_faceted, width=9, height=12, dpi=300)
################################################################################
################################################################################