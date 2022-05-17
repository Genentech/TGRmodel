library(ggplot2)
library(ggrepel)
library(gridExtra)
library(TGRmodel)

PK_para = data.frame(k_a = 2, k_e = 0.15, MW = 450, VF = 1)
GR_para = data.frame(GR_inf = -0.5, GEC50 = 0.5, h_GR = 2)
Schedule = "QD"
Duration_1 = 7 #days
Duration_2 = 14 #days

Dose_1 = 10 #mg/kg
Dose_2 = 1 #mg/kg

conc_profile_1 = PK_to_conc_profile(PK_para, Dose_1, Schedule, Duration_1)
conc_1 = data.frame("Time" = conc_profile_1$Time, "Conc" = conc_profile_1$Conc(conc_profile_1$Time))
relk_1 = relk_over_time(GR_para, conc_profile_1)
TGR_predicted_1 = GR_inVitro_integration(conc_profile_1, GR_para, int_method = 'mean')

conc_profile_2 = PK_to_conc_profile(PK_para, Dose_2, Schedule, Duration_2, conc_1$Conc[nrow(conc_1)])
conc_2 = data.frame("Time" = conc_profile_2$Time, "Conc" = conc_profile_2$Conc(conc_profile_2$Time))
relk_2 = relk_over_time(GR_para, conc_profile_2)
TGR_predicted_2 = GR_inVitro_integration(conc_profile_2, GR_para, int_method = 'mean')

conc_2 = conc_2[2:nrow(conc_2),]
conc_2$Time = conc_2$Time + Duration_1

df_conc = rbind(conc_1, conc_2)

relk_avg = mean(relk_fct(df_conc$Conc, GR_para))
TGR_eff = (2^relk_avg)-1

k0 = 1/14
growth_untreat = 2^(k0*df_conc$Time)
growth_1 = 2^(k0*mean(relk_1$relk)*relk_1$Time)
growth_2 = growth_1[[length(growth_1)]]*2^(k0*mean(relk_2$relk)*relk_2$Time[-1])

growth_unt_plot = data.frame(Time = df_conc$Time, rel_vol = growth_untreat)
growth_plot = data.frame(Time = df_conc$Time, rel_vol = c(growth_1, growth_2))

final_plot = ggplot(mapping = aes(x = Time)) +
  geom_line(growth_unt_plot, mapping = aes(y = rel_vol), color = "red", size = 1) +
  geom_line(growth_plot, mapping = aes(y = rel_vol), color = "black", size = 1) +
  geom_line(data = df_conc, mapping = aes(y = Conc/10)) +
  labs(x = "Time (days)", y = "log2 Relative Volume") +
  theme_bw() + 
  scale_y_continuous("Relative tumor volume", sec.axis = sec_axis(~ . *10, name = "Serum Concentration [µM]")) +
  coord_cartesian(ylim = c(0,3.5)) +
  annotate('text', x = Duration_1+Duration_2, y = growth_plot$rel_vol[nrow(growth_plot)] + .05, hjust = 1, label = sprintf('TGR = %.2f', TGR_eff)) + 
  annotate('rect', xmin = 0, xmax = Duration_1, ymin = 0, ymax = -.07, size = 5, fill = 'blue', alpha = .2) +
  annotate('rect', xmin = Duration_1, xmax = Duration_1+Duration_2, ymin = 0, ymax = -.07, size = 5, fill = 'green', alpha = .2) +
  annotate('text', x = Duration_1/2, y = -.01, hjust = .5, vjust = 1, label = sprintf('%.1f mg/kg', Dose_1)) +
  annotate('text', x = Duration_1 + Duration_2/2, y = -.01, hjust = .5, vjust = 1, label = sprintf('%.1f mg/kg', Dose_2))
  #theme(text = element_text(family = "sans"), axis.title.x = element_text(size = rel(3)), axis.title.y = element_text(size = rel(3)), axis.text.x=element_text(size=rel(3)), axis.text.y=element_text(size=rel(3)), plot.title = element_text(hjust = 0.5), legend.title = element_text(color = 'black'))

print(final_plot)


