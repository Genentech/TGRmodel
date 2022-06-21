library(ggplot2)
library(ggrepel)
library(gridExtra)
library(TGRmodel)


PK_para <- data.frame(k_a = 1.5, k_e = 0.2, MW = 450, VF = 1)
GR_para <- data.frame(GR_inf = -0.5, GEC50 = 0.3, h_GR = 1.5)

Schedule <- "QD"
Duration <- 21

Doses <- c(0.2, .5, 2, 5, 15)
TGR_predicted <- array(NA, length(Doses), dimnames = list(Doses))
dose_plots <- list()

k0 = 1/10

for (id in seq_along(Doses)) {

  conc_profile <- PK_to_conc_profile(PK_para, Doses[id], Schedule, Duration)
  df_conc <- data.frame("Time" = conc_profile$Time, "Conc" = conc_profile$Conc(conc_profile$Time))
  relk <- relk_over_time(GR_para, conc_profile)
  TGR_predicted[id] <- GR_inVitro_integration(conc_profile, GR_para, int_method = 'mean')

  growth_plot <- data.frame(Time = df_conc$Time, rel_vol = 2^(k0*mean(relk$relk)*df_conc$Time))
  growth_unt_plot <- data.frame(Time = df_conc$Time, rel_vol = 2^(k0*df_conc$Time))
  
  dose_plots[[id]] <- ggplot(mapping = aes(x = Time)) +
    geom_line(data = df_conc, mapping = aes(y = Conc/5), color = 'cyan', size = .5) +
    geom_line(growth_unt_plot, mapping = aes(y = rel_vol), color = "red", size = 1) +
    geom_line(growth_plot, mapping = aes(y = rel_vol), color = "black", size = 1) +
    labs(x = "Time (days)", y = "Relative Volume") +
    ggtitle(sprintf('Dose = %.1f mg/kg, TGR = %.2f',  Doses[id], TGR_predicted[id])) +
    theme_bw() + 
    scale_y_continuous("Relative tumor volume", sec.axis = sec_axis(~ . *5, name = "Serum Concentration [µM]")) +
    coord_cartesian(ylim = c(0,4))
}

print(grid.arrange(grobs = c(dose_plots, list(
  ggplot(data.frame(Dose = Doses, TGR = TGR_predicted),aes(Dose, TGR)) +
    geom_hline(yintercept = 1, alpha = .5) +
    geom_hline(yintercept = 0, alpha = .5) +
    geom_point() +
    scale_x_log10() +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = "Dose (mg/kg)", y = "TGR") +
    theme_bw() +
    ggtitle('TGR Dose-response curve')
))))
