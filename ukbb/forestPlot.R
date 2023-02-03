################
# Forest plots #
################
library(ggplot2)

## Results for time at baseline address >=3 years prior to baseline 
cancer_pm2.5 = read.csv(file = "AirPoll_timeAtBaselineBronchus_and_lung_CIs_noInteraction_noPackyears.csv", stringsAsFactors = F, header = T)

# if p value is less than 0.001, replace with p<0.001
cancer_pm2.5$Significance = NA
cancer_pm2.5[cancer_pm2.5$pval < 0.05, "Significance"] = "p<0.05"
cancer_pm2.5[cancer_pm2.5$pval >= 0.05, "Significance"] = "ns"
cancer_pm2.5$pval <- sprintf("%.3f", round(cancer_pm2.5$pval, digits = 3))
cancer_pm2.5[cancer_pm2.5$pval == "0.000", "pval"] <- "<0.001"

## clean coefficient names
cancer_pm2.5$coefficient = gsub(pattern = "Genetic_principal_components", replacement = "gPC", x = cancer_pm2.5$coefficient)
cancer_pm2.5[cancer_pm2.5$coefficient == "Smoking.status.0.0Previous", "coefficient"] = "Previous smoker"
cancer_pm2.5[cancer_pm2.5$coefficient == "Smoking.status.0.0Current", "coefficient"] = "Current smoker"
cancer_pm2.5[cancer_pm2.5$coefficient == "Air.pollution.pm2.5", "coefficient"] = "PM2.5"
cancer_pm2.5[cancer_pm2.5$coefficient == "Sex.0.0Male", "coefficient"] = "Male sex"
cancer_pm2.5[cancer_pm2.5$coefficient == "Tobacco.smoke.exposure.home.0.0", "coefficient"] = "Passive smoking"
cancer_pm2.5[cancer_pm2.5$coefficient == "Educational.attainment>= degree", "coefficient"] = "Education (degree/professional)"
cancer_pm2.5[cancer_pm2.5$coefficient == "Total.household.income.before.tax.0.0>= 31k", "coefficient"] = "Household income (>=£31k)"
cancer_pm2.5[cancer_pm2.5$coefficient == "Packyears.0.0", "coefficient"] = "Packyears"
cancer_pm2.5[cancer_pm2.5$coefficient == "BMI.0.0", "coefficient"] = "BMI"

cancer_pm2.5$coefficient = factor(cancer_pm2.5$coefficient,
                                  levels = cancer_pm2.5[order(cancer_pm2.5$HR, decreasing = F), "coefficient"])

cancer_pm2.5$logHR = log2(cancer_pm2.5$HR)
cancer_pm2.5$logLowerCI = log2(cancer_pm2.5$LowerCI)
cancer_pm2.5$logUpperCI = log2(cancer_pm2.5$UpperCI)

### forest plot with log x-scale
ggplot(data = cancer_pm2.5) +
  geom_point(mapping = aes(x = logHR, y = coefficient, group = Significance, color = Significance), size = 4, shape = 19) +
  geom_errorbarh(aes(xmin = logLowerCI, xmax = logUpperCI, y = coefficient), height = 0.2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate(geom = "text", x = rep(3.6, length(cancer_pm2.5$coefficient)), 
           y = cancer_pm2.5$coefficient, label = cancer_pm2.5$pval) +
  coord_cartesian(xlim = c(-1, 2.8), clip = "off") +
  xlab("log2 Hazard ratio (95% CI)") + ylab("Variables\n") +
  theme(axis.title.x = element_text(size = 16, vjust = 0, face = "bold"),
        axis.text = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"),
        axis.line = element_line(color = "black"), panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(1,6,1,1), "lines"), title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.position = "top") 
dev.off()
