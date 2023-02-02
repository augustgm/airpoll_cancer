############################################
# Merge ukbb data with radon data from PHE #
############################################
library(dplyr)

## Read ukbb data 
ukbb_decoded <- read.csv("UKBB_decoded_cancer_reg_hosp_diag_data.csv", stringsAsFactors = F, header = T)

## Read in cancer pre-baseline data
ukbb_cancer_prebaseline = read.csv(header = T, stringsAsFactors = F, file = "UKBB_decoded_cancers_before_baseline.csv")
ukbb_decoded = merge(x = ukbb_decoded, 
                     y = ukbb_cancer_prebaseline[, c("Patient_ID", "cancer_before_baseline_binary")], 
                     by = "Patient_ID")

## Read coordinate data
ukbb_coord = read.csv(file = "ukbb_address_data.csv", stringsAsFactors = F, header = T)
ukbb_decoded = merge(x = ukbb_decoded, y = ukbb_coord, by = "Patient_ID")
rm(ukbb_coord)

## Remove unnecessary relevant columns only
ukbb_decoded = ukbb_decoded[, !(colnames(ukbb_decoded) %in% grep("Diagnoses.ICD10", colnames(ukbb_decoded), value = T))]
ukbb_decoded = ukbb_decoded[, !(colnames(ukbb_decoded) %in% grep("Date.of.first.IP.diagnosis", colnames(ukbb_decoded), value = T))]
ukbb_decoded = ukbb_decoded[, !(colnames(ukbb_decoded) %in% grep("Cancer.record.origin", colnames(ukbb_decoded), value = T))]
ukbb_decoded = ukbb_decoded[, !(colnames(ukbb_decoded) %in% grep("Cancer.record.format", colnames(ukbb_decoded), value = T))]

## Read PHE air poll data
air_poll = read.csv(file = "all_pollutants.csv", stringsAsFactors = F, header = T)

### Merge by rounded coordinates to nearest 1000 (for all_pollutants.csv)
air_poll$eastings_round = signif(air_poll$eastings, digits = 3)
air_poll$northings_round = signif(air_poll$northings, digits = 3)

## Exclude those with cancers diagnosed before baseline
ukbb_decoded = ukbb_decoded[ukbb_decoded$cancer_before_baseline_binary == 0, ]

## Match based on string concatenation
ukbb_decoded$east_north_concat = paste(ukbb_decoded$Home.location.assessment.east.coord.rounded.0.0, ukbb_decoded$Home.location.assessment.north.coord.rounded.0.0, sep = "-")
air_poll$east_north_concat = paste(air_poll$eastings_round, air_poll$northings_round, sep = "-")

## Summarise 
find_mode <- function(x) { 
  u <- unique(x)
  tab <- tabulate(match(x, u))
  modal_value = u[tab == max(tab)]
  
  ## edit to return maximum Radon potential class if multiple values are modal
  if (length(modal_value) > 1) {
    output = max(modal_value)
  } else {
    output = modal_value
  }
  return(output)}

air_poll_sum = air_poll %>%
  group_by(east_north_concat) %>%
  summarize(mode_Radon = find_mode(Radon))

## perform merge
desired_cols = c("Patient_ID", "Ever.smoked.0.0", "Smoking.status.0.0", "Tobacco.smoke.exposure.home.0.0", "Sex.0.0",
                 "BMI.0.0", "Total.household.income.before.tax.0.0", grep("Educational.attainment.0", colnames(ukbb_decoded), value = T),
                 "PM2.5.0.0", "Average.weekly.red.wine.intake.0.0", "Average.weekly.champagne.plus.white.wine.intake.0.0",
                 "Average.weekly.beer.plus.cider.intake.0.0", "Average.weekly.spirits.intake.0.0", "Average.weekly.fortified.wine.intake.0.0",  
                 grep("Genetic.principal.components", colnames(ukbb_decoded), value = T)[1:15],
                 "Home.location.assessment.east.coord.rounded.0.0", "Home.location.assessment.north.coord.rounded.0.0", "east_north_concat")

ukbb_poll = merge(x = ukbb_decoded[, desired_cols], 
                  y = air_poll_sum, by = "east_north_concat")

write.csv(ukbb_poll, row.names = F,
          file = "ukbb_phe_air_poll_merge_rounded_coords.csv")
