#########################################################################
# Check non-PH assumptions by plotting output of coxzph() and KM curves #
#########################################################################
library(mice)
library(survival)
library(survminer)

## ukbb cancer before baseline - for alcohol 
ukbb_cancer_prebaseline = read.csv(header = T, stringsAsFactors = F, file = "UKBB_decoded_cancers_before_baseline.csv")

## read Non-PH vars
non_ph_vars = readRDS(file = "nonPH_vars.rds")
non_ph_cancers = unique(unlist(lapply(strsplit(non_ph_vars, split = "-imp"), FUN = function(x) {x[1]})))

binary_vars = c("Smoking.status.0.0", "Sex.0.0", "Total.household.income.before.tax.0.0", "Educational.attainment")
mod_coefs = c("Smoking.status.0.0", "Packyears", "Tobacco.smoke.exposure.home.0.0",
              "Sex.0.0Male", "Age", "BMI.0.0", "Total.household.income.before.tax.0.0", "Educational.attainment",                
              "Air.pollution.pm2.5", "Genetic_principal_components.0.1", "Genetic_principal_components.0.2",
              "Genetic_principal_components.0.3", "Genetic_principal_components.0.4", "Genetic_principal_components.0.5",               
              "Genetic_principal_components.0.6", "Genetic_principal_components.0.7", "Genetic_principal_components.0.8",               
              "Genetic_principal_components.0.9" , "Genetic_principal_components.0.10", "Genetic_principal_components.0.11",               
              "Genetic_principal_components.0.12", "Genetic_principal_components.0.13", "Genetic_principal_components.0.14",              
              "Genetic_principal_components.0.15")

for (i in 1:length(non_ph_cancers)) {
  curr_region = non_ph_cancers[i]
  
  test_alcohol = F
  if (curr_region %in% c("Lip_oral_cavity_pharynx", "Larynx")) {
    test_alcohol = T
    orig_form_str = "Surv(time, curr_cancer_binary) ~ Smoking.status.0.0 + Packyears.0.0 + Tobacco.smoke.exposure.home.0.0 + Sex.0.0 + 
                Age + BMI.0.0 + Total.household.income.before.tax.0.0 + Educational.attainment +  
                Air.pollution.pm2.5 + Average.weekly.red.wine.intake.0.0 + Average.weekly.champagne.plus.white.wine.intake.0.0 +
                Average.weekly.beer.plus.cider.intake.0.0 + Average.weekly.spirits.intake.0.0 + Average.weekly.fortified.wine.intake.0.0 + 
                Genetic_principal_components.0.1 + Genetic_principal_components.0.2 + Genetic_principal_components.0.3 + 
                Genetic_principal_components.0.4 + Genetic_principal_components.0.5 + Genetic_principal_components.0.6 +  
                Genetic_principal_components.0.7 + Genetic_principal_components.0.8 + Genetic_principal_components.0.9 + 
                Genetic_principal_components.10 + Genetic_principal_components.11 + Genetic_principal_components.12 + 
                Genetic_principal_components.13 + Genetic_principal_components.14 + Genetic_principal_components.15"
    
    ## create data frame to hold VIF
    vif_df = data.frame(matrix(NA, nrow = 30, ncol = 15))
    colnames(vif_df) = paste0("Imp", 1:15)
    vif_df$Region = curr_region
    rownames(vif_df) = c("Smoking.status.0.0Current", "Smoking.status.0.0Previous", "Packyears", "Tobacco.smoke.exposure.home.0.0", 
                         "Sex.0.0Male", "Age", "BMI.0.0", "Total.household.income.before.tax.0.0>= 31k", "Educational.attainment>= degree",                
                         "Air.pollution.pm2.5", "Average.weekly.red.wine.intake.0.0", "Average.weekly.champagne.plus.white.wine.intake.0.0",
                         "Average.weekly.beer.plus.cider.intake.0.0", "Average.weekly.spirits.intake.0.0", "Average.weekly.fortified.wine.intake.0.0", 
                         "Genetic_principal_components.0.1", "Genetic_principal_components.0.2",
                         "Genetic_principal_components.0.3", "Genetic_principal_components.0.4", "Genetic_principal_components.0.5",               
                         "Genetic_principal_components.0.6", "Genetic_principal_components.0.7", "Genetic_principal_components.0.8",               
                         "Genetic_principal_components.0.9" , "Genetic_principal_components.0.10", "Genetic_principal_components.0.11",               
                         "Genetic_principal_components.0.12", "Genetic_principal_components.0.13", "Genetic_principal_components.0.14",              
                         "Genetic_principal_components.0.15")
  } else {
    orig_form_str = "Surv(time, curr_cancer_binary) ~ Smoking.status.0.0 + Packyears.0.0 + Tobacco.smoke.exposure.home.0.0 + Sex.0.0 + 
                Age + BMI.0.0 + Total.household.income.before.tax.0.0 + Educational.attainment + Air.pollution.pm2.5 +
                Genetic_principal_components.0.1 + Genetic_principal_components.0.2 + Genetic_principal_components.0.3 + 
                Genetic_principal_components.0.4 + Genetic_principal_components.0.5 + Genetic_principal_components.0.6 +  
                Genetic_principal_components.0.7 + Genetic_principal_components.0.8 + Genetic_principal_components.0.9 + 
                Genetic_principal_components.10 + Genetic_principal_components.11 + Genetic_principal_components.12 + 
                Genetic_principal_components.13 + Genetic_principal_components.14 + Genetic_principal_components.15"
                
    ## create data frame to hold VIF
    vif_df = data.frame(matrix(NA, nrow = 25, ncol = 15))
    colnames(vif_df) = paste0("Imp", 1:15)
    vif_df$Region = curr_region
    rownames(vif_df) = c("Smoking.status.0.0Current", "Smoking.status.0.0Previous", "Packyears", "Tobacco.smoke.exposure.home.0.0",
                         "Sex.0.0Male", "Age", "BMI.0.0", "Total.household.income.before.tax.0.0>= 31k", "Educational.attainment>= degree",                
                         "Air.pollution.pm2.5", "Genetic_principal_components.0.1", "Genetic_principal_components.0.2",
                         "Genetic_principal_components.0.3", "Genetic_principal_components.0.4", "Genetic_principal_components.0.5",               
                         "Genetic_principal_components.0.6", "Genetic_principal_components.0.7", "Genetic_principal_components.0.8",               
                         "Genetic_principal_components.0.9" , "Genetic_principal_components.0.10", "Genetic_principal_components.0.11",               
                         "Genetic_principal_components.0.12", "Genetic_principal_components.0.13", "Genetic_principal_components.0.14",              
                         "Genetic_principal_components.0.15")
  }
  
  if (curr_region == "Brain") {curr_region = "GBM"}
  print(paste0("######## current cancer = ", curr_region, " ################"))
 
  ## read in imputed dataset object
  imp_path = paste0("imp_cox_analysis_baselineAge_incl_gPCs_", curr_region, "_outcomes.rds")
  curr_cancer_data = readRDS(file = imp_path)
  
  ## iterate through each imputation which has allegedly non-ph vars and confirm visually - loop needs to be run interactively by changing the index variable to covariate of interest
  for (j in 1:15) { 
    ukbb_clean = complete(curr_cancer_data, j)
    print(paste0("###### ", curr_region, ": imp = ", j, " ############"))
    
    if (test_alcohol) {
      ## if testing alcohol ##
      ukbb_clean = merge(x = ukbb_clean, by = "Patient_ID",
                         y = ukbb_cancer_prebaseline[, c("Patient_ID", grep("Average.weekly", colnames(ukbb_cancer_prebaseline), value = T))])
      ukbb_clean[, "Average.weekly.intake.of.other.alcoholic.drinks.0.0"] <- NULL      
      ukbb_clean$Average.weekly.beer.plus.cider.intake.0.0 = as.numeric(ukbb_clean$Average.weekly.beer.plus.cider.intake.0.0)
      ukbb_clean$Average.weekly.champagne.plus.white.wine.intake.0.0 = as.numeric(ukbb_clean$Average.weekly.champagne.plus.white.wine.intake.0.0)
      ukbb_clean$Average.weekly.red.wine.intake.0.0 = as.numeric(ukbb_clean$Average.weekly.red.wine.intake.0.0)
      ukbb_clean$Average.weekly.fortified.wine.intake.0.0 = as.numeric(ukbb_clean$Average.weekly.fortified.wine.intake.0.0)
      ukbb_clean$Average.weekly.spirits.intake.0.0 = as.numeric(ukbb_clean$Average.weekly.spirits.intake.0.0)
      
      # Remove NAs
      index <- unique(c(which(is.na(ukbb_clean$Average.weekly.beer.plus.cider.intake.0.0)),
                        which(is.na(ukbb_clean$Average.weekly.champagne.plus.white.wine.intake.0.0)),
                        which(is.na(ukbb_clean$Average.weekly.red.wine.intake.0.0)),
                        which(is.na(ukbb_clean$Average.weekly.fortified.wine.intake.0.0)),
                        which(is.na(ukbb_clean$Average.weekly.spirits.intake.0.0))))
      ukbb_clean = ukbb_clean[-index, ]
    }
    ## rename columns to match with variable names
    colnames(ukbb_clean)[colnames(ukbb_clean) == paste0(curr_region, "_binary")] <- "curr_cancer_binary"
    colnames(ukbb_clean)[colnames(ukbb_clean) == paste0(curr_region, "_time")] <- "time"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.10"] <- "Genetic_principal_components.10"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.11"] <- "Genetic_principal_components.11"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.12"] <- "Genetic_principal_components.12"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.13"] <- "Genetic_principal_components.13"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.14"] <- "Genetic_principal_components.14"
    colnames(ukbb_clean)[colnames(ukbb_clean) == "Genetic_principal_components.0.15"] <- "Genetic_principal_components.15"
    
    # Initial model without any time-dependent variables to check proportional hazards assumption
    form_str = orig_form_str
    model = coxph(formula = as.formula(form_str),data = ukbb_clean, na.action = na.omit, x = T)
    
    ## Test cox proportional hazards (PH) assumption
    ph_test = cox.zph(model)
    
    #### confirm results from test - need to repeat manually for each variable of interest by changing index value
    index = 21
    png(paste0("coxzph_", curr_region, "_", mod_coefs[index], "_imp", j, ".png"), height = 666, width = 818)
    plot(ph_test[index], lwd = 2, col = "red")    
    abline(0, 0, col = "black", lty = 3, lwd = 2)
    abline(#h = model$coef[names(model$coef) == "Sex.0.0Male"],
      #h = model$coef[names(model$coef) == "Age"],
      #h = model$coef[names(model$coef) == "Educational.attainment>= degree"],                 
      #h = model$coef[names(model$coef) == "Smoking.status.0.0Previous"],
      #h = model$coef[names(model$coef) == "Smoking.status.0.0Current"],
      #h = model$coef[names(model$coef) == "BMI.0.0"],
      #h = model$coef[names(model$coef) == "Packyears.0.0"],
      #h = model$coef[names(model$coef) == "Air.pollution.pm2.5"],
      h = model$coef[names(model$coef) == "Genetic_principal_components.12"],
      #h = model$coef[names(model$coef) == "Total.household.income.before.tax.0.0>= 31k"],
      #h = model$coef[names(model$coef) == "Tobacco.smoke.exposure.home.0.0"],
      col = "blue", lty = 2, lwd = 2)
    legend("bottomright",
           legend = c("Reference line for null effect", "Average hazard over time", "Time-varying hazard"),
           lty = c(3, 2, 1), col = c("black", "blue", "red"))
    dev.off()
  }
}

#### Code to plot KM curves - e.g. for sex
ukbb_clean$sex_km = NA
ukbb_clean[ukbb_clean$Sex.0.0 == "Female", "sex_km"] = 0
ukbb_clean[ukbb_clean$Sex.0.0 == "Male", "sex_km"] = 1

fit = survfit(Surv(ukbb_clean$time, ukbb_clean$curr_cancer_binary) ~ ukbb_clean$sex_km)
summary(fit)

## base R graphics - change ylim to better visualise KM curves as needed
plot(fit, lty = c("solid", "dashed"), col = c("black", "grey", "red"),
     xlab = "Survival Time In Days", ylab = "Survival Probabilities",
     ylim = c(0.9, 1.00), 
     main = paste0(curr_region, ": Sex: imp", j))
legend("bottomright", c("Female", "Male"), lty = c("solid", "dashed"), col = c("black", "red"))
