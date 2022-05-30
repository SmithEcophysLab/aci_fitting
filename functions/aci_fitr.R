# aci_fitr script
  ## code to fit aci curves
  ## author: nick smith
  ## contact: nick.smith@ttu.edu

  ## notes
  # calculation is based on partial pressure
  # must input elevation for atmospheric pressure calculation
  # the thought is that you can source this script and then have all the functions you'll need
  # might need to adjust the previous thought? not sure of the best practice here
  # check TPU calculation
  # add in Rd estimation

  ## functions
    ### calculate atmospheric pressure
    calc_patm <- function(z) {
      
      kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
      kTo <- 298.15   # base temperature, K (Prentice, unpublished)
      kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
      kG <- 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
      kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
      kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
      
      patm <- kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
      
      patm
      
    }

    ### calculate gammastar
    calc_gammastar <- function(tleaf, z) {
      
      patm <- calc_patm(z)
      rat <- calc_patm(z) / calc_patm(0)
      
      gammastar25 <- 4.332 * rat  # Pa
      Hgm <- 37830 # J mol-1
      R <- 8.314        # J K-1 mol-1
    
      tleaf_k <- 273.15+ tleaf
      
      gammastar <- gammastar25*exp((Hgm/R)*(1/298.15-1/tleaf_k))
      
      gammastar
      
    }
    
    ### calculate km
    calc_km <- function(tleaf, z) {
      
      patm <- calc_patm(z) 
      rat <- patm / calc_patm(0)
      
      R <- 8.314        
      O2 <- 2.09476e5      
      Kc25 <- 41.03 * rat 
      Ko25 <- 28210 * rat 
      Hkc <- 79430  
      Hko <- 36380 
      
      tleaf_k <- 273.15 + tleaf
      
      Kc_pa <-Kc25 * exp(Hkc * ((tleaf_k - 298.15) / (298.15 * R * tleaf_k)))
      Ko_pa <-Ko25* exp(Hko * ((tleaf_k - 298.15) / (298.15 * R * tleaf_k)))
      
      O2_pa <- O2 * (1e-6) * patm 
      
      km <- Kc_pa * (1 + O2_pa/Ko_pa)
      
      km 
      
    }

    
    ### aci_fitr
    aci_fitr <- function(a = a, ci = ci, par = par, tleaf = tleaf, 
                         ci_transition_jmax, ci_transition_tpu, z = 0, data){
      
      patm <- calc_patm(z)
      tleaf <- mean(data$tleaf, na.rm = T)
      par <- mean(data$par, na.rm = T)
      gammastar <- calc_gammastar(tleaf, z)
      km = calc_km(tleaf, z)
      data$ci_pa <- data$ci * 1e-6 * patm
      
      ac_data <- subset(data, ci <= ci_transition_jmax)
      aj_data <- subset(data, ci > ci_transition_jmax & ci <= ci_transition_tpu)
      at_data <- subset(data, ci > ci_transition_tpu)
      
      vcmax_nls <- nls(a ~ vcmax * ((ci_pa - gammastar) * (1 / (ci_pa + km))) - rd, 
                       start = list(vcmax = 100, rd = 2), data = ac_data)
      vcmax <- coef(summary(vcmax_nls))[1, 1]
      vcmax_se <- coef(summary(vcmax_nls))[1, 2]
      rd <- coef(summary(vcmax_nls))[1, 2]
      rd_se <- coef(summary(vcmax_nls))[2, 2]
      
      j_nls <- nls(a ~ (j/4) * ((ci_pa - gammastar) * (1 / (ci_pa + 2 * gammastar))) - rd, 
                   start = list(j = 50), data = aj_data)
      j <- coef(summary(j_nls))[1]
      j_se <- coef(summary(j_nls))[2]
      
      # jmax_nls <- nls()
      # jmax <- coef(summary(jmax_nls))[1]
      # jmax_se <- coef(summary(jmax_nls))[2]
      
      # tpu <- mean(at_data$a - rd, na.rm = T) ## check formula!
      # tpu_se <- sd(at_data$a - rd, na.rm = T) / sqrt(length(mean(at_data$a, na.rm = T)))
      
      fit_params <- as.data.frame(cbind(vcmax, vcmax_se, j, j_se, 
                                        # tpu, tpu_se, 
                                        rd, rd_se,
                                     patm, z, gammastar, km, tleaf, par))
      colnames(fit_params) <- c('vcmax', 'vcmax_se', 'j', 'j_se', 
                                # 'tpu', 'tpu_se', 
                                'rd', 'rd_se',
                             'patm', 'z', 'gammastar', 'km', 'tleaf', 'par')
      
      modeled_ac <- fit_params$vcmax * 
        (data$ci_pa - fit_params$gammastar) * (1 / (data$ci_pa + fit_params$km))
      modeled_aj <- (fit_params$j / 4) * 
        (data$ci_pa - fit_params$gammastar) * (1 / (data$ci_pa + 2 * fit_params$gammastar))
      # modeled_at <- fit_params$tpu
      
      modeled_results <- data.frame(cbind(modeled_ac, modeled_aj))
      modeled_results$modeled_anet <- pmin(modeled_ac, modeled_aj) - fit_params$rd
      modeled_results$ci_pa <- data$ci_pa
      
      fit_statistics <- cor(data$a, modeled_results$modeled_anet)**2
      
      results <- list(fit_params, modeled_results, fit_statistics)

    }

    
    
    ### plotting function
    aci_plotr <- function(data, # the raw data
                          fit_params, # the fitted parameters
                          modeled_results # the modeled results
                          ){
      
      plot(modeled_results$modeled_anet ~ data$ci, type = 'l', col = 'grey', lwd = 3, ylim = c(-5, 20))
      points(data$a ~ data$ci, data = data, pch = 1, col= 'black', cex = 1)
      
      ac <- modeled_results$modeled_ac - fit_params$rd
      lines(ac ~ data$ci, col = 'red', lty = 2, lwd =2)
      
      aj <- modeled_results$modeled_aj  - fit_params$rd
      lines(aj ~ data$ci, col = 'blue', lty = 2, lwd = 2)
      
      # at <- fit_params$tpu  - fit_params$rd
      # abline(a = at, b = 0, col = 'orange', lty = 2, lwd = 2)
      
    }

