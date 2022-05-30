# test_aci_fitr.R
  ## test the aci_fitr and associated functions

  ## source the function
  source('functions/aci_fitr.R')
  
  ## test the functions
    ### create data frame
    test_a <- c(0, 4, 8, 12, 13.5, 15, 16, 17, 18, 18.5, 18.6, 18.8, 19, 19, 19)
    test_ci <- seq(30, 1200, 1200/15)
    test_df <- as.data.frame(cbind(test_a, test_ci, 25, 1800))
    colnames(test_df) <- c('a', 'ci', 'tleaf', 'par')
    
    ### plot data
    plot(a ~ci, data = test_df)
    
    ### test aci_fitr.R
    test_results <- aci_fitr(data = test_df, ci_transition_jmax = 300, ci_transition_tpu = 1100)
    test_results
    
    ### test ac-_plotr
    aci_plotr(data = test_df, fit_params = test_results[[1]], modeled_results = test_results[[2]])
    
    ### other plots?
    plot(test_df$a, test_results[[2]]$modeled_anet)
    abline(0, 1)
  
  