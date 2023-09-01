#Testing data poisson took ~30 mins to complete for 2K itns
# 1: There were 2689 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems

#Clean output for Nos-GAL4 and NA driver + construct but Act5c and Tubp very wide variance (posterior credible interval basically (-3, +3))

#Separating into Nos-GAL4 or NA (including construct only orientations)

traceplot(model_output, c("base_reproduction", "theta","lp__"))
traceplot(model_output, c("beta_vals"))

summ_model_output <- summary(model_output)

sum(summ_model_output$summary[,10] > 1.01)

plot(summ_model_output$summary[1:10,])

plot(model_output, pars = paste0("beta_vals[",1:6,"]"))
plot(model_output, pars = paste0("beta_vals[",7:8,"]"))

plot(model_output, pars = c("base_reproduction",paste0("beta_vals[",1:8,"]")))

plot(model_output, pars = "alpha")
plot(model_output, pars = "etas")

plot(model_output, pars = paste0("etas[",sort(
                                 unique(eta_fx$eta_fctr[
                                   which(eta_fx$driver_id %in% c(2,4) & 
                                           eta_fx$eta_fctr > 1)] - 1)
                                 ),
                                 "]")) +
  geom_vline(xintercept = 0, lwd = 2, lty = 2)
