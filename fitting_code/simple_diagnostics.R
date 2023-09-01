traceplot(model_output, c("base_reproduction", "theta","lp__"))
traceplot(model_output, c("beta_vals"))



summ_model_output <- summary(model_output)

sum(summ_model_output$summary[,10] > 1.01)

plot(summ_model_output$summary[1:10,])

plot(model_output, pars = paste0("beta_vals[",1:6,"]"))
plot(model_output, pars = paste0("beta_vals[",7:8,"]"))

plot(model_output, pars = c("base_reproduction",paste0("beta_vals[",1:8,"]")))
