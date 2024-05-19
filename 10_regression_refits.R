# scatter plot of diff
lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                     d_trans = (log(rowMeans(rate_class_b_rate)) - 
                       log(rowMeans(rate_class_a_rate))),
                     d_turns = (log(rowMeans(rate_class_b_turn)) - 
                       log(rowMeans(rate_class_a_turn))))

init_fit = phylolm(d_turns ~ d_trans , data=lm_dat, phy=phy_bb, boot = 1000)
summary(init_fit)

plot(lm_dat)
abline(init_fit, col="red")
abline(a = 0, b = 1, type = 2)

all_fits <- list()
for(j in 1:10000){
  print(j)
  no_to_swap <- sample(1:49, 1)
  to_swap <- sample(1:49, no_to_swap, FALSE)
  d_turns <- d_trans <- c()
  for(i in 1:49){
    if(i %in% to_swap){
      d_trans <- (c(d_trans, log(rowMeans(rate_class_a_rate)[i]) - 
                      log(rowMeans(rate_class_b_rate)[i])))
      d_turns <- (c(d_turns, log(rowMeans(rate_class_a_turn)[i]) - 
                      log(rowMeans(rate_class_b_turn)[i])))
    }else{
      d_trans <- (c(d_trans, log(rowMeans(rate_class_b_rate)[i]) - 
                      log(rowMeans(rate_class_a_rate)[i])))
      d_turns <- (c(d_turns, log(rowMeans(rate_class_b_turn)[i]) - 
                      log(rowMeans(rate_class_a_turn)[i])))
    }
  }
  
  lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                       d_trans = d_trans,
                       d_turns = d_turns)
  fit = phylolm(d_turns ~ d_trans, data=lm_dat, phy=phy_bb, boot = 0)
  # plot(lm_dat, ylim=c(-6, 6), xlim=c(-6, 6))
  # abline(fit)
  # abline(v=0)
  # abline(h=0)
  all_fits[[j]] <- fit
  # Sys.sleep(.2)
}


tmp <- summary(fit)
tmp$coefficients
all_ps <- do.call(rbind, lapply(all_fits, function(x) summary(x)$coefficients[,4]))
all_est <- do.call(rbind, lapply(all_fits, function(x) summary(x)$coefficients[,1]))

colMeans(all_est)
hist(all_est[,2])
length(which(all_ps[,2] <= 0.05))

hist(all_ps[,2])
abline(v = 0.05, col = "red")
