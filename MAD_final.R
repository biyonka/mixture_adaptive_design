
random_seed = 22
set.seed(random_seed)

#2 arms, beta prior, bernoulli outcomes, p1=0.8, p2=0.2
#at each timepoint, estimate probabilitiy arm 1 > arm 2
# use that as plug in for propensity and calculate confidence region and ATE (0.6)
true_p1 = 0.8
true_p0 = 0.7
alpha = 0.05
t=10000
arm1_params = c(1, 1)
arm0_params = c(1, 1)
t_star = t
#set eta according to recommendation from Waubdy-Smith 2023
eta = sqrt(-2*log(alpha) +log(-2*log(alpha)+1))/sqrt(t_star)
delta = 1/seq(1, t)**0.24

compute_TS_Bern_CSs_standard = function(true_p1, true_p0, alpha, t, arm1_params, arm0_params){
  ates_i = c()
  ate = c()
  arm_path = c()
  sig2_hats = c()
  rewards = c()
 # V_is = c()
for (i in seq(1:t)){
  #calculate/estimate P(arm 1 > arm 2)
  if (i == 1){
    p_1 = 0.5
  }else{
    p_1 = sum(rbeta(10000, arm1_params[1], arm1_params[2]) >= rbeta(10000, arm0_params[1], arm0_params[2]))/10000
  }
  p_0 = 1-p_1
  
  #select arm
  if(rbeta(1, arm1_params[1], arm1_params[2]) >= rbeta(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
  arm_path = append(arm_path, arm)
  #observe reward and compute ATE
  if(arm==1){
    y_i = rbinom(1,1, true_p1)
    ate_i = y_i/p_1
    sigma2_hat_i = (y_i)**2/(p_1)**2
  }else{
    y_i = rbinom(1,1, true_p0)
    ate_i = -y_i/p_0
    sigma2_hat_i = (y_i)**2/(p_0)**2
  } 
  rewards = append(rewards, y_i)
  ates_i = append(ates_i, ate_i)
  ate_i = (1/i)*sum(ates_i)
  ate = append(ate, ate_i)
  #compute CS
  sig2_hats = append(sig2_hats, sigma2_hat_i)
  #s_hat = sum(s_hats)
  #print(s_hat)
  #V_i = (1/i)*sqrt((((s_hat * eta**2)+1)/eta**2) * log(((s_hat*eta**2) + 1)/alpha**2))
  #V_is = append(V_is, V_i)
  #update parameters
  if(arm==1){
    arm1_params[1] = arm1_params[1]+y_i
    arm1_params[2] = arm1_params[2]+1-y_i 
  }else{
    arm0_params[1] = arm0_params[1]+ y_i
    arm0_params[2] = arm0_params[2]+1-y_i}
}
  return(list(  'ate' = ate,
                'arm_path' = arm_path,
                'sig_hats' = sig2_hats,
                'rewards' = rewards
         ))
}


compute_CS_eta = function(eta, s_hats, alpha){
  V_is = (1/seq(1, t))*sqrt((((cumsum(s_hats) * eta**2)+1)/eta**2) * log(((cumsum(s_hats) * eta**2) + 1)/alpha**2))
  return(V_is)
}

# 

compute_TS_Bern_CSs_MAD = function(true_p1, true_p0, alpha, t, arm1_params, arm0_params, delta){
  ates_i_MAD = c()
  ate_MAD = c()
  arm_path_MAD = c()
  sig2_hats_MAD = c()
  rewards_MAD = c()
 # V_is_MAD = c()
  for (i in seq(1:t)){
    #calculate/estimate P(arm 1 > arm 2)
    if (i == 1){
      p_1 = 0.5
    }else{
      p_1 = sum(rbeta(10000, arm1_params[1], arm1_params[2]) >= rbeta(10000, arm0_params[1], arm0_params[2]))/10000
    }
    p_0 = 1-p_1
    p_1_MAD = (delta[i])*0.5 + (1-delta[i])*p_1
    p_0_MAD = 1-p_1_MAD 
    #select arm
    #if we go with Bernoulli design
    if(rbinom(1,1,delta[i]) == 1){
      if(rbinom(1,1,0.5)==1){arm = 1}else{arm=0}
    }else{
      if(rbeta(1, arm1_params[1], arm1_params[2]) >= rbeta(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
    }
    arm_path_MAD = append(arm_path_MAD, arm)
    #observe reward and compute ATE
    if(arm==1){
      y_i = rbinom(1,1, true_p1)
      ate_i = y_i/p_1_MAD
      sigma2_hat_i = (y_i)**2/(p_1_MAD)**2
    }else{
      y_i = rbinom(1,1, true_p0)
      ate_i = -y_i/p_0_MAD
      sigma2_hat_i = (y_i)**2/(p_0_MAD)**2
    } 
    rewards_MAD = append(rewards_MAD, y_i)
    ates_i_MAD = append(ates_i_MAD, ate_i)
    ate_i = (1/i)*sum(ates_i_MAD)
    ate_MAD = append(ate_MAD, ate_i)
    #compute CS
    sig2_hats_MAD = append(sig2_hats_MAD, sigma2_hat_i)
    #update parameters
    if(arm==1){
      arm1_params[1] = arm1_params[1]+y_i
      arm1_params[2] = arm1_params[2]+1-y_i 
    }else{
      arm0_params[1] = arm0_params[1]+ y_i
      arm0_params[2] = arm0_params[2]+1-y_i}
  }
  
  return(list(  'ate' = ate_MAD,
                'arm_path' = arm_path_MAD,
                'sig_hats' = sig2_hats_MAD,
                'rewards' = rewards_MAD
  ))
}


N=100
#compute CS's for both methods many times and show all CS paths
testing = lapply(seq(1, N), function(i){
  standard_results = compute_TS_Bern_CSs_standard(true_p1, true_p0, alpha, t, arm1_params, arm0_params)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  sig2_hats_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  V_i_standard = compute_CS_eta(eta, sig2_hats_standard, alpha)
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'sig_hats' = sig2_hats_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
  })

#compute CS's for MAD for both methods many times and show all CS paths
testing_MAD = lapply(seq(1, N), function(i){
  standard_results = compute_TS_Bern_CSs_MAD(true_p1, true_p0, alpha, t, arm1_params, arm0_params, delta)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  sig2_hats_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  V_i_standard = compute_CS_eta(eta, sig2_hats_standard, alpha)
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'sig_hats' = sig2_hats_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
})


#saveRDS(object = testing_MAD, file = paste0(data_path,'/MAD_CS_TS.rds'))


#plot average reward between two approaches with Bernoulli design as baseline
#compute CS's for MAD for both methods many times and show all CS paths
testing_bernoulli_des = lapply(seq(1, N), function(i){
  standard_results = compute_TS_Bern_CSs_MAD(true_p1, true_p0, alpha, t, arm1_params, arm0_params, delta = 1/seq(1, t)**0)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  sig2_hats_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  V_i_standard = compute_CS_eta(eta, sig2_hats_standard, alpha)
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'sig_hats' = sig2_hats_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
})


#saveRDS(object = testing_bernoulli_des, file = paste0(data_path,'/bernoulli_CS_TS.rds'))



#plotting sequence of CS's
#generate dataframe of rewards from each method
#reward over the 400 timesteps, across the 100 iterations
bern_ate = as.vector(matrix(sapply(seq(1, N), function(i){testing_bernoulli_des[[i]]$ate}), byrow = FALSE, nrow = 1))
MAD_ate = as.vector(matrix(sapply(seq(1, N), function(i){testing_MAD[[i]]$ate}), byrow = FALSE, nrow = 1))
standard_ate = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$ate}), byrow = FALSE, nrow = 1))
bern_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing_bernoulli_des[[i]]$CS}), byrow = FALSE, nrow = 1))
MAD_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing_MAD[[i]]$CS}), byrow = FALSE, nrow = 1))
standard_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$CS}), byrow = FALSE, nrow = 1))
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
CS_ate_df = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                             'ate' = bern_ate, 'Vt' = bern_Vt, 'replicate' = iteration),
                  data.frame('Method' = c('MAD'), 'Time' = time, 
                             'ate' = MAD_ate, 'Vt' = MAD_Vt,  'replicate' = iteration),
                  data.frame('Method' = c('Standard Design'), 'Time' = time, 
                             'ate' = standard_ate, 'Vt' = standard_Vt,  'replicate' = iteration))

#calculate number of standard design CS's which still include 0 by end of time
(sum(between(rep(0, length(CS_ate_df[CS_ate_df$Method == 'Standard Design' & CS_ate_df$Time == 10000, ]$ate)), CS_ate_df[CS_ate_df$Method == 'Standard Design' & CS_ate_df$Time == 10000, ]$ate -CS_ate_df[CS_ate_df$Method == 'Standard Design'& CS_ate_df$Time == 10000, ]$Vt,
CS_ate_df[CS_ate_df$Method == 'Standard Design'& CS_ate_df$Time == 10000, ]$ate +CS_ate_df[CS_ate_df$Method == 'Standard Design' & CS_ate_df$Time == 10000, ]$Vt), na.rm = TRUE))/length(CS_ate_df[CS_ate_df$Method == 'Standard Design' & CS_ate_df$Time == 10000, ]$ate)

#plot confidence sequences across the replicates
CS_plot_all_ham_log10 = ggplot(CS_ate_df, aes(x = Time, group=replicate, fill=Method)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.05, show.legend=FALSE)+ #ylim(c(-100, 100))+
 # scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
  #scale_x_continuous(expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+ xlab('Time (log10 Scale)')+
  geom_hline(aes(yintercept=true_p1-true_p0, linetype = 'True ATE'), alpha = 0.5) +scale_x_log10(expand = c(0, 0)) +
   scale_linetype_manual(name = "", values = c(1)) + #+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 

CS_plot_all_ham_lines = ggplot(CS_ate_df, aes(x = Time, group=replicate, color=Method)) + geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
 # geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.05, show.legend=FALSE)+ #ylim(c(-100, 100))+
  # scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
  #scale_x_continuous(expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+xlab('Time (log10 Scale)')+
  geom_hline(aes(yintercept=true_p1-true_p0, linetype = 'True ATE'), alpha = 0.5) +scale_x_log10(expand = c(0, 0)) +   scale_linetype_manual(name = "", values = c(1)) +
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 


CS_ate_df['upper'] = sapply(seq(1, 3000000), function(i){
  if(CS_ate_df$ate[i] + CS_ate_df$Vt[i] > 1 | is.na(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] + CS_ate_df$Vt[i])){1}else{CS_ate_df$ate[i] + CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
CS_ate_df['lower'] = sapply(seq(1, 3000000), function(i){if(CS_ate_df$ate[i] - CS_ate_df$Vt[i] < -1 | is.na(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] - CS_ate_df$Vt[i])){-1}else{CS_ate_df$ate[i] - CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)

CS_plot_all_ham_zoom = ggplot(CS_ate_df, aes(x = Time, group=replicate, fill=Method)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=upper, ymin = lower), alpha=0.05, show.legend=FALSE)+ #ylim(c(0, 1))+
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
 # scale_x_continuous(expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=true_p1-true_p0, linetype = 'True ATE')) +scale_x_log10(expand = c(0, 0)) + xlab('Time (log10 Scale)')+
  scale_linetype_manual(name = "", values = c(1)) +#+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 

ggsave(filename = paste0('CS_plot_all_ham_zoom.png'),
       plot = CS_plot_all_ham_zoom , path = '~/Documents/Research/MAD/plots/', bg = 'transparent')


#TODO: plot ratio of bernoulli width and others
ratio_plotting = rbind(data.frame('Method' = c('MAD'), 'Time' = time, 
                                  'Ratio' = MAD_Vt/bern_Vt, 'replicate' = iteration),
                       data.frame('Method' = c('Standard Design'), 'Time' = time, 
                                  'Ratio' = standard_Vt/bern_Vt, 'replicate' = iteration))

ratio_plot = ggplot(ratio_plotting, aes(x = Time, group=replicate, color=Method)) + geom_line(aes(y=Ratio)) +#+ geom_line(aes(y=ate-Vt)) +
 # geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.3, show.legend=FALSE)+ #ylim(c(-100, 100))+
  # scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE), scales = 'free')+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +scale_x_log10() #+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)

ggsave(filename = paste0('ratio_plot.png'),
       plot = ratio_plot, path = '~/Documents/Research/MAD/plots/', bg = 'white')

ggsave(filename = paste0('ratio_plot.eps'),
       plot = ratio_plot , path = '~/Documents/Research/MAD/plots/', bg = 'white')




#plot empirical coverage of both methods: empirical proportion of times each respective confidence sequence covers 
#the true treatment effect for each t, for our choice of alpha.

#coverage over time: proportion of replicates where the true treatment effect is within the CS
CS_ate_df['coverage'] = sapply(seq(1, length(CS_ate_df$ate)), function(i){if(between(true_p1-true_p0, CS_ate_df$lower[i], CS_ate_df$upper[i])){1}else{0}})
cov_df = CS_ate_df%>% group_by(Method, Time) %>% summarize(avg_cov = mean(coverage), sd = sd(coverage)/sqrt(N))

coverage = ggplot(cov_df, aes(x = Time, y = avg_cov, color=Method, fill=Method)) + geom_line() +
  geom_ribbon(aes(ymax=avg_cov + 2*sd, ymin = avg_cov-2*sd), alpha=0.4, show.legend=FALSE, colour=NA)+ theme_minimal() + ylab('Coverage')+
  geom_hline(aes(yintercept=1-alpha), linetype = 'dashed', alpha = 0.5) +
  scale_x_log10(expand=c(0, 0)) + 
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
      legend.title = element_text(size = 25))+ xlab('Time (log10 Scale)') #+ 
  #scale_linetype_manual(name = expression(paste(1-alpha)), values = c(2))#+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)

#generate dataframe of rewards from each method
#reward over the 400 timesteps, across the 100 iterations
bern_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing_bernoulli_des[[i]]$reward}), byrow = FALSE, nrow = 1))
MAD_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing_MAD[[i]]$reward}), byrow = FALSE, nrow = 1))
standard_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$reward}), byrow = FALSE, nrow = 1))
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
rewards = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                           'reward' = bern_rewards, 'replicate' = iteration),
                data.frame('Method' = c('MAD'), 'Time' = time, 
                           'reward' = MAD_rewards, 'replicate' = iteration),
                data.frame('Method' = c('Standard Design'), 'Time' = time, 
                           'reward' = standard_rewards, 'replicate' = iteration))




rewards_plotting = rewards %>% group_by(Method, Time) %>% summarize(avg_reward = mean(reward))%>% mutate(cum_avg_reward = cumsum(avg_reward), time_avg_reward = (1/Time)*cumsum(avg_reward))
#rewards %>% group_by(Method, Time) %>% summarize(se_reward = sd(reward))

rewards_ses = rewards %>% group_by(Method, replicate) %>% mutate(cumsum_reward_per_it = cumsum(reward),
                                                                 time_avg_reward_per_it = cumsum(reward)/Time)%>% 
  group_by(Method, Time) %>% summarize(se_cum_reward = sd(cumsum_reward_per_it), se_time_avg_reward = sd(time_avg_reward_per_it))

rewards_plotting['se_cum_reward'] = rewards_ses$se_cum_reward/sqrt(N)
rewards_plotting['se_time_avg_reward'] = rewards_ses$se_time_avg_reward/sqrt(N)

#average cumulative reward with SEs
ham_total_rewards_plot_log10=ggplot(rewards_plotting, aes(x = Time, y = cum_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=cum_avg_reward + 2*se_cum_reward, ymin = cum_avg_reward - 2*se_cum_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Cumulative Average Reward')  +scale_x_log10()+  
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + xlab('Time (log10 Scale)')

#time average reward with SEs
ham_time_avg_rewards_plot_log10=ggplot(rewards_plotting, aes(x = Time, y = time_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=time_avg_reward + 2*se_time_avg_reward, ymin = time_avg_reward - 2*se_time_avg_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Time Average Reward') +scale_x_log10() +  
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + xlab('Time (log10 Scale)')


ham_total_rewards_plot=ggplot(rewards_plotting, aes(x = Time, y = cum_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=cum_avg_reward + 2*se_cum_reward, ymin = cum_avg_reward - 2*se_cum_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Cumulative Average Reward')  + # theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) 
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 
#time average reward with SEs
ham_time_avg_rewards_plot=ggplot(rewards_plotting, aes(x = Time, y = time_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=time_avg_reward + 2*se_time_avg_reward, ymin = time_avg_reward - 2*se_time_avg_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Time Average Reward') + # theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) 
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 
#make plot for slide showing howard and Ham alone
CS_standard = CS_ate_df[CS_ate_df$Method == 'Standard Design',]

#count up number of standard CS's which by end of time, has a CS within [-1, 1]
counting = 0
for(rep in 1:100){
  rep_only = CS_standard[CS_standard$replicate == rep & CS_standard$Time == 10000,]
  if (sum((rep_only$ate + rep_only$Vt) < 1 & (rep_only$ate - rep_only$Vt) > -1, na.rm = TRUE) > 0){counting = counting + 1}
  #  if (sum(rep_only$upper < 1 & rep_only$lower > -1, na.rm = TRUE) > 0){counting = counting + 1}
}



CS_standard_alone = ggplot(CS_standard, aes(x = Time, group=replicate)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.3, show.legend=FALSE, fill = '#619cff')+ #ylim(c(-100, 100))+
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+#+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +#   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  scale_x_log10() +  theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) + xlab('Time (log10 Scale)')


CS_standard_alone_zoom = ggplot(CS_standard, aes(x = Time, group=replicate)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=upper, ymin= lower), alpha=0.03, show.legend=FALSE, fill = '#619cff')+ #ylim(c(-100, 100))+
  # scale_y_continuous(limits = c(-20, 20), expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+#+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +#   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  scale_x_log10() +  theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) + xlab('Time (log10 Scale)')


ggsave(filename = paste0('CS_standard_alone.png'),
       plot = CS_standard_alone, path = '~/Documents/Research/MAD/plots/', bg = 'transparent')


ggsave(filename = paste0('CS_standard_alone_zoom.png'),
       plot = CS_standard_alone_zoom, path = '~/Documents/Research/MAD/plots/', bg = 'transparent')





#plot stopping times of all methods
CS_ate_df = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                             'ate' = bern_ate, 'Vt' = bern_Vt, 'replicate' = iteration),
                  data.frame('Method' = c('MAD'), 'Time' = time, 
                             'ate' = MAD_ate, 'Vt' = MAD_Vt,  'replicate' = iteration),
                  data.frame('Method' = c('Standard'), 'Time' = time, 
                             'ate' = standard_ate, 'Vt' = standard_Vt,  'replicate' = iteration))

standard_stoppings = sapply(seq(1, N), function(i){
  subset = CS_ate_df[CS_ate_df$Method == 'Standard' & CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
  })

MAD_stoppings = sapply(seq(1, N), function(i){
  subset = CS_ate_df[CS_ate_df$Method == 'MAD' & CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
})

bernoulli_stoppings = sapply(seq(1, N), function(i){
  subset = CS_ate_df[CS_ate_df$Method == 'Bernoulli Design' & CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
})

#
mean(standard_stoppings == -1)

stoppings_df = rbind(data.frame(data.frame('Method' = c('Bernoulli Design'), 'replicate' = seq(1, N), 'stopping' = bernoulli_stoppings)),
data.frame(data.frame('Method' = c('MAD'), 'replicate' = seq(1, N), 'stopping' = MAD_stoppings)),
data.frame(data.frame('Method' = c('Standard Design'), 'replicate' = seq(1, N), 'stopping' = standard_stoppings))) 


stoppings_plot_log10 = ggplot(stoppings_df) + geom_histogram(aes(stopping, fill = Method), show.legend = FALSE) +theme_minimal() +
  facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE)) + xlab('Stopping Time') + scale_x_log10() +xlim(c(0, 10000))+
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + xlab('Time')

ggsave(filename = paste0('stoppings_plot_log10.png'),
       plot = stoppings_plot_log10, path = '~/Documents/Research/MAD/plots/', bg = 'transparent')

mean(standard_stoppings == -1)
mean(bernoulli_stoppings == -1)
mean(MAD_stoppings == -1)

