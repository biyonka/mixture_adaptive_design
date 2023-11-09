
random_seed = 22
set.seed(random_seed)
N=100
#2 arms, beta prior, bernoulli outcomes, p1=0.8, p2=0.2
true_p1 = 0.8
true_p0 = 0.4
alpha = 0.05
t=10000
arm1_params = c(1, 1)
arm0_params = c(1, 1)
delta = pmax(1/seq(1, t)**0.24, 1/100)

#parameters for CS in Corollary 2 of Howard paper
#suggested choosing eta close to 1
eta1 = 1.1
m = 10
s=1.1 #s>1
#use h suggested from paper
h =function(k,s){return(eta1^(s*k)/(1-eta1^(-s)))}
l0 = 1
l = function(v,m,eta1,alpha,l0){log(h(log(v/m, base = eta1), s)) + log(l0/alpha)}
#calcuated gamma exponential mixture boundary
calculate_boundary =function(v, l, c, eta1, alpha, l0, m){
  k1 = (eta1^0.25 + eta1^(-0.25))/sqrt(2)
  k2 = (sqrt(eta1)+1)/2
  return(sqrt(k1^2 * max(v, m) *l(max(v, m),m,eta1,alpha,l0) + k2^2 * c^2 * l(max(v, m),m,eta1,alpha,l0)^2) + k2*c*l(max(v, m),m,eta1,alpha,l0))}


compute_TS_howard_MAD = function(true_p1, true_p0, alpha, t, arm1_params, arm0_params, C, delta){
  ates_i = c()
  ate = c()
  arm_path = c()
  sig2_hats = c()
  rewards = c()
  CSs = c()
  for (i in seq(1:t)){
    #calculate/estimate P(arm 1 > arm 2)
    if (i == 1){
      p_1 = 0.5
    }else{
      #estimate probability one beta is bigger than another
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
    arm_path = append(arm_path, arm)
    if(arm==1){
      y_i = rbinom(1,1, true_p1)
      ate_i = y_i/p_1_MAD
    }else{
      y_i = rbinom(1,1, true_p0)
      ate_i = -y_i/p_0_MAD
    }
    rewards = append(rewards, y_i)
    ates_i = append(ates_i, ate_i)
    ate_i = (1/i)*sum(ates_i)
    Vi = sum((ates_i-(1/i)*sum(ates_i))**2)
    ate = append(ate, ate_i)
    CS = calculate_boundary(Vi, l, C, eta, alpha, l0, m)/i
    CSs = append(CSs, CS)
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
                'CSs' = CSs,
                'rewards' = rewards
  ))
}

compute_TS_howard = function(true_p1, true_p0, alpha, t, arm1_params, arm0_params, C){
  ates_i = c()
  ate = c()
  arm_path = c()
  sig2_hats = c()
  rewards = c()
  CSs = c()
  for (i in seq(1:t)){
    #calculate/estimate P(arm 1 > arm 2)
    if (i == 1){
      p_1 = 0.5
    }else{
      p_1 = sum(rbeta(10000, arm1_params[1], arm1_params[2]) >= rbeta(10000, arm0_params[1], arm0_params[2]))/10000
    }
    p_0 = 1-p_1
    #select arm
    #select arm
    if(rbeta(1, arm1_params[1], arm1_params[2]) >= rbeta(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
    arm_path = append(arm_path, arm)
    #observe reward and compute ATE
    if(arm==1){
      y_i = rbinom(1,1, true_p1)
      ate_i = y_i/p_1
    }else{
      y_i = rbinom(1,1, true_p0)
      ate_i = -y_i/p_0
    }
    rewards = append(rewards, y_i)
    ates_i = append(ates_i, ate_i)
    ate_i = (1/i)*sum(ates_i)
    Vi = sum((ates_i-(1/i)*sum(ates_i))**2)
    ate = append(ate, ate_i)
    CS = calculate_boundary(Vi, l, C, eta, alpha, l0, m)/i
    CSs = append(CSs, CS)
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
                'CSs' = CSs,
                'rewards' = rewards
  ))
}


#' compute standard howard results with heuristic p_min
pmin = 1/t
C=2/pmin
#compute CS's for both methods many times and show all CS paths
howard_standard = lapply(seq(1, N), function(i){
  standard_results = compute_TS_howard(true_p1, true_p0, alpha, t, arm1_params, arm0_params, C)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  V_i_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
})
#' 
#' set p_min for MAD
pmin = 0.5/100
C=2/pmin
howard_MAD = lapply(seq(1, N), function(i){
  standard_results = compute_TS_howard_MAD(true_p1, true_p0, alpha, t, arm1_params, arm0_params, C, delta)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  V_i_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
})

p_min = 0.5
C = 2/p_min
howard_Bernoulli = lapply(seq(1, N), function(i){
  standard_results = compute_TS_howard_MAD(true_p1, true_p0, alpha, t, arm1_params, arm0_params, C, delta= 1/seq(1, t)**0)
  ate_standard = standard_results[[1]]
  arm_path_standard = standard_results[[2]]
  V_i_standard = standard_results[[3]]
  rewards_standard = standard_results[[4]]
  print(i)
  return(list(
    'ate' = ate_standard,
    'arm_path' = arm_path_standard,
    'rewards' = rewards_standard,
    'CS' = V_i_standard
  ))
})
#' 
#' 
#' 
how_MAD_ate = as.vector(matrix(sapply(seq(1, N), function(i){howard_MAD[[i]]$ate}), byrow = FALSE, nrow = 1))
howard_ate = as.vector(matrix(sapply(seq(1, N), function(i){howard_standard[[i]]$ate}), byrow = FALSE, nrow = 1))
how_bern_ate = as.vector(matrix(sapply(seq(1, N), function(i){howard_Bernoulli[[i]]$ate}), byrow = FALSE, nrow = 1))
how_MAD_Vt = as.vector(matrix(sapply(seq(1, N), function(i){howard_MAD[[i]]$CS}), byrow = FALSE, nrow = 1))
howard_Vt = as.vector(matrix(sapply(seq(1, N), function(i){howard_standard[[i]]$CS}), byrow = FALSE, nrow = 1))
how_bern_Vt = as.vector(matrix(sapply(seq(1, N), function(i){howard_Bernoulli[[i]]$CS}), byrow = FALSE, nrow = 1))
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
howard_CS_ate_df = rbind(data.frame('Method' = c('Standard Design  (Howard et al., 2021)'), 'Time' = time,
                             'ate' = howard_ate, 'Vt' = howard_Vt, 'replicate' = iteration),
                  data.frame('Method' = c('MAD  (Howard et al., 2021)'), 'Time' = time,
                             'ate' = how_MAD_ate, 'Vt' = how_MAD_Vt,  'replicate' = iteration),
                  data.frame('Method' = c('Bernoulli Design  (Howard et al., 2021)'), 'Time' = time,
                             'ate' = how_bern_ate, 'Vt' = how_bern_Vt,  'replicate' = iteration))


howard_CS_ate_df['upper'] = sapply(seq(1, 3000000), function(i){
  if(howard_CS_ate_df$ate[i] + howard_CS_ate_df$Vt[i] > 1 |is.infinite(howard_CS_ate_df$ate[i] + howard_CS_ate_df$Vt[i])| is.nan(howard_CS_ate_df$ate[i] + howard_CS_ate_df$Vt[i])|is.na(howard_CS_ate_df$ate[i] + howard_CS_ate_df$Vt[i])){1}else{howard_CS_ate_df$ate[i] + howard_CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
howard_CS_ate_df['lower'] = sapply(seq(1, 3000000), function(i){if(howard_CS_ate_df$ate[i] - howard_CS_ate_df$Vt[i] < -1|is.infinite(howard_CS_ate_df$ate[i] - howard_CS_ate_df$Vt[i] < -1)|is.nan(howard_CS_ate_df$ate[i] - howard_CS_ate_df$Vt[i] < -1)|is.na(howard_CS_ate_df$ate[i] - howard_CS_ate_df$Vt[i] < -1)){-1}else{howard_CS_ate_df$ate[i] - howard_CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)

#plot confidence sequences across the replicates
howard_CS_plot = ggplot(howard_CS_ate_df, aes(x = Time, group=replicate, fill=Method)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.1, show.legend=FALSE)+ #ylim(c(-100, 100))+
  # scale_y_continuous(limits = c(-20, 20), expand = c(0, 0))+
  scale_y_continuous(labels = label_comma()) +
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +#   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  scale_x_log10() +
  theme(axis.title =element_text(size=32), axis.text = element_text(size = 30)) 


howard_CS_plot_zoom = ggplot(howard_CS_ate_df, aes(x = Time, group=replicate, fill=Method)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=upper, ymin = lower), alpha=0.05, show.legend=FALSE)+ #ylim(c(0, 1))+
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
  # scale_x_continuous(expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=true_p1-true_p0, linetype = 'True ATE')) +scale_x_log10(expand = c(0, 0)) + xlab('Time (log10 Scale)')+
  scale_linetype_manual(name = "", values = c(1)) +#+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 



standard_stoppings_howard = sapply(seq(1, N), function(i){
  subset = howard_CS_ate_df[howard_CS_ate_df$Method == 'Standard Design' & howard_CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
})

MAD_stoppings_howard = sapply(seq(1, N), function(i){
  subset = howard_CS_ate_df[CS_ate_df$Method == 'MAD' & howard_CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
})

bernoulli_stoppings_howard = sapply(seq(1, N), function(i){
  subset = howard_CS_ate_df[howard_CS_ate_df$Method == 'Bernoulli Design' & howard_CS_ate_df$replicate == i,]
  x = which(!between(rep(0, length(subset$ate)), subset$ate -subset$Vt, subset$ate + subset$Vt))
  if(length(x) == 0){-1}else{x[1]}
})

mean(standard_stoppings_howard == -1)
mean(bernoulli_stoppings_howard == -1)
mean(MAD_stoppings_howard == -1)


howard_standard_df = howard_CS_ate_df[howard_CS_ate_df$Method == 'Standard Design',]
#count number of times Howard CS with standard design excluded -1, 1
counting = 0
for(rep in 1:100){
  rep_only = howard_standard_df[howard_standard_df$replicate == rep,]
  if (sum((rep_only$ate + rep_only$Vt) < 1 & (rep_only$ate - rep_only$Vt) > -1, na.rm = TRUE) > 0){counting = counting + 1}
  #  if (sum(rep_only$upper < 1 & rep_only$lower > -1, na.rm = TRUE) > 0){counting = counting + 1}
}


#generate dataframe of rewards from each method
#reward over the 400 timesteps, across the 100 iterations
how_MAD_rewards= as.vector(matrix(sapply(seq(1, N), function(i){howard_MAD[[i]]$rewards}), byrow = FALSE, nrow = 1))
howard_rewards = as.vector(matrix(sapply(seq(1, N), function(i){howard_standard[[i]]$rewards}), byrow = FALSE, nrow = 1))
how_bern_rewards = as.vector(matrix(sapply(seq(1, N), function(i){howard_Bernoulli[[i]]$rewards}), byrow = FALSE, nrow = 1))
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
howard_rewards = rbind(data.frame('Method' = c('Standard Design'), 'Time' = time, 
                                  'reward' = howard_rewards, 'replicate' = iteration),
                       data.frame('Method' = c('MAD'), 'Time' = time, 
                                  'reward' = how_MAD_rewards, 'replicate' = iteration),
                       data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                                  'reward' = how_bern_rewards, 'replicate' = iteration))




howard_rewards_plotting = howard_rewards %>% group_by(Method, Time) %>% summarize(avg_reward = mean(reward))%>% mutate(cum_avg_reward = cumsum(avg_reward), time_avg_reward = (1/Time)*cumsum(avg_reward))
#rewards %>% group_by(Method, Time) %>% summarize(se_reward = sd(reward))

howard_rewards_ses = howard_rewards %>% group_by(Method, replicate) %>% mutate(cumsum_reward_per_it = cumsum(reward),
                                                                               time_avg_reward_per_it = cumsum(reward)/Time)%>% 
  group_by(Method, Time) %>% summarize(se_cum_reward = sd(cumsum_reward_per_it), se_time_avg_reward = sd(time_avg_reward_per_it))

howard_rewards_plotting['se_cum_reward'] = howard_rewards_ses$se_cum_reward/sqrt(N)
howard_rewards_plotting['se_time_avg_reward'] = howard_rewards_ses$se_time_avg_reward/sqrt(N)

#average cumulative reward with SEs
ggplot(howard_rewards_plotting, aes(x = Time, y = cum_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=cum_avg_reward + 2*se_cum_reward, ymin = cum_avg_reward - 2*se_cum_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Cumulative Average Reward') +  theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) 

#time average reward with SEs
howard_rewards_plot_log10 = ggplot(howard_rewards_plotting, aes(x = Time, y = time_avg_reward, color = Method, fill = Method)) + geom_line() +
  geom_ribbon(aes(ymax=time_avg_reward + 2*se_time_avg_reward, ymin = time_avg_reward - 2*se_time_avg_reward), alpha=0.3, colour = NA)+ 
  theme_minimal() +  ylab('Time Average Reward')  +
  scale_x_log10()+
  theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) 



#make plot for slide showing howard and Ham alone
howard_CS_standard = howard_CS_ate_df[howard_CS_ate_df$Method == 'Standard Design  (Howard et al., 2021)',]

library(scales)
howard_CS_plot_alone = ggplot(howard_CS_standard, aes(x = Time, group=replicate)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=ate+Vt, ymin = ate-Vt), alpha=0.3, show.legend=FALSE, fill = '#619cff')+ #ylim(c(-100, 100))+
  scale_y_continuous(labels = label_comma()) +# scale_y_continuous(limits = c(-20, 20), expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+#+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +#   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  scale_x_log10() + theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) + xlab('Time (log10 Scale)')



ggsave(filename = paste0('howard_CS_plot_alone.png'),
       plot = howard_CS_plot_alone, path = '~/Documents/Research/MAD/plots/', bg = 'transparent')



howard_CS_plot_alone_zoom = ggplot(howard_CS_standard, aes(x = Time, group=replicate)) +# geom_line(aes(y=ate+Vt)) + geom_line(aes(y=ate-Vt)) +
  geom_ribbon(aes(ymax=upper, ymin = lower), alpha=0.05, show.legend=FALSE, fill = '#619cff')+ #ylim(c(-100, 100))+
  # scale_y_continuous(limits = c(-20, 20), expand = c(0, 0))+
  theme_minimal() +  ylab('CS')+#+facet_wrap(~Method, ncol = 1, dir = 'h', labeller = labeller(.multi_line = FALSE))+
  geom_hline(aes(yintercept=0.2), alpha = 0.5) +#   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
  scale_x_log10() + theme(aspect.ratio = 1, axis.title =element_text(size=32), axis.text = element_text(size = 30)) + xlab('Time (log10 Scale)')


ggsave(filename = paste0('howard_CS_plot_alone_zoom.png'),
       plot = howard_CS_plot_alone_zoom, path = '~/Documents/Research/MAD/plots/', bg = 'transparent')




