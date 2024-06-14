library(dplyr)
library(ggplot2)
random_seed = 21
set.seed(random_seed)

#2 arms, beta prior, bernoulli outcomes, p1=0.8, p2=0.2
#at each timepoint, estimate probabilitiy arm 1 > arm 2
# use that as plug in for propensity and calculate confidence region and ATE (0.6)
N=100
true_mu1 = 0.5
true_mu0 = 0.5
mus = c(true_mu0, true_mu1)
sigma=2 #1,5,10,20,30
alpha = 0.05
t=10000
reward = 'Bern'
alg = 'UCB'
sigma0 = 1
mu0 = 1
arm1_params = c(mu0, sigma0) #initialize priors on reward distribution
arm0_params = c(mu0, sigma0) #initialize priors on reward distribution
t_star = t
#set eta according to recommendation from Waubdy-Smith 2023
eta = sqrt(-2*log(alpha) +log(-2*log(alpha)+1))/sqrt(t_star)
delta = 1/seq(1, t)**0.24 
#delta=pmax(1/seq(1, t)**0.24, 0.2) #for clipping
#delta = 1/seq(1, t)**0.2 #for t dist and cauchy things

compute_TS_CSs_standard = function(mus, alpha, t, arm1_params, arm0_params, reward='Norm'){
  ates_i = c()
  ate = c()
  arm_path = c()
  sig2_hats = c()
  rewards = c()
  for (i in seq(1:t)){
    #calculate/estimate P(arm 1 > arm 2)
    if (i == 1){
      p_1 = 1
    }else{
      if(reward %in% c('Norm', 'Cauchy', 'tdist')){
        p_1 = sum(rnorm(10000, arm1_params[1], arm1_params[2]) >= rnorm(10000, arm0_params[1], arm0_params[2]))/10000
      }else if(reward == 'Bern'){
        p_1 = sum(rbeta(10000, arm1_params[1], arm1_params[2]) >= rbeta(10000, arm0_params[1], arm0_params[2]))/10000
      } 
    }
    p_0 = 1-p_1
    
    #select arm
    if(reward %in% c('Cauchy', 'Norm', 'tdist')){
      if(rnorm(1, arm1_params[1], arm1_params[2]) >= rnorm(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
    }else if(reward == 'Bern'){
      if(rbeta(1, arm1_params[1], arm1_params[2]) >= rbeta(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
    } 
    #arm1draws = append(arm1draws,arm1_params[2])
    arm_path = append(arm_path, arm)
    #observe reward and compute ATE
    if(reward=='Norm'){
      y_i = rnorm(1,mean=mus[arm+1], sd=sigma)
    }else if(reward=='Bern'){
      y_i = rbinom(1,1, mus[arm+1])
    }else if(reward=='Cauchy'){ #TS updates using normal-normal but actual rewards drawn from cauchy
      y_i = rcauchy(1, mus[arm+1], scale=sigma)
    }else if(reward=='tdist'){ #TS updates using normal-normal but actual rewards drawn from cauchy
      y_i = rt(1, ncp=mus[arm+1], df=sigma)
    }
    if(arm==1){
      ate_i = y_i/p_1
      sigma2_hat_i = (y_i)**2/(p_1)**2
    }else{
      ate_i = -y_i/p_0
      sigma2_hat_i = (y_i)**2/(p_0)**2
    } 
    rewards = append(rewards, y_i)
    ates_i = append(ates_i, ate_i)
    ate_i = (1/i)*sum(ates_i)
    ate = append(ate, ate_i)
    #compute CS
    sig2_hats = append(sig2_hats, sigma2_hat_i)
    #update parameters
    if(reward %in% c('Norm', 'Cauchy', 'tdist')){
      if(arm==1){
        arm1_params[1] = (arm1_params[2]^2/(arm1_params[2]^2 + sigma^2))*y_i + (sigma^2/(arm1_params[2]^2 + sigma^2))*arm1_params[1]
        arm1_params[2] = 1/((1/sigma^2) + (1/arm1_params[2]^2))
      }else{
        arm0_params[1] = (arm0_params[2]^2/(arm0_params[2]^2 + sigma^2))*y_i + (sigma^2/(arm0_params[2]^2 + sigma^2))*arm0_params[1]
        arm0_params[2] = 1/((1/sigma^2) + (1/arm0_params[2]^2))
      }
    }else if (reward == 'Bern'){
      if(arm==1){
        arm1_params[1] = arm1_params[1]+y_i
        arm1_params[2] = arm1_params[2]+1-y_i 
      }else{
        arm0_params[1] = arm0_params[1]+ y_i
        arm0_params[2] = arm0_params[2]+1-y_i
      }
    }
  }
  return(list(  'ate' = ate,
                'arm_path' = arm_path,
                'sig_hats' = sig2_hats,
                'rewards' = rewards#,
                # 'arm1draws' = arm1draws
  ))
}

compute_CS_eta = function(eta, s_hats, alpha){
  V_is = (1/seq(1, t))*sqrt((((cumsum(s_hats) * eta**2)+1)/eta**2) * log(((cumsum(s_hats) * eta**2) + 1)/alpha**2))
  return(V_is)
}

# if alg==TS, armi_params is a list of prior parameters for reward distribution, if alg==UCB, armi_params is single float of the UCB for arm i

compute_TS_CSs_MAD = function(mus, alpha, t, arm1_params, arm0_params, delta, alg = 'TS', reward='Norm', alpha_ucb = 0.05){
  ates_i_MAD = c()
  ate_MAD = c()
  arm_path_MAD = c()
  sig2_hats_MAD = c()
  rewards_MAD = c()
  cum_rewards = c(0, 0)
  cum_pulls = c(0, 0)
  # V_is_MAD = c()
  for (i in seq(1:t)){
    #calculate/estimate P(arm 1 > arm 2)
    if(alg=='TS'){
      if (i == 1){
        p_1 = 0.5
      }else{
        if(reward %in% c('Norm', 'Cauchy', 'tdist')){
          p_1 = sum(rnorm(10000, arm1_params[1], arm1_params[2]) >= rnorm(10000, arm0_params[1], arm0_params[2]))/10000
        }else if(reward == 'Bern'){
          p_1 = sum(rbeta(10000, arm1_params[1], arm1_params[2]) >= rbeta(10000, arm0_params[1], arm0_params[2]))/10000
        } 
      }
    }else if(alg=='UCB'){
      #UCB
      if (i == 1){
        p_1 = 1
      }else if (i==2){
        p_1 = 0
      }else{
        p_1 = 1*(arm1_params >= arm0_params)
      }
    }
    p_0 = 1-p_1
    p_1_MAD = (delta[i])*0.5 + (1-delta[i])*p_1
    p_0_MAD = 1-p_1_MAD 
    #select arm
    #if we go with Bernoulli design
    if(rbinom(1,1,delta[i]) == 1){
      if(rbinom(1,1,0.5)==1){arm = 1}else{arm=0}
    }else{
      if(alg=='TS'){
        if(reward %in% c('Norm', 'Cauchy', 'tdist')){
          if(rnorm(1, arm1_params[1], arm1_params[2]) >= rnorm(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
        }else{
          if(rbeta(1, arm1_params[1], arm1_params[2]) >= rbeta(1, arm0_params[1], arm0_params[2])){arm = 1}else{arm=0}
        }
      }else if (alg=='UCB'){
        if(p_1==1){arm=1}else{arm=0}
      }
    }
    arm_path_MAD = append(arm_path_MAD, arm)
    
    if(reward=='Norm'){
      y_i = rnorm(1,mean=mus[arm+1], sd=sigma)
    }else if (reward =='Bern'){
      y_i = rbinom(1,1, mus[arm+1])
    }else if(reward=='Cauchy'){ #TS updates using normal-normal but actual rewards drawn from cauchy
      y_i = rcauchy(1, mus[arm+1], scale=sigma)
    }else if(reward=='tdist'){ #TS updates using normal-normal but actual rewards drawn from cauchy
      y_i = rt(1, ncp=mus[arm+1], df=sigma)
    }
    
    #observe reward and compute ATE
    if(arm==1){
      # y_i = rnorm(1,mean=mus[arm+1], sd=sigma)
      ate_i = y_i/p_1_MAD
      sigma2_hat_i = (y_i)**2/(p_1_MAD)**2
      cum_rewards[2] = cum_rewards[2] + y_i
      cum_pulls[2] = cum_pulls[2] + 1
    }else{
      #y_i = rnorm(1,mean=mus[arm+1], sd=sigma)
      ate_i = -y_i/p_0_MAD
      sigma2_hat_i = (y_i)**2/(p_0_MAD)**2
      cum_rewards[1] = cum_rewards[1] + y_i
      cum_pulls[1] = cum_pulls[1] + 1
    } 
    rewards_MAD = append(rewards_MAD, y_i)
    ates_i_MAD = append(ates_i_MAD, ate_i)
    ate_i = (1/i)*sum(ates_i_MAD)
    ate_MAD = append(ate_MAD, ate_i)
    #compute CS
    sig2_hats_MAD = append(sig2_hats_MAD, sigma2_hat_i)
    #update parameters
    if(alg == 'TS'){
      if(reward %in% c('Norm', 'Cauchy', 'tdist')){
        if(arm==1){
          arm1_params[1] = (arm1_params[2]^2/(arm1_params[2]^2 + sigma^2))*y_i + (sigma^2/(arm1_params[2]^2 + sigma^2))*arm1_params[1]
          arm1_params[2] = 1/((1/sigma^2) + (1/arm1_params[2]^2))
        }else{
          arm0_params[1] = (arm0_params[2]^2/(arm0_params[2]^2 + sigma^2))*y_i + (sigma^2/(arm0_params[2]^2 + sigma^2))*arm0_params[1]
          arm0_params[2] = 1/((1/sigma^2) + (1/arm0_params[2]^2))
        }
      } else if (reward=='Bern'){
        if(arm==1){
          arm1_params[1] = arm1_params[1]+y_i
          arm1_params[2] = arm1_params[2]+1-y_i 
        }else{
          arm0_params[1] = arm0_params[1]+ y_i
          arm0_params[2] = arm0_params[2]+1-y_i}
      }
    } else if (alg == 'UCB'){
      if(t>=3){
        #if UCB, arm1_params is just the UCB for arm 1 
        arm1_params = cum_rewards[2] + (alpha_ucb * log(i))/(2*cum_pulls[2])
        arm0_params = cum_rewards[1] + (alpha_ucb * log(i))/(2*cum_pulls[1])
      }
    }
  }
  
  return(list(  'ate' = ate_MAD,
                'arm_path' = arm_path_MAD,
                'sig_hats' = sig2_hats_MAD,
                'rewards' = rewards_MAD
  ))
}

# #compute_TS_Norm_CSs_standard(mus, alpha, t, arm1_params, arm0_params)
# x = compute_TS_CSs_standard(mus, alpha, t, arm1_params, arm0_params, reward='Bern')
# q = compute_TS_CSs_MAD(mus, alpha, t, arm1_params, arm0_params, delta, alg = 'TS', reward='Bern', alpha_ucb = 0.05)
# mean(x$arm_path) #essentially only drew from one arm! so ofc we sometimes get Type I error inflation, we can't actually estiamte the ATE if we only draw from one arm
# mean(q$arm_path)

if(alg == 'TS'){
  #compute CS's for both methods many times and show all CS paths
  testing = lapply(seq(1, N), function(i){
    print(i)
    standard_results = compute_TS_CSs_standard(mus, alpha, t, arm1_params, arm0_params, reward)
    ate_standard = standard_results[[1]]
    arm_path_standard = standard_results[[2]]
    sig2_hats_standard = standard_results[[3]]
    rewards_standard = standard_results[[4]]
    V_i_standard = compute_CS_eta(eta, sig2_hats_standard, alpha)
    return(list(
      'ate' = ate_standard,
      'arm_path' = arm_path_standard,
      'sig_hats' = sig2_hats_standard,
      'rewards' = rewards_standard,
      'CS' = V_i_standard
    ))
  })
}


#compute CS's for MAD for both methods many times and show all CS paths
testing_MAD = lapply(seq(1, N), function(i){
  standard_results = compute_TS_CSs_MAD(mus, alpha, t, arm1_params, arm0_params, delta, alg, reward, alpha_ucb = 0.05)
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
  standard_results = compute_TS_CSs_MAD(mus, alpha, t, arm1_params, arm0_params, delta = 1/seq(1, t)**0, alg = alg, reward=reward    , alpha_ucb = 0.05)
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
bern_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing_bernoulli_des[[i]]$CS}), byrow = FALSE, nrow = 1))
MAD_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing_MAD[[i]]$CS}), byrow = FALSE, nrow = 1))

if(alg == 'TS'){
  standard_ate = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$ate}), byrow = FALSE, nrow = 1))
  standard_Vt = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$CS}), byrow = FALSE, nrow = 1))
}
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
if(alg == 'TS'){
  CS_ate_df = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                               'ate' = bern_ate, 'Vt' = bern_Vt, 'replicate' = iteration),
                    data.frame('Method' = c('MAD'), 'Time' = time, 
                               'ate' = MAD_ate, 'Vt' = MAD_Vt,  'replicate' = iteration),
                    data.frame('Method' = c('Standard Design'), 'Time' = time, 
                               'ate' = standard_ate, 'Vt' = standard_Vt,  'replicate' = iteration))
}else{
  CS_ate_df = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                               'ate' = bern_ate, 'Vt' = bern_Vt, 'replicate' = iteration),
                    data.frame('Method' = c('MAD'), 'Time' = time, 
                               'ate' = MAD_ate, 'Vt' = MAD_Vt,  'replicate' = iteration))
}
ate = true_mu1 - true_mu0 
#TODO: SAVE DATA
write.csv(CS_ate_df, paste0('new_CS_ate_df_reward_', reward,'_df_',sigma,'_alg_', alg, '_ate_', ate, '_truemu1_',true_mu1, '_truemu0_', true_mu0,'.csv'))


#Generate upper and lower bounds on CSs
#For Bernoulli Rewards
#width and ATE can end up being -inf and inf because of how TS behaves, so need to adjust for that

CS_ate_df['wid'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){
  if( is.na(CS_ate_df$Vt[i]) | is.nan(CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$Vt[i])){1}else{CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
#CS_ate_df['width'] = sapply(seq(1, 3*N*t), function(i){if( is.na(CS_ate_df$Vt[i]) | is.nan( CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$Vt[i])){-1}else{ CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
CS_ate_df['ate_reg'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){
  if( is.na(CS_ate_df$ate[i]) | is.nan(CS_ate_df$ate[i]) | is.infinite(CS_ate_df$ate[i])){0}else{CS_ate_df$ate[i]}})




if(reward=='Bern'){
  #width plotting for Bernoulli
  CS_ate_df$plottingwidth = sapply(CS_ate_df$wid, function(i){if(i>1){1}else{i}})
  
  plotting_df = CS_ate_df %>% group_by(Method, Time) %>%
    summarize(ate_se = sd(ate_reg, na.rm = TRUE)/sqrt(N), mean_ate = mean(ate_reg, na.rm=TRUE), width = mean(wid, na.rm=TRUE), width_se = sd(wid, na.rm=TRUE)/sqrt(N), plotting_width = mean(plottingwidth, na.rm=TRUE), plotting_width_se = sd(plottingwidth, na.rm=TRUE)/sqrt(N))
  
  width_plot = ggplot(plotting_df, aes(x = Time, y = plotting_width, fill=Method, color=Method)) + geom_line() + 
    #geom_errorbar(aes(ymin=width-2*width_se, ymax = width+2*width_se)) + 
    geom_ribbon(aes(ymin=plotting_width-2*plotting_width_se, ymax = plotting_width+2*plotting_width_se), alpha=0.4, show.legend=FALSE, colour=NA)+ #ylim(c(-100, 100))+
    # scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
    #scale_x_continuous(expand = c(0, 0))+
    theme_minimal() + 
    xlab('Time (log10 Scale)')+scale_x_log10(expand = c(0, 0)) +  ylab('Average Width') + ylim(c(0, 1))+
    scale_linetype_manual(name = "", values = c(1)) + #+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
    theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
          panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 
  width_plot
}else{
  plotting_df = CS_ate_df %>% group_by(Method, Time) %>%
    summarize(ate_se = sd(ate_reg, na.rm = TRUE)/sqrt(N), mean_ate = mean(ate_reg, na.rm=TRUE), width = mean(wid, na.rm=TRUE), width_se = sd(wid, na.rm=TRUE)/sqrt(N) )
  
  width_plot = ggplot(plotting_df, aes(x = Time, y = width, fill=Method, color=Method)) + geom_line() + 
    #geom_errorbar(aes(ymin=width-2*width_se, ymax = width+2*width_se)) + 
    geom_ribbon(aes(ymin=width-2*width_se, ymax = width+2*width_se), alpha=0.4, show.legend=FALSE, colour=NA)+ #ylim(c(-100, 100))+
    # scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))+
    #scale_x_continuous(expand = c(0, 0))+
    theme_minimal() + 
    xlab('Time (log10 Scale)')+scale_x_log10(expand = c(0, 0)) +  ylab('Average Width') +
    scale_linetype_manual(name = "", values = c(1)) + #+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
    theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
          panel.border = element_blank(), legend.title = element_text(size = 22),  strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) 
  width_plot
}


#Generate upper and lower bounds on CSs
#For Bernoulli Rewards
if(reward == 'Bern'){
  CS_ate_df['upper'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){
    if(CS_ate_df$ate[i] + CS_ate_df$Vt[i] > 1 | is.na(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] + CS_ate_df$Vt[i])){1}else{CS_ate_df$ate[i] + CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
  CS_ate_df['lower'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){if(CS_ate_df$ate[i] - CS_ate_df$Vt[i] < -1 | is.na(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] - CS_ate_df$Vt[i])){-1}else{CS_ate_df$ate[i] - CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
}else{
  # For Normal Rewards
  CS_ate_df['upper'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){
    if( is.na(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] + CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] + CS_ate_df$Vt[i])){1e3}else{CS_ate_df$ate[i] + CS_ate_df$Vt[i]}})
  CS_ate_df['lower'] = sapply(seq(1, length(CS_ate_df$Vt)), function(i){if(is.na(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] - CS_ate_df$Vt[i])){1e3}else{CS_ate_df$ate[i] - CS_ate_df$Vt[i]}})##CS_ate_df$ate -CS_ate_df$Vt#sapply(seq(1, 3*N*t), function(i){if(CS_ate_df$ate[i] - CS_ate_df$Vt[i] < -1 | is.na(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.nan(CS_ate_df$ate[i] - CS_ate_df$Vt[i]) | is.infinite(CS_ate_df$ate[i] - CS_ate_df$Vt[i])){-1}else{CS_ate_df$ate[i] - CS_ate_df$Vt[i]}})# CS_ate_df$ate + CS_ate_df$Vt, CS_ate_df$ate + CS_ate_df$Vt >= 1)
  
}


#Stopping Time: proportion of replicates where 0 is within CS
CS_ate_df['stopped'] = sapply(seq(1, length(CS_ate_df$ate)), function(i){if(between(0, CS_ate_df$lower[i], CS_ate_df$upper[i])){1}else{0}})
stopped_df = CS_ate_df%>% group_by(Method, Time) %>% summarize(avg_stop = mean(stopped), sd = sd(stopped)/sqrt(N))

plotting_df['avg_stop'] = stopped_df$avg_stop
plotting_df['stopped_se'] = stopped_df$sd
#can get long runs of just one arm that skews your ATE estimation drastically towards the mean of that arm

stopped = ggplot(plotting_df, aes(x = Time, y = avg_stop, color=Method, fill=Method)) + geom_line() +
  #  geom_errorbar(aes(ymin=avg_stop-2*sd, ymax = avg_stop+2*sd)) + 
  geom_ribbon(aes(ymin=avg_stop-2*stopped_se, ymax = avg_stop+2*stopped_se), alpha=0.4, show.legend=FALSE, colour=NA)+
  theme_minimal() + ylab('Average Proportion Stopped')+
  #geom_hline(aes(yintercept=alpha), linetype = 'dashed', alpha = 0.5) +
  scale_x_log10(expand=c(0, 0)) + 
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))+ xlab('Time (log10 Scale)') #+ 
#scale_linetype_manual(name = expression(paste(1-alpha)), values = c(2))#+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
stopped

#coverage over time: proportion of replicates where the true treatment effect is within the CS
CS_ate_df['coverage'] = sapply(seq(1, length(CS_ate_df$ate)), function(i){if(between(mus[2]-mus[1], CS_ate_df$lower[i], CS_ate_df$upper[i])){1}else{0}})
#check proportion of replicates such that all true ATE are inside the CS up until t
cov_df = CS_ate_df%>% group_by(Method, replicate) %>% mutate(emp_cov = cumsum(coverage))
CS_ate_df['emp_coverage'] = 1*(cov_df$emp_cov == cov_df$Time)

cov_df_1 = CS_ate_df%>% group_by(Method, Time) %>% summarize(avg_cov = mean(emp_coverage), sd = sd(emp_coverage)/sqrt(N))


plotting_df['avg_cov'] = cov_df_1$avg_cov
plotting_df['cov_se'] = cov_df_1$sd

coverage = ggplot(plotting_df, aes(x = Time, y = avg_cov, color=Method, fill=Method)) + geom_line() +
  geom_ribbon(aes(ymax=pmin(avg_cov + 2*cov_se, 1), ymin = avg_cov-2*cov_se), alpha=0.4, show.legend=FALSE, colour=NA)+ theme_minimal() + ylab('Coverage')+
  geom_hline(aes(yintercept=1-alpha), linetype = 'dashed', alpha = 0.5) +
  scale_x_log10(expand=c(0, 0)) + 
  theme(axis.title =element_text(size=25), axis.text = element_text(size = 20),legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))+ xlab('Time (log10 Scale)') #+ 
#scale_linetype_manual(name = expression(paste(1-alpha)), values = c(2))#+   geom_hline(aes(yintercept=1)) +   geom_hline(aes(yintercept=0), alpha=0.5)
coverage



#generate dataframe of rewards from each method
#reward over the 400 timesteps, across the 100 iterations
bern_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing_bernoulli_des[[i]]$reward}), byrow = FALSE, nrow = 1))
MAD_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing_MAD[[i]]$reward}), byrow = FALSE, nrow = 1))
if(alg=='TS'){
  standard_rewards = as.vector(matrix(sapply(seq(1, N), function(i){testing[[i]]$reward}), byrow = FALSE, nrow = 1))
}
time = rep(seq(1, t, 1), N)
iteration = rep(seq(1, N), each = t)
if(alg=='TS'){
  rewards = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                             'reward' = bern_rewards, 'replicate' = iteration),
                  data.frame('Method' = c('MAD'), 'Time' = time, 
                             'reward' = MAD_rewards, 'replicate' = iteration),
                  data.frame('Method' = c('Standard Design'), 'Time' = time, 
                             'reward' = standard_rewards, 'replicate' = iteration))
}else{
  rewards = rbind(data.frame('Method' = c('Bernoulli Design'), 'Time' = time, 
                             'reward' = bern_rewards, 'replicate' = iteration),
                  data.frame('Method' = c('MAD'), 'Time' = time, 
                             'reward' = MAD_rewards, 'replicate' = iteration))
}



rewards_plotting = rewards %>% group_by(Method, Time) %>% summarize(avg_reward = mean(reward), se_avg_reward = sd(reward)/sqrt(N))%>% mutate(cum_avg_reward = cumsum(avg_reward), time_avg_reward = (1/Time)*cumsum(avg_reward))
#rewards %>% group_by(Method, Time) %>% summarize(se_reward = sd(reward))

rewards_ses = rewards %>% group_by(Method, replicate) %>% mutate(cumsum_reward_per_it = cumsum(reward),
                                                                 time_avg_reward_per_it = cumsum(reward)/Time)%>% 
  group_by(Method, Time) %>% summarize(se_cum_reward = sd(cumsum_reward_per_it), se_time_avg_reward = sd(time_avg_reward_per_it))

rewards_plotting['se_avg_reward'] = rewards_ses$se_cum_reward/sqrt(N)
rewards_plotting['se_cum_reward'] = rewards_ses$se_cum_reward/sqrt(N)
rewards_plotting['se_time_avg_reward'] = rewards_ses$se_time_avg_reward/sqrt(N)



plotting_df = cbind(plotting_df, rewards_plotting[,seq(3, length(names(rewards_plotting)))])


write.csv(plotting_df, paste0('~/Documents/Research/MAD/data/plotting_df_', reward,'_df_',sigma,'_alg_', alg, '_ate_', ate, '_truemu1_',true_mu1, '_truemu0_', true_mu0,'.csv'))


