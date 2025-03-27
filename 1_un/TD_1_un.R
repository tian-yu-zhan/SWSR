
setwd("~/TD/1_un")
n.cluster = 120

# setwd("~/Dropbox/Research/AbbVie/Time_drfit/Final_v3/TD_1/")
# n.cluster = 8

# library(flipr)
library(segmented)
library(doParallel)
library(MASS)
library(msme)
library(robustreg)
library(rpart)
library(rpart.plot)
library(glmnet)
library(dplyr)
library(splines)

one.alpha = 0.025
mu1 = 0
sd.1 = 0.4
sd.2 = 0.2
trt.prop = 3/4
trt.diff.mag = 0.13
n = 600
n.block.size = 40
n.block.rand = n/n.block.size 
n.cross = 5
n.block = n/n.cross
n.itt = 10^5
n.B = 10^3

bs.knots.cand = c(1, 1, 5, 5)
bs.degree.cand = c(1, 2, 2, 3)

# bs.knots.cand = c(40)
# bs.degree.cand = c(5)

n.scen = 8
table.out = matrix(NA, nrow = n.scen, ncol = 27)
table.name.out = rep(NA, n.scen)
time.start = Sys.time()

for (scen.ind in c(1:8)){
  print(scen.ind)
  if (scen.ind==1){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0
    }
  } else if (scen.ind==2){
    trt.diff = trt.diff.mag
    mu.time.func = function(time, seed){
      0
    }
  } else if (scen.ind==3){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0.3/(n-1)*(time-1)
    }
  } else if (scen.ind==4){
    trt.diff = trt.diff.mag
    mu.time.func = function(time, seed){
      0.3/(n-1)*(time-1)
    }
  } else if (scen.ind==5){
    trt.diff = 0

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.002))
      time.temp = (c(0, cumsum(time.temp.1)))
      time.temp[time]
    }
    
  }   else if (scen.ind==6){
    trt.diff = trt.diff.mag

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.002))
      time.temp = (c(0, cumsum(time.temp.1)))
      time.temp[time]
    }
    
  } else if (scen.ind==7){
    trt.diff = 0

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.004))
      time.temp = (c(0, cumsum(time.temp.1)))
      time.temp[time]
    }
    
  }   else if (scen.ind==8){
    trt.diff = trt.diff.mag

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.004))
      time.temp = (c(0, cumsum(time.temp.1)))
      time.temp[time]
    }
    
  } 
  
  
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  p.vec = foreach(itt = 1:n.itt) %dopar% {
    
    set.seed(itt + scen.ind*n.itt)
    
    # library(flipr)
    library(segmented)
    library(MASS)
    library(msme)
    library(robustreg)
    library(rpart)
    library(glmnet)
    library(dplyr)
    library(splines)
    
    data.temp = data.frame("ind" = 1:n,
                           "grp" = 1)
    
    data.temp$grp[sort(sample(1:n, trt.prop*n))] = 2
  
    data.temp$mu_base = mu1
    data.temp$mu_base[data.temp$grp==2] = mu1+trt.diff
    
    data.temp$mu_time = mu.time.func(data.temp$ind, itt + scen.ind*n.itt)
    
    data.temp$mu = data.temp$mu_base+data.temp$mu_time
    data.temp$sd = sd.1
    data.temp$sd[data.temp$grp==2] = sd.2
    
    data.temp$x = rnorm(n, data.temp$mu, data.temp$sd)
    
    ## Method 1.1: Welch t test
    t.test.func = function(data.1, data.2){
      t.fit = t.test(data.1,
                     data.2,
                     alternative = "less")
      return(as.numeric(c(t.fit$p.value,
               t.fit$estimate[2]-t.fit$estimate[1],
               t.fit$stderr)))
    }
    
    p.t.out = t.test.func(data.1 = data.temp$x[data.temp$grp==1],
                          data.2 = data.temp$x[data.temp$grp==2])
    
    ## Method 1.2: Wilcox t test
    wil.test.func = function(data.1, data.2){
      wil.fit = wilcox.test(data.2,
                     data.1,
                     alternative = "greater", conf.int = TRUE)
      wil.two.fit = wilcox.test(data.2,
                            data.1,conf.int = TRUE)
      return(as.numeric(c(wil.fit$p.value,
               wil.fit$estimate,
               wil.two.fit$conf.int)))
    }
    
    p.wil.out = wil.test.func(data.1 = data.temp$x[data.temp$grp==1],
                          data.2 = data.temp$x[data.temp$grp==2])

    ## Method 2: Simple linear regression 
    lm.test.func = function(data.in){
      lm.fit = lm(x~grp+ind, data = data.in)

      return(as.numeric(c(1-pnorm(summary(lm.fit)$coefficients[2, 3]),
               summary(lm.fit)$coefficients[2, 1:2]
               )))

    }
    
    p.lm.out = lm.test.func(data.in = data.temp)
    
    ## Method 3:  weighted linear regression 
    wlm.test.func = function(data.in){
      lm.fit.init = lm(x~grp+ind, data = data.in)
      # data.in$lm_weight = (c(1/mean((lm.fit.init$residuals[data.in$grp==1])^2),
      #                          1/mean((lm.fit.init$residuals[data.in$grp==2])^2)))
      
      data.in$lm_weight[data.in$grp==1] = 
        1/mean((lm.fit.init$residuals[data.in$grp==1])^2)
      data.in$lm_weight[data.in$grp==2] = 
        1/mean((lm.fit.init$residuals[data.in$grp==2])^2)
      
      wlm.fit = lm(x~grp+ind, data = data.in, weights = lm_weight)
      
      return(as.numeric(c(1-pnorm(summary(wlm.fit)$coefficients[2, 3]),
                          summary(wlm.fit)$coefficients[2, 1:2])))
    }
    
    p.wlm.out = wlm.test.func(data.in = data.temp)
    
    ## Method 4:  iterative weighted least square 
    rlm.test.func = function(data.in){
      
      rlm.fit = rlm(x~grp+ind, data = data.in)
      
      return(as.numeric(c(1-pnorm(summary(rlm.fit)$coefficients[2, 3]),
                          summary(rlm.fit)$coefficients[2, 1:2])))
    }
    
    p.rlm.out = rlm.test.func(data.in = data.temp)
    
    ## Method 5: randomization test
    rt.test.func = function(data.in, bs.knots.in, bs.degree.in){
      
      per.fit = sapply(1:n.B, function(itt.B){
        grp.vec.temp = rep(1, n)
        grp.vec.temp[sort(sample(1:n, n*(trt.prop)))] = 2
        
        return(c(-t.test(data.in$x[grp.vec.temp==1],
                         data.in$x[grp.vec.temp==2],
                         alternative = "less")$statistic,
                 -wilcox.test(data.in$x[grp.vec.temp==1],
                              data.in$x[grp.vec.temp==2],
                              alternative = "less")$statistic
               # summary(lm.2.fit)$coefficients[2, 3]
               ))
        
      })
      
      t.obs.fit = t.test(data.in$x[data.in$grp==1],
                         data.in$x[data.in$grp==2],
                               alternative = "less")
      wil.obs.fit = wilcox.test(data.in$x[data.in$grp==1],
                         data.in$x[data.in$grp==2],
                         alternative = "less")
      
      return(c(mean((-t.obs.fit$statistic)<=per.fit[1, ]),
               mean((-wil.obs.fit$statistic)<=per.fit[2, ])
               # mean((summary(lm.obs.fit)$coefficients[2, 3])<=per.fit[3, ])
               ))
    }
    
    p.rt.out = rt.test.func(data.in = data.temp,
                            bs.knots.in = bs.knots,
                            bs.degree.in = bs.degree)
    
    ## Method 6: rlm + spline + bootstrap
    cust.test.func = function(data.in, method.ind){
      
    if (method.ind==1){
        out.cross = matrix(NA, nrow = length(bs.knots.cand), ncol = 5)
        for (cross.feature.itt in 1:length(bs.knots.cand)){
          bs.knots.temp.in = bs.knots.cand[cross.feature.itt]
          bs.degree.temp.in = bs.degree.cand[cross.feature.itt]
          
          for (cross.itt in 1:n.cross){
            
            cross.ind = sample(rep(1:n.cross, n/n.cross))
            
            data.cross.all = data.in
            
            rlm.fit.init = lm(x~grp+bs(ind, knots = 
                                         as.numeric(round(quantile(1:n, prob = (1:bs.knots.temp.in)/(bs.knots.temp.in+1)))),
                                       degree = bs.degree.temp.in), 
                              data = data.cross.all)
            # data.cross.all$rlm_weight = 
            #   (c(1/mean((rlm.fit.init$residuals[data.cross.all$grp==1])^2),
            #      1/mean((rlm.fit.init$residuals[data.cross.all$grp==2])^2)))
            
            data.cross.all$rlm_weight[data.cross.all$grp==1] =
              1/mean((rlm.fit.init$residuals[data.cross.all$grp==1])^2)
            data.cross.all$rlm_weight[data.cross.all$grp==2] =
              1/mean((rlm.fit.init$residuals[data.cross.all$grp==2])^2)
            
            data.cross.train = data.cross.all[!cross.ind==cross.itt, ]
            data.cross.val = data.cross.all[cross.ind==cross.itt, ]
            
            lm.2.fit = lm(x~grp+bs(ind, knots = 
                                     as.numeric(round(quantile(1:n, prob = (1:bs.knots.temp.in)/(bs.knots.temp.in+1)))),                         
                                   degree = bs.degree.temp.in), 
                          data = data.cross.train, 
                          weights = rlm_weight)
            
            lm.2.pred = predict(lm.2.fit, data.cross.val, weights = rlm_weight)
            out.cross[cross.feature.itt, cross.itt] = 
              mean((lm.2.pred-data.cross.val$x)^2)
            
          }
        }
        
        ## post cross-validation
        cross.feature.final = which.min(apply(out.cross, 1 , mean))
        
        bs.knots.final = bs.knots.cand[cross.feature.final]
        bs.degree.final = bs.degree.cand[cross.feature.final]
      } else if (method.ind==2){
        bs.knots.final = 1
        bs.degree.final = 1
      } else if (method.ind==3){
        bs.knots.final = 100
        bs.degree.final = 3
      }
      
      rlm.fit.init = lm(x~grp+bs(ind, knots = 
                                   as.numeric(round(quantile(1:n, prob = (1:bs.knots.final)/(bs.knots.final+1)))),                                  
                                 degree = bs.degree.final), 
                        data = data.in)
      # data.in$rlm_weight = (c(1/mean((rlm.fit.init$residuals[data.in$grp==1])^2),
      #                         1/mean((rlm.fit.init$residuals[data.in$grp==2])^2)))
      
      data.in$rlm_weight[data.in$grp==1] = 
        1/mean((rlm.fit.init$residuals[data.in$grp==1])^2)
      data.in$rlm_weight[data.in$grp==2] = 
        1/mean((rlm.fit.init$residuals[data.in$grp==2])^2)
      
      lm.2.fit = lm(x~grp+bs(ind, knots = 
                               as.numeric(round(quantile(1:n, prob = (1:bs.knots.final)/(bs.knots.final+1)))),  
                             degree = bs.degree.final), 
                    data = data.in, 
                    weights = rlm_weight)
      
      return(as.numeric(c(1-pnorm(summary(lm.2.fit)$coefficients[2, 3]),
                          summary(lm.2.fit)$coefficients[2, 1:2])))
  
      
    
      
    }
    
    
    
    p.cust.out.1 = cust.test.func(data.in = data.temp, method.ind=1)
    # p.cust.out.2 = cust.test.func(data.in = data.temp, method.ind=2)
    # p.cust.out.3 = cust.test.func(data.in = data.temp, method.ind=3)

    ## output
    return(c(p.t.out[1]<=one.alpha, 
             p.wil.out[1]<=one.alpha, 
             p.lm.out[1]<=one.alpha, 
             p.wlm.out[1]<=one.alpha, 
             p.rlm.out[1]<=one.alpha, 
             p.rt.out<=one.alpha, 
             p.cust.out.1[1]<=one.alpha,
             
             ## t-test
             p.t.out[2] - trt.diff,
             (p.t.out[2]+qnorm(1-one.alpha)*p.t.out[3]>=trt.diff)&
               (p.t.out[2]-qnorm(1-one.alpha)*p.t.out[3]<=trt.diff),
             
             ## wil-test
             p.wil.out[2] - trt.diff,
             (p.wil.out[3]<=trt.diff)&
               (p.wil.out[4]>=trt.diff),
             
             ## lm
             p.lm.out[2] - trt.diff,
             (p.lm.out[2]+qnorm(1-one.alpha)*p.lm.out[3]>=trt.diff)&
               (p.lm.out[2]-qnorm(1-one.alpha)*p.lm.out[3]<=trt.diff),
             
             ## wlm
             p.wlm.out[2] - trt.diff,
             (p.wlm.out[2]+qnorm(1-one.alpha)*p.wlm.out[3]>=trt.diff)&
               (p.wlm.out[2]-qnorm(1-one.alpha)*p.wlm.out[3]<=trt.diff),
             
             ## rlm
             p.rlm.out[2] - trt.diff,
             (p.rlm.out[2]+qnorm(1-one.alpha)*p.rlm.out[3]>=trt.diff)&
               (p.rlm.out[2]-qnorm(1-one.alpha)*p.rlm.out[3]<=trt.diff),
             
             ## cust
             p.cust.out.1[2] - trt.diff,
             (p.cust.out.1[2]+qnorm(1-one.alpha)*p.cust.out.1[3]>=trt.diff)&
               (p.cust.out.1[2]-qnorm(1-one.alpha)*p.cust.out.1[3]<=trt.diff)
             
             ))
  }
  stopCluster(cl)
  
  p.vec.out = matrix(unlist(p.vec),nrow = n.itt, ncol = 20, byrow = TRUE)
  
  # print(mean((p.vec.out[, 8])))
  # print(sd((p.vec.out[, 8])))
  # print(hist(p.vec.out[, 6]))
  
  table.name.out[scen.ind] = paste0(deparse(mu.time.func), collapse = "")
  table.out[scen.ind,] =   c(trt.diff,
              apply(p.vec.out, 2, function(x){mean(x)}),
              apply(p.vec.out[, c(9, 11, 13, 15, 17, 19)], 2, function(x){sd(x)}))

  print(table.out)
}

print(Sys.time()-time.start)

colnames(table.out) = c("trt_diff", "t_test", "wil_test", 
                        "lm", "wlm", "rlm",
                        "rt_t", "rt_wil", "cust.1",
                        "t_bias", "t_cp",
                        "wil_bias", "wil_cp",
                        "lm_bias",  "lm_cp",
                        "wlm_bias",  "wlm_cp",
                        "rlm_bias",  "rlm_cp",
                        "cust_bias",  "cust_cp", 
                        "t_sd", "wil_sd", "lm_sd", "wlm_sd",
                        "rlm_sd", "cust_sd")

write.csv(cbind(table.name.out, table.out), 
          paste0("Time_drift_sd1_", sd.1, "_sd2_", sd.2, 
                 "_trtprop_", trt.prop, 
                 ".csv", sep=""))















