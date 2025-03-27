
setwd("~/TD/5_real")
n.cluster = 120

# setwd("~/Dropbox/Research/AbbVie/Time_drfit/Final_v3/TD_5/")
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

# data.obs = data.frame("time" = c(0, 1, 12, 30),
#                       "AN" = c(0.46, 0.23, 0.22, 0.32))
# data.obs$time2 = (data.obs$time)^2
# data.obs$time22 = log(35-data.obs$time)
# data.obs$time3 = (data.obs$time)^1.3
# data.obs$logtime = log(data.obs$time+10)
# data.pred = data.frame("time" = seq(0, 30, length.out = 1000))
# data.pred$time2 = (data.pred$time)^2
# data.pred$time22 = log(35-data.pred$time)
# data.pred$time3 = (data.pred$time)^1.3
# data.pred$logtime = log(data.pred$time+10)
# 
# ## quadratic
# fit.1 = lm(AN~time + time2, data = data.obs)
# 
# print(summary(fit.1))
# plot(data.pred$time, predict(fit.1, newdata = data.pred))
# points(data.obs$time, data.obs$AN)
# 
# ## cubic
# fit.2 = lm(AN~time + time2+time3, data = data.obs)
# 
# print(summary(fit.2))
# plot(data.pred$time, predict(fit.2, newdata = data.pred))
# print(summary(predict(fit.2, newdata = data.pred)))
# points(data.obs$time, data.obs$AN)
# 
# ## order of 4
# fit.3 = lm(AN~ time + logtime + time2, data = data.obs)
# 
# print(summary(fit.3))
# plot(data.pred$time, predict(fit.3, newdata = data.pred))
# points(data.obs$time, data.obs$AN)



# ## piecewise linear
# fit.3 = function(x){
#   (x<=1.30)*(0.46-0.23*x)+(x>=1.30)*(0.1/18*x+0.22-12*0.1/18)
# }
# 
# plot(data.pred$time, fit.3(data.pred$time))
# print(round(fit.3(data.obs$time), 3))

# mu.time.func.1 = function(time){
#   0.36-0.021*time+0.00064*time^2
# }
#
# mu.time.func.2 = function(time){
#   0.46-0.257*time+0.0273*time^2-0.00063*time^3
# }
#
# mu.time.func.3 = function(time){
#   26.57+0.863*time-11.34*log(time+10)-0.0114*time^2
# }
#
# plot(seq(0, 30, length.out = 1000), mu.time.func.2(seq(0, 30, length.out = 1000)))
# points(seq(0, 30, length.out = 1000), mu.time.func.3(seq(0, 30, length.out = 1000)))
# points(seq(0, 30, length.out = 1000), mu.time.func.1(seq(0, 30, length.out = 1000)))


#################################

one.alpha = 0.025
mu1 = 0
trt.prop = 2/4
trt.diff.mag = 0.12
n = 400
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

n.scen = 12
table.out = matrix(NA, nrow = n.scen, ncol = 16)
table.name.out = rep(NA, n.scen)
time.start = Sys.time()

for (scen.ind in c(1:12)){
  print(scen.ind)
  if (scen.ind==1){
    trt.diff = 0; sd.1 = 0.3; sd.2 = 0.3
    mu.time.func = function(time){
      0.36-0.021*time+0.00065*time^2
    }
  } else if (scen.ind==2){
    trt.diff = trt.diff.mag; sd.1 = 0.3; sd.2 = 0.3
    mu.time.func = function(time){
      0.36-0.021*time+0.00065*time^2
    }
  } else if (scen.ind==3){
    trt.diff = 0; sd.1 = 0.3; sd.2 = 0.3
    
    mu.time.func = function(time){
      0.46-0.507*time-0.00977*time^2+0.287*time^1.3
    }
    
  }   else if (scen.ind==4){
    trt.diff = trt.diff.mag; sd.1 = 0.3; sd.2 = 0.3
    
    mu.time.func = function(time){
      0.46-0.507*time-0.00977*time^2+0.287*time^1.3
    }
    
  } else if (scen.ind==5){
    trt.diff = 0; sd.1 = 0.3; sd.2 = 0.3

    mu.time.func = function(time){
      26.57+0.863*time-11.34*log(time+10)-0.0114*time^2
    }
    
  } else if (scen.ind==6){
    trt.diff = trt.diff.mag; sd.1 = 0.3; sd.2 = 0.3

    mu.time.func = function(time){
      26.57+0.863*time-11.34*log(time+10)-0.0114*time^2
    }
    
  } else if (scen.ind==7){
    trt.diff = 0; sd.1 = 0.4; sd.2 = 0.2
    mu.time.func = function(time){
      0.36-0.021*time+0.00065*time^2
    }
  } else if (scen.ind==8){
    trt.diff = trt.diff.mag; sd.1 = 0.4; sd.2 = 0.2
    mu.time.func = function(time){
      0.36-0.021*time+0.00065*time^2
    }
  } else if (scen.ind==9){
    trt.diff = 0; sd.1 = 0.4; sd.2 = 0.2
    
    mu.time.func = function(time){
      0.46-0.507*time-0.00977*time^2+0.287*time^1.3
    }
    
  }   else if (scen.ind==10){
    trt.diff = trt.diff.mag; sd.1 = 0.4; sd.2 = 0.2
    
    mu.time.func = function(time){
      0.46-0.507*time-0.00977*time^2+0.287*time^1.3
    }
    
  } else if (scen.ind==11){
    trt.diff = 0; sd.1 = sd.1 = 0.4; sd.2 = 0.2
    
    mu.time.func = function(time){
      26.57+0.863*time-11.34*log(time+10)-0.0114*time^2
    }
    
  } else if (scen.ind==12){
    trt.diff = trt.diff.mag; sd.1 = 0.4; sd.2 = 0.2
    
    mu.time.func = function(time){
      26.57+0.863*time-11.34*log(time+10)-0.0114*time^2
    }
    
  } 
  
  
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  p.vec = foreach(itt = 1:n.itt) %dopar% {
    
    set.seed(itt + scen.ind*n.itt*2)
    
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
                           "time" = (0:(n-1))/(n-1)*30, 
                           "block.ind" = rep(1:n.block.rand, each = n.block.size),
                           "grp" = 1)
    
    data.temp$grp[sort(sample(1:n, trt.prop*n))] = 2
  
    data.temp$mu_base = mu1
    data.temp$mu_base[data.temp$grp==2] = mu1+trt.diff
    
    data.temp$mu_time = mu.time.func(data.temp$time)
    
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
    
    
    ## Method 5: randomization test
    rt.test.func = function(data.in, bs.knots.in, bs.degree.in){
      per.fit = sapply(1:n.B, function(itt.B){
        grp.vec.temp = rep(1, n)
        grp.vec.temp[sort(sample(1:n, n*(trt.prop)))] = 2
        
        # grp.vec.temp = rep(1, n)
        # for (i in 1:n.block.rand){
        #   grp.vec.temp[sort(sample(1:(n.block.size), n.block.size*trt.prop))+
        #                   (i-1)*n.block.size] = 2
        # }
        
        # data.rand = data.in
        # data.rand$grp = grp.vec.temp 
        # rlm.fit.init = lm(x~grp+bs(ind, knots = bs.knots.in, 
        #                            degree = bs.degree.in), 
        #                   data = data.rand)
        # data.rand$rlm_weight = (c(1/mean((rlm.fit.init$residuals[data.rand$grp==1])^2),
        #               1/mean((rlm.fit.init$residuals[data.rand$grp==2])^2)))
        # 
        # lm.2.fit = lm(x~grp+bs(ind, knots = bs.knots.in, degree = bs.degree.in), 
        #               data = data.rand, 
        #               weights = rlm_weight)
        
        return(c(-t.test(data.in$x[grp.vec.temp==1],
                         data.in$x[grp.vec.temp==2],
                         alternative = "less")$statistic
               # summary(lm.2.fit)$coefficients[2, 3]
               ))
        
      })
      
      t.obs.fit = t.test(data.in$x[data.in$grp==1],
                         data.in$x[data.in$grp==2],
                               alternative = "less")
      
      # rlm.fit.init = lm(x~grp+bs(ind, knots = bs.knots.in, 
      #                            degree = bs.degree.in), 
      #                   data = data.in)
      # data.in$rlm_weight = (c(1/mean((rlm.fit.init$residuals[data.in$grp==1])^2),
      #                           1/mean((rlm.fit.init$residuals[data.in$grp==2])^2)))
      # 
      # lm.obs.fit = lm(x~grp+bs(ind, knots = bs.knots.in, degree = bs.degree.in), 
      #               data = data.in, 
      #               weights = rlm_weight)
      
      return(c(mean((-t.obs.fit$statistic)<=per.fit)
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

    ## output
    return(c(p.t.out[1]<=one.alpha, 
             p.wlm.out[1]<=one.alpha, 
             p.rt.out<=one.alpha, 
             p.cust.out.1[1]<=one.alpha,
             
             ## t-test
             p.t.out[2] - trt.diff,
             (p.t.out[2]+qnorm(1-one.alpha)*p.t.out[3]>=trt.diff)&
               (p.t.out[2]-qnorm(1-one.alpha)*p.t.out[3]<=trt.diff),
             
             ## wlm
             p.wlm.out[2] - trt.diff,
             (p.wlm.out[2]+qnorm(1-one.alpha)*p.wlm.out[3]>=trt.diff)&
               (p.wlm.out[2]-qnorm(1-one.alpha)*p.wlm.out[3]<=trt.diff),
             
             ## cust
             p.cust.out.1[2] - trt.diff,
             (p.cust.out.1[2]+qnorm(1-one.alpha)*p.cust.out.1[3]>=trt.diff)&
               (p.cust.out.1[2]-qnorm(1-one.alpha)*p.cust.out.1[3]<=trt.diff)
             
             ))
  }
  
  stopCluster(cl)
  p.vec.out = matrix(unlist(p.vec),nrow = n.itt, ncol = 10, byrow = TRUE)
  
  # print(mean((p.vec.out[, 8])))
  # print(sd((p.vec.out[, 8])))
  # print(hist(p.vec.out[, 6]))
  
  table.name.out[scen.ind] = paste0(deparse(mu.time.func), collapse = "")
  table.out[scen.ind,] =   c(trt.diff, sd.1, sd.2, 
              apply(p.vec.out, 2, function(x){mean(x)}),
              apply(p.vec.out[, c(5, 7, 9)], 2, function(x){sd(x)}))

  print(table.out)
}

print(Sys.time()-time.start)

colnames(table.out) = c("trt_diff", "sd1", "sd2", "t_test", "wlm", 
                        "rt_t",  "cust.1",
                        "t_bias", "t_cp",
                        "wlm_bias",  "wlm_cp",
                        "cust_bias",  "cust_cp", 
                        "t_sd", "wlm_sd","cust_sd")

write.csv(cbind(table.name.out, table.out), 
          paste0("Real_time_drift_trtprop_", trt.prop, 
                 ".csv", sep=""))















