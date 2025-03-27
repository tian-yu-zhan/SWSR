
setwd("~/TD/6_error")
n.cluster = 120

# setwd("~/Dropbox/Research/AbbVie/Time_drfit/Final_v3/TD_6/")
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
n.itt = 10^6
n.B = 10^3

bs.knots.cand = c(1, 1, 5, 5)
bs.degree.cand = c(1, 2, 2, 3)

# bs.knots.cand = c(40)
# bs.degree.cand = c(5)

n.scen = 6
table.out = matrix(NA, nrow = n.scen, ncol = 4)
table.name.out = rep(NA, n.scen)
time.start = Sys.time()

for (scen.ind in c(1:6)){
  print(scen.ind)
  if (scen.ind==1){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0.3 - 0.3/(n-1)*(time-1)
    }
  } else if (scen.ind==2){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0.5*(time-1)/(n-1)
    }
  } else if (scen.ind==3){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0.3*(time^2-1)/(n^2-1)
    }
  } else if (scen.ind==4){
    trt.diff = 0
    mu.time.func = function(time, seed){
      0.3*(time^3-1)/(n^3-1)
    }
  } else if (scen.ind==5){
    trt.diff = 0

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.003))
      time.temp = (c(0, cumsum(time.temp.1)))
      time.temp[time]
    }
    
  }   else if (scen.ind==6){
    trt.diff = 0

    mu.time.func = function(time, seed){
      set.seed(seed)
      time.temp.1 = rnorm(n = n-1, 0, sd = sqrt(0.005))
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
                          summary(lm.2.fit)$coefficients[2, 1:2]),
                        min(apply(out.cross, 1 , mean))
                        ))
  
      
  
      
    }

    p.cust.out.1 = cust.test.func(data.in = data.temp, method.ind=1)

    ## output
    return(c(
             p.cust.out.1[1]<=one.alpha,
             
             ## cust
             p.cust.out.1[2] - trt.diff,
             (p.cust.out.1[2]+qnorm(1-one.alpha)*p.cust.out.1[3]>=trt.diff)&
               (p.cust.out.1[2]-qnorm(1-one.alpha)*p.cust.out.1[3]<=trt.diff),
             p.cust.out.1[3]
             
             ))
  }
  stopCluster(cl)
  
  p.vec.out = matrix(unlist(p.vec),nrow = n.itt, ncol = 4, byrow = TRUE)
  
  # print(mean((p.vec.out[, 8])))
  # print(sd((p.vec.out[, 8])))
  # print(hist(p.vec.out[, 6]))
  
  table.name.out[scen.ind] = paste0(deparse(mu.time.func), collapse = "")
  table.out[scen.ind,] = apply(p.vec.out, 2, function(x){mean(x)})
  print(table.out)
}

print(Sys.time()-time.start)

colnames(table.out) = c( "cust.1",
                        "cust_bias",  "cust_cp",
                        "val_MSE")

write.csv(cbind(table.name.out, table.out), 
          paste0("Time_drift_sd1_", sd.1, "_sd2_", sd.2, 
                 "_trtprop_", trt.prop, 
                 ".csv", sep=""))















