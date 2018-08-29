setwd('~/d/sci/src/prnp_onset')
options(stringsAsFactors=FALSE)
library(survival)
library(sqldf)
library(reshape2)
library(binom)

top_hp_muts = c('P102L','D178N','E200K')
supp_hp_muts = c('5-OPRI','6-OPRI','P105L','A117V')
rapid_muts = c('E200K','D178N')
slow_muts = c('P102L','A117V','5-OPRI','6-OPRI','P105L')

mutparms = data.frame(mut=c('E200K','P102L','D178N','6-OPRI','5-OPRI','A117V','P105L'),
                      col=c('#7570B3','#1B9E77','#D95F02','#676752','#979770','#120112','#66A61E'),
                      pri=c('main','main','main','supp','supp','supp','supp'))
rownames(mutparms) = mutparms$mut

# coalesce courtesy of Kevin Jin: http://www.cureffi.org/2013/05/02/r-equivalent-of-sql-coalesce/#comment-2157150738
coalesce = function(..., default = NA) apply(cbind(..., default), 1, function(x) x[which(!is.na(x))[1]])

expand_range = function(x, by=.5) {
  return ( c(min(x)-by,max(x)+by) )
}

# percent function from https://github.com/macarthur-lab/exac_2015/blob/master/exac_constants.R
percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

fix_colnames = function (x) {
  return ( gsub('[^a-z0-9_]+','_',tolower(colnames(x))) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

# extract p value from a kaplan-meier test (survdiff)
km_pval = function(km_model) {
  chisq = km_model$chisq
  df = length(km_model$n) - 1
  return (1-pchisq(q=chisq, df=df))
}

smooth_gaussian = function(x, s=3, maxwidth=15) {
  smoothed_x = numeric(length(x))
  for (i in 1:length(x)) {
    total_density = 0
    xtemp = 0
    for (j in (i-floor(maxwidth/2)):(i+floor(maxwidth/2))) {
      if (j > 0 & j < length(x)) {
        total_density = total_density + dnorm(j-i, mean=0, sd=s)
        xtemp = xtemp + dnorm(j-i, mean=0, sd=s)*x[j]
      }
    }
    xtemp = xtemp / total_density
    smoothed_x[i] = xtemp
  }
  return (smoothed_x)
}

# abbreviations used are described in http://ocw.jhsph.edu/courses/FundEpi/PDFs/Lecture8.pdf
create_life_table = function(df) {
  lt = data.frame(age=1:100, lt=as.numeric(NA), dt=as.numeric(NA), wt=as.numeric(NA))
  for (i in 1:nrow(lt)) {
    lt$lt[i] = sum(df$surv_age >= lt$age[i], na.rm=T)
    lt$dt[i] = sum(df$surv_age == lt$age[i] & df$surv_status == 1, na.rm=T)
    lt$wt[i] = sum(df$surv_age == lt$age[i] & df$surv_status == 0, na.rm=T)
  }
  max_age = max(lt$age[lt$lt > 0])
  lt$lt1 = lt$lt - lt$wt/2 # avg number observed at midpoint of time interval
  lt$qt = lt$dt/lt$lt1 # P(event in this time interval)
  lt$qt[lt$lt1==0] = 0.0
  lt$qt[lt$age > max_age] = NA
  lt$qt_smooth[lt$age > max_age] = NA
  lt$qt_smooth[lt$age <= max_age] = smooth_gaussian(lt$qt[lt$age <= max_age])
  lt$pt = 0.0 # Pt = cumulative survival. after max_age it is 0, so may as well initialize with zero
  lt$pt[1] = 1.0
  for (i in 2:max_age) {
    lt$pt[i] = (1 - lt$qt[i]) * lt$pt[i-1]
  }
  lt$cm = NA # cm = conditional median, i.e. median age of onset conditioned on already surviving to age t
  for (i in 1:max_age) {
    lt$cm[i] = which(lt$pt / lt$pt[i] < .5)[1]
  }
  
  return (lt)
}

# function to undo the above one - useful to explore bootstrapping the confidence intervals
reverse_life_table = function(lt) {
  orig_data = data.frame(surv_age = integer(lt$lt[1]), surv_status = integer(lt$lt[1]))
  j = 1
  for (i in 1:nrow(lt)) {
    k = 0
    while (k < lt$dt[i]) {
      orig_data$surv_age[j] = lt$age[i]
      orig_data$surv_status[j] = 1
      j = j + 1
      k = k + 1
    }
    k = 0
    while (k < lt$wt[i]) {
      orig_data$surv_age[j] = lt$age[i]
      orig_data$surv_status[j] = 0
      j = j + 1
      k = k + 1
    } 
  }
  return (orig_data)
}




create_duration_table = function(df) {
  dt = data.frame(months=0:360, alive=as.numeric(NA), dead=as.numeric(NA), censored=as.numeric(NA))
  df$duration_mo = round(df$duration_d/30.44) # average number of days per month
  df = subset(df, !is.na(duration_mo) & !is.na(death_surv_status))
  for (i in 1:nrow(dt)) {
    dt$alive[i] = sum(df$duration_mo >= dt$months[i], na.rm=T)
    dt$dead[i] = sum(df$duration_mo == dt$months[i] & df$death_surv_status == 1) # 1 = had event (==dead) in this month
    dt$censored[i] = sum(df$duration_mo == dt$months[i] & df$death_surv_status == 0) # 0 = lost to followup, not dead, in this month
  }

  return (dt)
}



find_quartiles = function(lt) {
  if (class(lt) == 'data.frame') {
    tile25 = lt$age[which(lt$pt < c(.75))[1]]
    tile50 = lt$age[which(lt$pt < c(.50))[1]]
    tile75 = lt$age[which(lt$pt < c(.25))[1]]
  } else { # if they just pass the pt vector itself
    pt = lt
    tile25 = which(pt < .75)[1]
    tile50 = which(pt < .50)[1]
    tile75 = which(pt < .25)[1]
  }
  return (c(tile25, tile50, tile75))
}


# compute 95%CIs for survival & hazard

bootstrap_confidence_intervals = function (lt, n_iter = 100, seed=1) {
  # create columns for upper & lower 95%CI for qt (hazard) and cy (conditional years to onset)
  lt$qt_l95 = as.numeric(NA)
  lt$qt_u95 = as.numeric(NA)
  
  # Rstudio's SQLite doesn't have percentile functions so will do sorting & indexing manually...
  l95_index = ceiling(n_iter * .025)
  u95_index = floor(n_iter * .975)
  
  odata = reverse_life_table(lt)
  
  set.seed(seed) # set seed for random number generator so results are reproducible
  # n_iter = 100 # now a parameter
  max_age = max(lt$age[lt$lt > 0])
  permuted = expand.grid(age=1:max_age, iteration=1:n_iter)
  permuted$qt_smooth = as.numeric(NA)
  permuted$cy = as.numeric(NA)
  
  for (i in 1:n_iter) {
    cat(paste("\rBootstrapping confidence intervals, iteration ",i,"...",sep=''))
    flush.console()
    perm_data = odata[sort(sample(1:nrow(odata), size=nrow(odata), replace=T)),]
    perm_lt = create_life_table(perm_data)
    for (j in 1:max_age) {
      permuted$qt_smooth[permuted$iteration==i & permuted$age==j] = perm_lt$qt_smooth[j]
    }
  }
  
  cat(paste("\nSelecting confidence intervals...",sep=''))
  flush.console()
  
  permuted_sorted_qt = sqldf("
                             select   age, qt_smooth, iteration
                             from     permuted
                             order by 1 asc, 2 asc
                             ;")
  
  for (j in 1:max_age) {
    lt$qt_l95[j] = permuted_sorted_qt$qt_smooth[permuted_sorted_qt$age==j][l95_index]
    if (is.na(lt$qt_l95[j])) {
      lt$qt_l95[j] = min(permuted_sorted_qt$qt_smooth[permuted_sorted_qt$age==j], na.rm=T)
    }
    lt$qt_u95[j] = permuted_sorted_qt$qt_smooth[permuted_sorted_qt$age==j][u95_index]
    if (is.na(lt$qt_u95[j])) {
      lt$qt_u95[j] = max(permuted_sorted_qt$qt_smooth[permuted_sorted_qt$age==j], na.rm=T)
    }
  }
  cat(paste("\nDone!\n",sep=''))
  
  return (lt)
}

# convert hazards into cumulative survival
# (this cumulative survival is adjusted for # observed at each timepoint, unlike
# if you just divided number currently living by number that started out)
qt_to_pt = function(qt) {
  pt = numeric(length(qt))
  pt[1] = 1
  for (i in 2:length(qt)) {
    pt[i] = pt[i-1] * (1-qt[i])
  }
  return (pt)
}

surv_pval = function(survdiff_object) {
  return ( 1-pchisq(q=survdiff_object$chisq,df=length(survdiff_object$n)-1) )
}


# not using the below code currently.
# based on Lee & Wang 2003 section 4.2.2 clinical life tables - pp. 100 - 103 (PDF numbering)
formulaic_life_table = function(df) {
  lt = data.frame(age=1:100, # 1. interval - e.g. 1 is age 1 to 2
                  midpoint=1:100+0.5, # 2. midpoint (for plotting, etc)
                  b=rep(1,100), # 3. width - 1 year
                  w=integer(100), # 4/5. total "withdrawn alive or lost to followup" (alive and well or died of unrelated cause)
                  d=integer(100), # 6. number having onset/death in this interval
                  n_enter=integer(100), # 7. number entering this interval
                  n=integer(100), # 8. number exposed to risk
                  q=numeric(100), # 9. conditional proportion having onset/death in this interval
                  p=numeric(100), # 10. conditional proportion surviving
                  S=numeric(100), # 11. cumulative proportion surviving
                  f=numeric(100), # 12. probability density function
                  h=numeric(100), # 13. hazard function
                  varh = numeric(100), # Var(h) per Lee & Wang equation 4.2.13
                  h_l95 = numeric(100), # 95%CI lower bound on hazard
                  h_u95 = numeric(100) # 95%CI upper bound on hazard
                  )
  for (i in 1:nrow(lt)) {
    lt$w[i] = sum(df$surv_age == lt$age[i] & df$surv_status == 0, na.rm=T)
    lt$d[i] = sum(df$surv_age == lt$age[i] & df$surv_status == 1, na.rm=T)
    lt$n_enter[i] = sum(df$surv_age >= lt$age[i], na.rm=T)
  }
  lt$n = lt$n_enter - 0.5*lt$w
  lt$q = lt$d / lt$n
  lt$p = 1 - lt$q
  lt$S[1] = 1.0
  for (i in 2:nrow(lt)) {
    lt$S[i] = lt$p[i-1] * lt$S[i - 1]
  }
  lt$f = lt$S * lt$q / lt$b
  lt$h = 2*lt$q / (lt$b * (1 + lt$p))
  lt$varh = ( (lt$h^2) / (lt$n * lt$q) ) * (1 - (0.5*lt$h*lt$b)^2 )
  lt$varh[lt$q==0] = 0 # resolve division by zero errors. say no variance where no onsets
  lt$h_l95 = lt$h - 1.96 * sqrt(lt$varh) / sqrt(lt$n)
  lt$h_u95 = lt$h + 1.96 * sqrt(lt$varh) / sqrt(lt$n)
  return (lt)
}

