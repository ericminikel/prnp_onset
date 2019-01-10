setwd('~/d/sci/src/prnp_onset')
library(WriteXLS)
source('src/shared_code.R')
utils::assignInNamespace("format.pval",function(x,...){x},"base") # Konrad's trick to get rid of the silly 2.2e-16 thing

# code to support the claim in the text that "a CpG transition... occurs by spontaneous 
# mutation 10-100X more often than other mutation types"
mu = read.table('data/fordist_1KG_mutation_rate_table.txt',header=T) # downloaded from https://raw.githubusercontent.com/macarthur-lab/exac_2015/master/data/recurrence/fordist_1KG_mutation_rate_table.txt
cpg = (mu$from %in% c('ACG','CCG','GCG','TCG') & mu$to %in% c('ATG','CTG','GTG','TTG')) | # plus strand CpGs
  (mu$from %in% c('CGA','CGC','CGG','CGT') & mu$to %in% c('CAA','CAC','CAG','CAT')) # minus strand CpGs
mean(mu$mu_snp[cpg])/min(mu$mu_snp) # ratio CpG vs. least likely mutation - 119X
mean(mu$mu_snp[cpg])/max(mu$mu_snp[!cpg]) # ratio CpG vs. most likely non-CpG - 7.9X
mean(mu$mu_snp[cpg])/mean(mu$mu_snp[!cpg]) # ratio CpG vs. all other - 25X


# Figure S1: prevalence of highly penetrant mutations
allmutprev = read.table('data/mutation_prevalence_all.tsv',sep='\t',header=T) # segregation/de novo status and case count of all mutations
allmutprev$lhp = allmutprev$segregation=='x' | allmutprev$de_novo=='x' # which ones are likely high penetrance (LHP)
sum(allmutprev$lhp) # count unique variants with evidence for high penetrance = 27
sum(allmutprev$case_count[allmutprev$lhp]) # the 25 account for 1176 cases total in case series
mutprev = allmutprev[allmutprev$lhp & allmutprev$case_count > 0,c('variant','case_count')] # spin off a new table, subsetting to high penetrance variants with at least 1 case in the case series
mutprev = mutprev[order(-mutprev$case_count),] # sort by prevalence
mutprev$rank = 1:nrow(mutprev) # add a rank column
mutprev$proportion = mutprev$case_count / sum(mutprev$case_count) # proportion of cases of each variant
mutprev$cumulative_proportion = cumsum(mutprev$proportion) # cumulative proportion of cases of top N variants
write.table(mutprev, 'data/mutation_prevalence.tsv', sep='\t', row.names=F, col.names=T, quote=F) # save for later if needed

# code to support the claim that "The next four most common mutations... are also much rarer, accounting for only 10% of cases with a high penetrance variant"
sum(allmutprev$case_count[allmutprev$lhp]) # 0.1003

pdf('figures/figure_s1.pdf',width=8,height=5)
col_left = '#570067'
col_right = '#C70123'
par(mar=c(5,5,1,6))
plot(mutprev$rank+.5, mutprev$case_count, type='h', col=col_left, lwd=20, lend=3, 
     ylim=c(0,650), xlim=c(0,max(mutprev$rank)+1), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
axis(side=1, at=(1:max(mutprev$rank))+.5, labels=mutprev$variant, cex.axis=.8, las=2, lwd=0, line=-.5)
axis(side=1, at=c(1,max(mutprev$rank))+.5, labels=c(1,max(mutprev$rank)), line=2, lwd=0, lwd.ticks=0, font=2)
axis(side=2, at=(0:6)*100, lwd=0, lwd.ticks=1, las=2, col.axis=col_left, col.ticks=col_left)
par(new=TRUE)
plot(mutprev$rank+.5, mutprev$cumulative_proportion, type='l', col=col_right, lwd=5, lend=3, 
     ylim=c(0,1.05), xlim=c(0,max(mutprev$rank)+1), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
abline(h=0, lwd=1)
axis(side=4, at=(0:4)/4, labels=percent((0:4)/4), lwd.ticks=1, lwd=0, las=2, col.axis=col_right, col.ticks=col_right)
mtext(side=4, line=4.5, text='cumulative proportion of cases with\nhigh-penetrance variants', col=col_right, font=2)
mtext(side=1, line=3, text='rank prevalence', font=2)
mtext(side=2, line=3, text='number of cases', col=col_left, font=2)
#mtext(side=3, line=1, text='relative prevalence of high-penetrance\nvariants in prion disease cases', font=2, cex=1.2)
dev.off()

# --------- reconstruct limited version of original dataset, from life tables

data = data.frame(mut=character(0), age=integer(0), event=integer(0))
lt_list = {}
for (mut in c(top_hp_muts, supp_hp_muts)) {
  lt = read.table(paste('data/lt_',mut,'.tsv',sep=''), sep='\t', header=T)
  lt_list[[mut]] = lt
  mutdata = reverse_life_table(lt)
  mutdata$mut = mut
  data = rbind(data, mutdata[,c('mut','surv_age','surv_status')])
} 
colnames(data) = c('mut','age','event')


# summary statistics

# descriptive stats
ao_stats = sqldf("
                 select   mut variant,
                 sum (case when event is not null then 1 else 0 end) age_known
                 from     data
                 group by 1
                 order by 2 asc
                 ;")

ao_stats$mean = 0.0
ao_stats$sd = 0.0
ao_stats$median = 0
ao_stats$tile25 = 0
ao_stats$tile75 = 0

for (i in 1:nrow(ao_stats)) {
  ao_stats$mean[i] = mean(data$age[data$event==1 & data$mut==ao_stats$variant[i]], na.rm=T)
  ao_stats$sd[i] = sd(data$age[data$event==1 & data$mut==ao_stats$variant[i]], na.rm=T)
  ao_stats$onset_n[i] = sum(!is.na(data$age[data$event==1 & data$mut==ao_stats$variant[i]]))
  quartiles = find_quartiles(lt_list[[ao_stats$variant[i]]])
  ao_stats$median[i] = quartiles[2]
  ao_stats$tile25[i] = quartiles[1]
  ao_stats$tile75[i] = quartiles[3]
  ao_stats$surv_n[i] = max(lt_list[[ao_stats$variant[i]]]$lt)
  ao_stats$range_min[i] = min(data$age[data$event==1 & data$mut==ao_stats$variant[i]], na.rm=T)
  whichmax = which(data$mut==ao_stats$variant[i] & data$age==max(data$age[data$mut==ao_stats$variant[i]],na.rm=T) & !is.na(data$age)) 
  ao_stats$range_max[i] = data$age[whichmax]
  ao_stats$max_status[i] = data$event[whichmax]
}

ao_disp = ao_stats
ao_disp$meansd = paste(formatC(ao_disp$mean,format='f',digits=1), '+-', formatC(ao_disp$sd,format='f',digits=1))
ao_disp$surv = paste(ao_disp$median, ' (', ao_disp$tile25, ' - ', ao_disp$tile75, ')', sep='')
ao_disp$star = ifelse(ao_disp$max_status==0, '*', '')
ao_disp$range = paste(ao_disp$range_min, ' - ', ao_disp$range_max, ao_disp$star, sep='')
ao_disp = ao_disp[,c('variant','meansd','onset_n','surv','range','surv_n')]

ao_disp

# this is Table 1 and Table S3. For copying into Word manuscript it is easiest to export to Excel
ao_disp$meansd = gsub('\\+-','\U00B1',ao_disp$meansd)
ao_disp = ao_disp[match(c('P102L','D178N','E200K','5-OPRI','6-OPRI','P105L','A117V'),ao_disp$variant),]
WriteXLS(ao_disp, 'figures/table_1_table_s3.xls')

# for ppt version only - add in some more details
ao_disp$prop = percent(mutprev$proportion[match(ao_disp$variant,mutprev$variant)])
ao_disp$cumprop = percent(mutprev$cumulative_proportion[match(ao_disp$variant,mutprev$variant)])
for (i in 1:nrow(ao_disp)) {
  lt = lt_list[[ao_disp$variant[i]]]
  rowno = min(80, max(lt$age[!is.na(lt$pt)]))
  ao_disp$pen80[i] = percent(1-lt$pt[rowno])
}

ao_disp = ao_disp[with(ao_disp,order(cumprop)),]
ao_disp


### NAPOLEON PLOT
# Inspired by Charles Joseph Minard's plot of Napoleon's invasion of Russia
# https://www.edwardtufte.com/tufte/minard

napoleon_plot = function(mutlist) {
  
  par(mar=c(4,6,1,6))
  plot(NA,NA,axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',
       xlim=c(0,80), ylim=c(-.025,.30))
  axis(side=1, at=(0:8)*10, lwd=0, lwd.ticks=1)
  axis(side=2, at=(0:5)*.05, labels=percent((0:5)*.05), lwd=0, lwd.ticks=1, las=2)
  abline(h=0, lwd=2)
  segments(x0=0,x1=0,y0=0,y1=.25,lwd=2)
  abline(h=(0:5)*.05,lty=1,lwd=.25,col='#777777')
  mtext(side=1, line=2.5, text='age', font=1, cex=0.9)
  mtext(side=2, line=3.5, text='annual probability of disease onset', font=1, cex=0.9)
  
  for (mut in mutlist) {
    mutcol = mutparms$col[mutparms$mut==mut]
    cicol = paste(mutcol,'80',sep='') # 50% transparency for 95%CI blocks
    
    lt = lt_list[[mut]]
    max_age = max(lt$age[lt$lt > 0])
    
    max_thickness = 100*mutprev$proportion[mutprev$variant==mut]
    
    if (mut %in% top_hp_muts) {
      mtext(side=4, at=lt$qt_smooth[min(max_age,80)], las=2, text=mut, col=mutcol, font=2, line=1)    
    } else if (mut %in% supp_hp_muts) {
      label_age = min(lt$age[lt$qt_smooth > 0.025], na.rm=T) + 20
      text(x=label_age, y=.27, pos=3, srt=45, label=mut, col=mutcol, font=2, cex=0.6)
    }
    
    points(c(0,1), c(0,0), type='l', lwd=max_thickness, col=mutcol)
    for (i in 2:max_age) {
      hex_alpha = format(as.hexmode(round(lt$pt[i]*255)), width=2) # '20' 
      cicol = paste(mutcol,hex_alpha,sep='')
      points(lt$age[(i-1):i], lt$qt_smooth[(i-1):i], type='l', lwd=max_thickness*lt$pt[i], col=mutcol)
    }

  }
  
}



pdf('figures/figure_1.pdf', width=8, height=4)
napoleon_plot(c('E200K','P102L','D178N'))
thickness_legend = data.frame(prop=c(.01,.03,.10))
thickness_legend$lwd = 100 * thickness_legend$prop * max(mutprev$proportion)
legend(x=0,y=.255,xpd=T,legend=percent(thickness_legend$prop),lwd=thickness_legend$lwd,title='proportion surviving',bty='n',bg='white',horiz=T)
dev.off()



# an alternate version at same aspect ratio but with all 7 mutations
# I use this sometimes in presentations of this work
pdf('figures/figure_1_alternate.pdf', width=8, height=4)
napoleon_plot(c('E200K','P102L','D178N',supp_hp_muts))
dev.off()

# code to make hazard plots (for figure S3)

plot_smoothed_hazard = function(mut, xlims=c(0,80), display_mutation=FALSE) {
  
  mutcol = mutparms$col[mutparms$mut==mut]
  cicol = paste(mutcol,'80',sep='') # 50% transparency for 95%CI blocks
  
  lt = lt_list[[mut]]
  max_age = max(lt$age[lt$lt > 0])
  
  par(mar=c(4,5,3,1))
  plot(lt$age[1:max_age], lt$qt_smooth[1:max_age], type='l', lwd=5, col=mutcol, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i',
       xlim=xlims, ylim=c(0,.25))
  axis(side=1, at=(0:8)*10, lwd=0, lwd.ticks=1)
  axis(side=2, at=(0:5)*.05, labels=percent((0:5)*.05), lwd=0, lwd.ticks=1, las=2)
  abline(h=0, lwd=2)
  abline(v=0, lwd=2)
  abline(h=(0:5)*.05,lty=1,lwd=.25,col='#777777')
  points(lt$age[1:max_age], lt$qt_smooth[1:max_age], type='l', lwd=5, col=mutcol)
  polygon(x=c(lt$age[1:max_age], rev(lt$age[1:max_age])), y=c(lt$qt_u95[1:max_age], rev(lt$qt_l95[1:max_age])), col=cicol, border=NA)
  abline(v=find_quartiles(lt),lty=3,lwd=1,col='red')
  text(x=find_quartiles(lt), y=.25, pos=2, srt=90, col='red', labels=c('25th percentile','median onset','75th percentile'))
  mtext(side=1, line=2.5, text='age', font=1, cex=0.9)
  mtext(side=2, line=3.5, text='annual probability of onset', font=1, cex=0.9)
  if(display_mutation) {
    mtext(side=3, line=1, cex=1.2, font=2, text=mut, col=mutcol) 
  }
}

# --------- Figure S3: survival curve for all 7 mutations and hazards for 4 supplementary mutations

pdf('figures/figure_s3.pdf',width=9,height=11)
par(oma=c(2,2,5,2))

# layout:
# A = full-width napoleon plot, all 7 mutations
# B = 2/3-width survival curve, all 7 mutations
# C-E = 1/3 width hazard plots, top 3 mutations
# F-I = 1/4 width hazard plots, next 4 mutations

layout_matrix = matrix(c(rep(1,12),rep(2,8),rep(3:5,each=4),rep(6:9,each=4)),nrow=4,byrow=T)
layout(layout_matrix)
#layout(matrix(c(rep(1,8),2,3,2,3,2,3,4,5,4,5,4,5),nrow=10,byrow=T))

panel = 1

# A: Napoleon plot with all 7 mutations
napoleon_plot(c('E200K','P102L','D178N',supp_hp_muts))
mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 0.3)
panel = panel + 1

# B: Conventional survival curve with all 7 mutations

par(mar=c(7,5,3,1))
plot(NA,NA,xlim=c(0,100),ylim=c(0,1.05),xaxs='i',yaxs='i',ann=FALSE,axes=FALSE)
axis(side=1, at=(0:10)*10, lwd=0, lwd.ticks=1, line=2.5)
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1, las=2)
abline(h=0)
abline(v=0)
for (eamutation in c(top_hp_muts,supp_hp_muts)) {
  mutcol = mutparms$col[mutparms$mut==eamutation]
  this_lt = lt_list[[eamutation]]
  rowmax = which(this_lt$lt == 0)[1] - 1 # don't plot after there are no people left
  points(this_lt$age[1:rowmax], this_lt$pt[1:rowmax], type='l', lwd=5, col=mutcol)
  mutparms$median[mutparms$mut==eamutation] = find_quartiles(this_lt)[2]
  mutparms$max[mutparms$mut==eamutation] = max(this_lt$age[this_lt$lt > 0])
}
text(x=mutparms$max, y=rep(0,6), labels=mutparms$mut, col=mutparms$col, font=2, srt=-90, xpd=NA, adj=0)
mtext(side=1, line=4.5, text='age',font=1,cex=0.9)
mtext(side=2, line=3, text='proportion asymptomatic',font=1,cex=0.9)
#mtext(side=3, line=1, font=2, text='age of onset in genetic prion disease',cex=1.2)

mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 0.3)
panel = panel + 1

# C-E at 1/3 width, x axis up to 80
for (eamutation in c(top_hp_muts)) {
  plot_smoothed_hazard(eamutation, display_mutation=TRUE)
  mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 0.3)
  panel = panel + 1
}

# F-I at 1/4 width, x axis up to 80
for (eamutation in c(supp_hp_muts)) {
  plot_smoothed_hazard(eamutation, display_mutation=TRUE)
  mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 0.3)
  panel = panel + 1
}

dev.off() # -- end figure S3



# -------------------- survival power calculations

# Code here is ported from code for: http://www.sample-size.net/sample-size-survival-analysis/
# Provided by Michael Kohn.
# That code is in turn based on this citation: 
# Schoenfeld DA. Sample-size formula for the proportional-hazards regression model. Biometrics. 1983
# Jun;39(2):499-503. PubMed PMID: 6354290
# https://www.ncbi.nlm.nih.gov/pubmed/6 354290 

# Schoenfeld's equaton #1: 
n_events_required = function(alpha, beta, proportion_treated, hazard_ratio) {
  z1_a = qnorm(1-alpha/2) # /2 because two-tailed; not explained in Schoenfeld reference but noted here: http://www.sample-size.net/sample-size-survival-analysis/
  zb = qnorm(1-beta) # need qnorm(1-beta) or abs(qnorm(beta)) because should be positive - again not explained by Schoenfeld but see http://www.sample-size.net/sample-size-survival-analysis/
  pa = proportion_treated
  pb = 1 - proportion_treated
  events_required = ((zb + z1_a)^2) / (pa * pb * log(hazard_ratio)^2)
  return (events_required)
}

# test:
# n_events_required(.05, .20, .5, .5)

median_survival = function(hazard) {
  return (log(2) / hazard)
}

# "runin" is the number of years thrown out before data collection begins
# e.g. runin=1 means that onsets within the 1st year of the trial are not counted in the results.
n_enroll_required = function(alpha, beta, proportion_treated, hazard_ratio, baseline_hazard, censor_rate, followup, runin=0) {
  n_events = n_events_required(alpha, beta, proportion_treated, hazard_ratio)
  untreated_hazard = baseline_hazard
  treated_hazard = baseline_hazard * hazard_ratio
  treated_cumul_event_rate = (treated_hazard / (treated_hazard + censor_rate)) * (1 - exp(-(treated_hazard + censor_rate)*followup))
  treated_cumul_usable_event_rate = (treated_hazard / (treated_hazard + censor_rate)) * (exp(-(treated_hazard + censor_rate)*runin) - exp(-(treated_hazard + censor_rate)*followup))
  untreated_cumul_event_rate = (untreated_hazard / (untreated_hazard + censor_rate)) * (1 - exp(-(untreated_hazard + censor_rate)*followup))
  untreated_cumul_usable_event_rate = (untreated_hazard / (untreated_hazard + censor_rate)) * (exp(-(untreated_hazard + censor_rate)*runin) - exp(-(untreated_hazard + censor_rate)*followup))
  overall_cumul_event_rate = treated_cumul_event_rate * proportion_treated + untreated_cumul_event_rate * (1-proportion_treated)
  overall_cumul_usable_event_rate = treated_cumul_usable_event_rate * proportion_treated + untreated_cumul_usable_event_rate * (1-proportion_treated)
  n_total = n_events / overall_cumul_usable_event_rate
  n_treated = n_total * proportion_treated
  n_untreated = n_total * (1-proportion_treated)
  events_treated = n_treated * treated_cumul_usable_event_rate
  events_untreated = n_untreated * untreated_cumul_usable_event_rate
  return(c(n_treated, n_untreated, events_treated, events_untreated))
}


# some sanity checks
# n_enroll_required(.05, .20, .5, .2, .075, 0.30, 5)
# n_enroll_required(.05, .20, .5, .5, .075, 0.30, 5)
# n_enroll_required(.05, .20, .5, .8, .075, 0.30, 5)
# n_enroll_required(.05, .20, .5, .5, .056, 0.152, 5, runin=0)
# n_enroll_required(.05, .20, .5, .5, .056, 0.152, 5, runin=1)
# n_enroll_required(.05, .20, .9, .5, .056, 0.152, 15, runin=1)
# n_events_required(.05, .20, .9, .5)


average_hazard = function(mutations, min_age, max_age) {
  weighted_hazards = expand.grid(mut=mutations,age=min_age:max_age)
  weighted_hazards$hazard = NA
  weighted_hazards$weight = NA
  for (i in 1:nrow(weighted_hazards)) {
    # hazard is simply hazard for that mutation at that age
    weighted_hazards$hazard[i] = lt_list[[weighted_hazards$mut[i]]]$qt_smooth[weighted_hazards$age[i]]
    # weight is product of proportion of people with that mutation who survive to that age, and proportion of cases who have that mutation
    weighted_hazards$weight[i] = lt_list[[weighted_hazards$mut[i]]]$pt[weighted_hazards$age[i]] * mutprev$proportion[mutprev$variant==weighted_hazards$mut[i]]
  }
  # note use of na.rm=TRUE here. if no individuals have been observed to survive to a given age, then hazard is NA, 
  # but that means the weight is also ~0 so fine to remove them.
  overall_av_hazard = weighted.mean(x=weighted_hazards$hazard, w=weighted_hazards$weight, na.rm=TRUE)
  return (overall_av_hazard)
}


# some tests: 
# average_hazard(top_hp_muts, 18, 80)
# average_hazard(top_hp_muts, 40, 80)
# average_hazard(supp_hp_muts, 40, 80)


calculated_power_table = function(mutlist, wrate, min_age, max_age, runin, followup, alpha=.05, beta=.20, prop_trt=.5, avhaz=NA) {
  
  # if average hazard rate is not provided, calculate it from other params
  if (is.na(avhaz)) {
    avhaz = average_hazard(mutlist, min_age, max_age)
  }
  cat(paste('assumed average hazard: ',percent(avhaz),'\n'))
  
  # now construct a table
  calc_surv_power = expand.grid(hr=(1:9)/10, yrs=followup)
  calc_surv_power$n_enroll = as.numeric(NA)

  for (i in 1:nrow(calc_surv_power)) {
    enrollment_numbers = n_enroll_required(alpha, beta, prop_trt, calc_surv_power$hr[i], avhaz, wrate, calc_surv_power$yrs[i], runin)
    calc_surv_power$n_enroll[i] = round(sum(enrollment_numbers[1:2]))
  }
  
  calc_sp = dcast(calc_surv_power, hr ~ yrs, value.var='n_enroll')
  
  # figure out how many years' average delay in onset does each hazard ratio correspond to
  yrs_added = expand.grid(hr=calc_sp$hr, mut=mutlist)
  yrs_added$mut = as.character(yrs_added$mut) # otherwise it automatically becomes a factor, which messes things up
  yrs_added$proportion = mutprev$proportion[match(yrs_added$mut, mutprev$variant)]
  for (i in 1:nrow(yrs_added)) {
    lt = lt_list[[yrs_added$mut[i]]]
    pt_treated_c = qt_to_pt(lt$qt_smooth * yrs_added$hr[i])
    new_median = find_quartiles(pt_treated_c)[2]
    pt_untreated_c = qt_to_pt(lt$qt_smooth)
    old_median = find_quartiles(pt_untreated_c)[2]
    yrs_added$benefit[i] = as.numeric(new_median - old_median)
  }
  
  calc_sp$wmean_benefit = as.numeric(NA)
  for (i in 1:nrow(calc_sp)) {
    undefined_proportion = sum(yrs_added$proportion[yrs_added$hr==calc_sp$hr[i] & is.na(yrs_added$benefit)]) / sum(yrs_added$proportion[yrs_added$hr==calc_sp$hr[i]])
    use_rows = yrs_added$hr == calc_sp$hr[i] & !is.na(yrs_added$benefit)
    if (undefined_proportion < .5) {
      calc_sp$wmean_benefit[i] = round(sum(yrs_added$benefit[use_rows] * yrs_added$proportion[use_rows])/sum(yrs_added$proportion[use_rows]), digits=0)
    } else {
      calc_sp$wmean_benefit[i] = NA
    }
  }
  
  calc_sp$n_deaths = as.integer(NA)
  for (i in 1:nrow(calc_sp)) {
    deaths_required = n_events_required(alpha, beta, prop_trt, calc_sp$hr[i])
    calc_sp$n_deaths[i] = round(deaths_required)
  }
  
  colnames(calc_sp) = c('hr','fu5','benefit','deaths')
  column_order = c('hr','benefit','deaths','fu5')
  calc_sp = calc_sp[,column_order]
  
  return ( calc_sp )
  
}

# different withdrawal rate assumptions
# for basis for these, see Table S4...
# and Table S4is adapted from this blog post: http://www.cureffi.org/2017/03/03/preventive-trial-withdrawal-rates/
median_w = .152
min_w = .069
max_w = .549

# assumptions for min & max age from text
min_age = 40
max_age = 80

# Table 2
table_2 = calculated_power_table(mutlist=top_hp_muts, wrate=median_w, min_age=min_age, max_age=max_age, runin=1, followup=5)
write.table(table_2, 'figures/table_2.tsv', sep='\t', row.names=F, col.names=F, quote=F)
table_2$deaths = formatC(table_2$deaths,format='d',big.mark=',')
table_2$fu5 = formatC(table_2$fu5,format='d',big.mark=',')
WriteXLS(table_2,'figures/table_2.xls')

# Table S5
ts5_1 = calculated_power_table(mutlist=c(top_hp_muts,supp_hp_muts), wrate=min_w, min_age=min_age, max_age=max_age, runin=0, followup=5)
ts5_2 = calculated_power_table(mutlist=c(top_hp_muts), wrate=max_w, min_age=min_age, max_age=max_age, runin=1, followup=5, avhaz=average_hazard(top_hp_muts,min_age=min_age,max_age=max_age)*.75)
ts5_3 = calculated_power_table(mutlist=c('5-OPRI','6-OPRI','A117V'),  wrate=median_w, min_age=min_age, max_age=max_age, runin=1, followup=5)
ts5_4 = calculated_power_table(mutlist=top_hp_muts, wrate=median_w, min_age=min_age, max_age=max_age, runin=1, followup=15)
ts5_5 = calculated_power_table(mutlist=top_hp_muts, wrate=0.0, min_age=min_age, max_age=max_age, runin=1, followup=5)

table_s5 = rbind(ts5_1, ts5_2, ts5_3, ts5_4, ts5_5)
write.table(table_s5, 'figures/table_s5.tsv', sep='\t', row.names=F, col.names=F, quote=F)
table_s5$deaths = formatC(table_s5$deaths,format='d',big.mark=',')
table_s5$fu5 = formatC(table_s5$fu5,format='d',big.mark=',')
WriteXLS(table_s5,'figures/table_s5.xls')

# ------- begin code for simulation of clinical trials

generate_trial_outcome = function(
  years_followup = 5,
  max_age = 80,
  min_age = 40,
  proportion_treated = .5,
  hazard_ratio = .5,
  withdrawal_rate = .152,
  n_participants = 100,
  stratify = TRUE,
  mutations = top_hp_muts,
  mutdist = 'prevalence'
) {
  
  # sample a distribution of mutations among trial participants
  if (mutdist == 'prevalence') {
    participant_mutations = sample(mutprev$variant[mutprev$variant %in% mutations], prob=mutprev$proportion[mutprev$variant %in% mutations], size=n_participants, replace=TRUE)
  } else if (mutdist == 'uniform') {
    participant_mutations = sample(mutations, size=n_participants, replace=TRUE)
  }
  
  # sample a distribution of ages for those participants, conditioned on number of survivors at each age for each mutation
  starting_ages = rep(as.numeric(NA),n_participants)
  for (mut in mutations) {
    this_mut_indices = which(participant_mutations==mut)
    n_this_mut = length(this_mut_indices)
    lt = lt_list[[mut]]
    this_max_age = min(max_age, max(lt$age[!is.na(lt$pt)]))
    starting_ages[this_mut_indices] = sample(lt$age[min_age:this_max_age], replace=TRUE, size=n_this_mut, prob=lt$pt[min_age:this_max_age])
  }
  # alternate uniform method:
  #  starting_ages = round(runif(min=min_age,max=max_age,n=n_participants))
  
  # create a data frame of trial participants
  # use starting age, mutation, follow-up year, and follow-up event (1=onset, 0=censored)
  trial_participants = data.frame(starting_age = starting_ages, mutation = participant_mutations, fu_year=0, fu_event=as.integer(NA))
  
  # randomize who gets drug
  if (proportion_treated == .5) {
    if (stratify) {
      trial_participants$drug = NA
      for (mut in mutations) {
        n_this_mut = sum(trial_participants$mutation==mut)
        trial_participants$drug[trial_participants$mutation==mut] = rep(c(TRUE,FALSE),n_this_mut/2+1)[1:n_this_mut] # the +1)[1:n_this_mut] makes it robust to odd numbers
      }
    } else {
      trial_participants$drug = rep(c(TRUE,FALSE),n_participants/2+1)[1:n_participants] #sample(c(TRUE,FALSE),prob=c(.5,.5),replace=TRUE,size=nrow(trial_participants))
    }
  } else if (proportion_treated == 1) {
    trial_participants$drug = TRUE
  } else {
    stop("generate_trial_outcome not yet implemented for proportion_treated other than 0.5 or 1.0")
  }
  
  # for every year of the trial
  for (y in 1:years_followup) {
    # for every participant
    for (j in 1:nrow(trial_participants)) {
      if (!is.na(trial_participants$fu_event[j])) {
        next # if person already had onset or withdrew, skip them
      }
      # roll the dice to see if they withdraw this year
      if (sample(c(TRUE,FALSE),size=1,prob=c(withdrawal_rate, 1-withdrawal_rate))) {
        # and if yes, record the previous year of the trial as their follow-up year (since we won't have follow-up from _this_ year)
        # and censored (0) as their event
        trial_participants$fu_year[j] = y - 1
        trial_participants$fu_event[j] = 0
        # and then skip the rest of this iteration of the loop
        next
      }
      # ok now calculate this person's probability of onset during this year
      # based on life table for their mutation, at their starting age plus current year of the trial
      p_onset = lt_list[[trial_participants$mutation[j]]]$qt_smooth[trial_participants$starting_age[j] + y]
      # for simulations where some individuals survive beyond the last datapoint, fill in the hazard from the
      # highest observed age as the hazard thereafter
      if (is.na(p_onset)) {
        this_mut_max_age = max(lt_list[[trial_participants$mutation[j]]]$age[lt_list[[trial_participants$mutation[j]]]$lt > 0])
        p_onset = lt_list[[trial_participants$mutation[j]]]$qt_smooth[this_mut_max_age]
      }
      # if they get the drug, reduce their p_onset by the hazard ratio being simulated
      if (trial_participants$drug[j]) {
        p_onset = p_onset * hazard_ratio
      }
      # now roll the dice to see if they have onset this year
      if (sample(c(TRUE,FALSE),size=1,prob=c(p_onset, 1-p_onset))) {
        # and if yes, record this year of the trial as their onset year
        trial_participants$fu_year[j] = y
        trial_participants$fu_event[j] = 1
      }
    }
  }
  
  # end of trial! for people who are still A&W and participating at end of time period, set years
  # of followup and censored status
  alive_and_well_and_participating = trial_participants$fu_year == 0 & is.na(trial_participants$fu_event)
  trial_participants$fu_year[alive_and_well_and_participating] = years_followup
  trial_participants$fu_event[alive_and_well_and_participating] = 0
  
  return (trial_participants)
}

# some tests

generate_trial_outcome(
  years_followup = 5,
  max_age = 65,
  min_age = 40,
  proportion_treated = .5,
  hazard_ratio = .5,
  withdrawal_rate = .152,
  n_participants = 100,
  stratify = TRUE,
  mutations = top_hp_muts
)



simulate_randomized_trial = function(
  years_followup = 5,
  max_age = 80,
  min_age = 40,
  hazard_ratio = .5,
  proportion_treated = .5,
  withdrawal_rate = .152,
  n_participants = 100,
  n_iter = 100,
  alpha = .05,
  stratify = TRUE,
  mutations = top_hp_muts,
  mutdist = 'prevalence',
  runin = 0
) {
  
  n_positive_results = 0 # this will be incremented in the loop
  
  for (i in 1:n_iter) {
    
    cat(paste("\rcurrently on iteration: ",i,sep=''))
    
    trial_participants = generate_trial_outcome(years_followup, max_age, min_age, proportion_treated, hazard_ratio, withdrawal_rate, n_participants, stratify, mutations, mutdist)
    
    trial_participants = subset(trial_participants, fu_year > runin) # apply run-in period
    
    # this if/else clause handles special cases to make this code robust to parameters with very low N or low HR
    # as well as stratified vs. non-stratified trials
    if (length(trial_participants)==0) { # if all participants withdrew
      # need to skip the cox test because result is undefined and will throw error
      next # skip to next iteration of loop, thus skipping the n_positive_results++ step
    } else if (length(unique(trial_participants$drug)) < 2) { # if all drug or all placebo participants withdrew
        # need to skip the cox test because result is undefined and will throw error
        next # skip to next iteration of loop, thus skipping the n_positive_results++ step
    } else if (sum(trial_participants$fu_event)==0) { # if no one had onset during the trial, trial failed.
      # need to skip the cox test because result is undefined and will throw error
      next # skip to next iteration of loop, thus skipping the n_positive_results++ step
    } else if (!stratify) { # non-stratified trial
      # simply do a log-rank test to see if the drug affected survival
      km_model = survdiff(Surv(fu_year,fu_event)~drug,data=trial_participants)
      drug_pval = km_pval(km_model)
    } else if (stratify) { # stratified trial
      # build a Cox Proportional Hazards model to test if the drug affected survival, while controlling for mutation
      # this requires a try/catch because sometimes at low N, there will be no onsets for one or more groups
      # (mutation+treatment combo), and this makes the Cox result undefined and can throw
      # warnings or errors. This is one weakness of the stratified model as explained in
      # supplementary discussion
      possibleError = tryCatch({
        survival_model = coxph(Surv(fu_year,fu_event)~drug+mutation,data=trial_participants)
      }, warning = function(w) { w
      }, error = function (e) { e
      }, finally = { 
      })
      if (inherits(possibleError, 'error')) {
        next 
      } else {
        drug_pval = summary(survival_model)$coefficients["drugTRUE","Pr(>|z|)"]
      }
        
    } 
      
    # if p value is less than specified alpha threshold, you have a positive result
    if (drug_pval < alpha) {
      n_positive_results = n_positive_results + 1
    }
  }
  
  cat(paste("\rdone!\n",sep=''))
  power = n_positive_results / n_iter
  #print(survival_model) # debugging
  #print(table(trial_participants[,c('mutation','drug','fu_event')]))
  
  
  return (power)
}

# test it out - example call
# simulate_randomized_trial(
#   years_followup = 5,
#   max_age = 80,
#   min_age = 40,
#   hazard_ratio = .1,
#   withdrawal_rate = .152,
#   n_participants = 500,
#   mutations = top_hp_muts,
#   alpha = .05,
#   n_iter = 500,
#   stratify = TRUE,
#   runin = 1
# )

# ---- Table S6
regenerate_table_s6 = FALSE # this takes hours so by default don't redo it every time the script is run
if (regenerate_table_s6) {
  set.seed(1)
  table_s6 = table_2 # start from table 2
  table_s6$fu5 = as.integer(gsub(',','',table_s6$fu5)) # convert follow-up years back to integer
  table_s6$calc_power = .80
  table_s6$sim_power_nostrat = as.numeric(NA)
  table_s6$sim_power_strat = as.numeric(NA)
  n_iter = 500
  for (i in 1:9) {
    
    cat(paste("\rnow re-generating table S6 row: ",i,'               \n',sep='')) # 
    flush.console()
    
    table_s6$sim_power_nostrat[i] = simulate_randomized_trial(years_followup = 5, max_age=80, min_age=40, hazard_ratio = table_s6$hr[i], proportion_treated = 0.5, withdrawal_rate = median_w, 
                                                            n_participants = table_s6$fu5[i], n_iter = n_iter, alpha = 0.05, stratify = FALSE, mutations = top_hp_muts, runin = 1)
    table_s6$sim_power_strat[i] = simulate_randomized_trial(years_followup = 5, max_age=80, min_age=40, hazard_ratio = table_s6$hr[i], proportion_treated = 0.5, withdrawal_rate = median_w, 
                                                          n_participants = table_s6$fu5[i], n_iter = n_iter, alpha = 0.05, stratify = TRUE, mutations = top_hp_muts, runin = 1)
  }
  table_s6

  write.table(table_s6, 'figures/table_s6.tsv', sep='\t', row.names=F, col.names=T, quote=F)
  WriteXLS(table_s6, 'figures/table_s6.xls')
}

# --- to explain why stratification doesn't help, how much variance is explained by mutation?
m = lm(age ~ mut, data = subset(data, event == 1 & data$mut %in% top_hp_muts))
summary(m)




simulate_natural_history_trial = function(
  years_followup = 5,
  max_age = 80,
  min_age = 40,
  hazard_ratio = .5,
  withdrawal_rate = .152,
  n_participants = 100,
  mutations = top_hp_muts,
  alpha = .05,
  n_iter = 100,
  stratify_cox = TRUE,
  runin = 1 # run-in period, the number of initial years of data to ignore
) {
  proportion_treated = 1 # always for this natural history design
  stratify = FALSE # always for this natural history design (the analysis can be stratified but there is no randomization in assigning drug status)
  n_positive_results = 0
  for (i in 1:n_iter) {
    cat(paste("\rcurrently on iteration: ",i,'               ',sep='')) # extra spaces to blank out leftover digits from earlier messages
    flush.console()
    tp = generate_trial_outcome(years_followup, max_age, min_age, proportion_treated, hazard_ratio, withdrawal_rate, n_participants, stratify, mutations)
    
    tp$end_age = tp$starting_age + tp$fu_year 
    tp$event = tp$fu_event
    
    
    # handle run-in period, e.g. if runin == 1, then only consider data from year 2 onward
    tp$starting_age = tp$starting_age + runin # ignore 
    tp = subset(tp, end_age >= starting_age + 1) # + 1 because for cox left-truncated model, cannot include 1st year withdrawals in analysis
    nh = subset(data, !is.na(age) & mut %in% mutations)
    nh$mutation = nh$mut
    nh$end_age = nh$age
    nh$starting_age = 0
    
    nh$drug = FALSE
    tpnh = rbind(tp[,c('mutation','starting_age','end_age','event','drug')], nh[,c('mutation','starting_age','end_age','event','drug')])
    if (nrow(tp)==0) { # if all participants withdrew
      # need to skip the cox test because result is undefined and will throw error
      next # skip to next iteration of loop, thus skipping the n_positive_results++ step
    } else if (stratify_cox) { # Cox PH model stratified by mutation
      coxobj = coxph(Surv(time=starting_age,time2=end_age,event=event,type='counting')~drug+mutation,data=tpnh)
    } else { # non-stratified Cox model
      coxobj = coxph(Surv(time=starting_age,time2=end_age,event=event,type='counting')~drug,data=tpnh)
    }
    p = summary(coxobj)$coefficients['drugTRUE','Pr(>|z|)']
    if (p < alpha) {
      n_positive_results = n_positive_results + 1
    }
  }
  cat("done!\n")
  power = n_positive_results / n_iter
  return(power)
}

redo_costly_analyses = FALSE
if (redo_costly_analyses) {
  
  # ---- test for bias
  # run a simulation with a large N but with HR==1 (i.e. an ineffective drug)
  # and the power (1-beta) should be ~= the statistical threshold (alpha)
  set.seed(1)
  simulate_natural_history_trial(
    years_followup = 5,
    max_age = 80,
    min_age = 40,
    hazard_ratio = 1,
    withdrawal_rate = .152,
    n_participants = 1000,
    mutations = top_hp_muts,
    alpha = .05,
    n_iter = 1000,
    stratify_cox = TRUE,
    runin = 1
  )
  
  
  # check if stratification matters at realistic N
  set.seed(1)
  simulate_natural_history_trial(
    years_followup = 15,
    max_age = 80,
    min_age = 40,
    hazard_ratio = 0.5,
    withdrawal_rate = median_w,
    n_participants = 156,
    mutations = top_hp_muts,
    alpha = .05,
    n_iter = 1000,
    stratify_cox = FALSE,
    runin = 1
  )
  # 90.6%
  simulate_natural_history_trial(
    years_followup = 15,
    max_age = 80,
    min_age = 40,
    hazard_ratio = 0.5,
    withdrawal_rate = median_w,
    n_participants = 156,
    mutations = top_hp_muts,
    alpha = .05,
    n_iter = 1000,
    stratify_cox = TRUE,
    runin = 1
  )
  # 94.1%
  binom.confint(906,1000,method='wilson')
  binom.confint(941,1000,method='wilson')
}



explore_nh_params = function(simparams, max_age = 80, min_age = 40, mutations = top_hp_muts, alpha = 0.05, runin = 1, n_iter = 100) {
  set.seed(1) # to make results of this function exactly reproducible
  simresults = simparams
  simresults$power = as.numeric(NA)
  cat(paste("number of scenarios to simulate: ",nrow(simresults),"\n"))
  for (i in 1:nrow(simresults)) {
    cat(paste("scenario ",i,"/",nrow(simresults)," -- current parameters: ",paste(simparams[i,],collapse=' '),"\n",sep=''))
    flush.console()
    simresults$power[i] = simulate_natural_history_trial(years_followup=simresults$followup[i],
                                                            max_age = max_age, min_age = min_age, withdrawal_rate = simresults$w[i],
                                                            hazard_ratio = simresults$hazard_ratio[i],
                                                            n_participants = simresults$n[i],
                                                            n_iter = n_iter,
                                                            alpha = alpha,
                                                            stratify_cox = simresults$stratify[i],
                                                            mutations = mutations,
                                                            runin = runin)
  }
  return (simresults)
  
}


#### Begin Table 3

regenerate_table3_data = FALSE
n_iter = 2000
if (regenerate_table3_data) {
  # iterate to figure out approximately where 80% power is reached; manually inspect results
  table3_params = rbind(expand.grid(hazard_ratio = 0.5, followup = c(5,15), stratify = FALSE, n=(5:30)*10, w=c(median_w)),
                        expand.grid(hazard_ratio = 0.5, followup = c(15), stratify = FALSE, n=(1:20)*10, w=c(0.0)))
  table3_data = explore_nh_params(table3_params, n_iter = n_iter)
  write.table(table3_data, 'figures/table3_data_temp.tsv', sep='\t', row.names=F, col.names=T, quote=F)
  # come back and do a finer-grained search right around the 80% power point for each 
  table3_addlparams = rbind(expand.grid(hazard_ratio = 0.5, followup = 5, stratify = FALSE, n=221:229, w=c(median_w)),
                            expand.grid(hazard_ratio = 0.5, followup = 15, stratify = FALSE, n=120:130, w=c(median_w)),
                            expand.grid(hazard_ratio = 0.5, followup = 15, stratify = FALSE, n=30:40, w=c(0.0)))
  table3_addldata = explore_nh_params(table3_addlparams, n_iter = n_iter)
  # combine the two, sort by followup and w and write to disk
  table3 = rbind(table3_data, table3_addldata)
  table3 = table3[with(table3, order(followup,-w,n)),]
  write.table(table3, 'figures/table3_data.tsv', sep='\t', row.names=F, col.names=T, quote=F)
} else {
  table3_data = read.table('figures/table3_data.tsv', sep='\t', header=T)
}



##### Begin Figure S6

regenerate_fig_s6_data = FALSE

## Figure S6A
# power over time for 2 N values, depending whether withdrawals are accounted for
n_iter = 2000
if (regenerate_fig_s6_data) {
  figs6a_params = expand.grid(hazard_ratio = 0.5, followup = 2:15, stratify = FALSE, n=c(156, 60), w=c(median_w,0.0))
  figs6a_data = explore_nh_params(figs6a_params, n_iter = n_iter)
  write.table(figs6a_data, 'figures/figure_s6a_data.tsv', sep='\t', row.names=F, col.names=T, quote=F)
} else {
  figs6a_data = read.table('figures/figure_s6a_data.tsv', sep='\t', header=T)
}



plotparms = unique(figs6a_data[,c("n","w")])
hicol = "#31A354"
locol = "#ADDD8E"
plotparms$col = ""
plotparms$col[plotparms$n==156] = hicol
plotparms$col[plotparms$n==60] = locol
plotparms$lty = as.integer(NA)
plotparms$lty[plotparms$w==median_w] = 1
dashed = 3
plotparms$lty[plotparms$w==0] = dashed

pdf('figures/figure_s6.pdf',width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,1,1),oma=c(1,2,1,1))

plot(NA, NA, xlim=c(0,15), ylim=c(0,1.1), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
axis(side=1, at=(0:4)*5, lwd=0, lwd.ticks=1)
axis(side=2, at=(0:10)/10, labels=percent((0:10)/10), lwd=0, lwd.ticks=1, las=2)
abline(h=0, lwd=2)
abline(v=0, lwd=2)
abline(h=1, lwd=1)
abline(h=.8, lwd=1, lty=2)
for (i in 1:nrow(plotparms)) {
  subs = subset(figs6a_data, n == plotparms$n[i] & w == plotparms$w[i])
  # because the run-in period is 1 year, there is 0 power after 1 year. manually add these points to the curve:
  points(c(1,subs$followup), c(0,subs$power), type='l',  lwd=3, lty=plotparms$lty[i], col=plotparms$col[i])
}
mtext(side=1, line=2.5, text='years of follow-up')
mtext(side=2, line=3, text='statistical power')
legend('bottomright',c('N = 156','N = 60','15.2% withdrawal','no withdrawal'),col=c(hicol, locol,"#000000",'#000000'),pch=c(15,15,NA,NA),lwd=c(0,0,3,3),lty=c(NA,NA,1,dashed), bty='n', cex=0.7)

mtext('A', side=3, cex=2, adj = -0.2, line = 0.3)
# End Figure S6A


### Figure S6B
# Number of people needed as a function of hazard ratio - three curves for three trial designs
n_iter = 100
n_to_simulate = c(seq(10,100,10),seq(120,300,30),seq(400,1000,100)) 
if (regenerate_fig_s6_data) {
  figs6b_params = expand.grid(hazard_ratio = 1:9/10, followup = 15, stratify = FALSE, n=n_to_simulate, w=c(median_w,0.0))
  figs6b_data = explore_nh_params(figs6b_params, n_iter = n_iter)  
  write.table(figs6b_data, 'figures/figure_s6b_data.tsv', sep='\t', row.names=F, col.names=T, quote=F)
} else {
  figs6b_data = read.table('figures/figure_s6b_data.tsv', sep='\t', header=T)
}



with_w = sqldf("
               select   hazard_ratio, min(n) min_n
               from     figs6b_data
               where    w = 0.152
               and      power > 0.80
               and      hazard_ratio < 0.8
               group by 1
               order by 1
               ;")

no_w = sqldf("
             select   hazard_ratio, min(n) min_n
             from     figs6b_data
             where    w = 0
             and      power > 0.80
             and      hazard_ratio < 0.8
             group by 1
             order by 1
             ;")


col_random = '#E34A33'
col_nh_w = '#253494'
col_nh_no = '#41B6C4'

# Begin plotting Figure S6B
plot(NA, NA, xlim=c(0,1.3), ylim=c(0,1000), axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=(0:10)/10, labels=NA, lwd=1, lwd.ticks=1)
axis(side=1, at=c(0,.5,1), labels=c("0.0","0.5","1.0"), lwd=0, lwd.ticks=1, line=0.5)
axis(side=2, at=(0:10)*100, lwd=0, lwd.ticks=1, las=2)
#abline(h=0,lwd=2)
abline(v=0,lwd=2)
num_t2 = table_2
num_t2$fu5 = as.integer(gsub(',','',num_t2$fu5))
points(num_t2$hr, num_t2$fu5, col=col_random, type='l', lwd=3)
text(num_t2$hr[5], num_t2$fu5[5], pos=4, labels='randomized 5-year', col=col_random, font=2, cex=0.8)
points(with_w$hazard_ratio, with_w$min_n, col=col_nh_w, type='l', lwd=3)
text(max(with_w$hazard_ratio), max(with_w$min_n), pos=4, labels='post-marketing\n15-year', col=col_nh_w, font=2, cex=0.8)
points(no_w$hazard_ratio, no_w$min_n, col=col_nh_no, type='l', lwd=3)
text(max(no_w$hazard_ratio), max(no_w$min_n), pos=4, labels='post-marketing\n15-year\nno withdrawal', col=col_nh_no, font=2, cex=0.8)
mtext(side=1, line=2.5, at=0.5, text='hazard ratio')
mtext(side=2, line=3, text='participants required')
#legend('topright',c('randomized 5-year\n','post-marketing 15 year\n','post-marketing 15 year\nno withdrawal'),
#       col=c(col_random, col_nh_w, col_nh_no), lwd=3, bg = 'white', bty='n', cex=0.8)


mtext('B', side=3, cex=2, adj = -0.2, line = 0.3)
# End Figure S6B



dev.off()