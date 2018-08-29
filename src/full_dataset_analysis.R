setwd('~/d/sci/src/prnp_onset')
source('src/shared_code.R')
library(WriteXLS)

mutprev = read.table('data/mutation_prevalence.tsv',sep='\t',header=T)

ucsf = read.table('../prnp_onset_private/data/ucsf.tsv',sep='\t',header=T,quote='',comment.char='')
german = read.table('../prnp_onset_private/data/german.tsv',sep='\t',header=T,quote='',comment.char='')
mrc = read.table('../prnp_onset_private/data/mrc.tsv',sep='\t',header=T,quote='',comment.char='')
e200k = read.table('../prnp_onset_private/data/e200k-updated.tsv',sep='\t',header=T,quote='',comment.char='')
bologna = read.table('../prnp_onset_private/data/bologna.tsv',sep='\t',header=T,quote='',comment.char='')
veneto = read.table('../prnp_onset_private/data/veneto.tsv',sep='\t',header=T,quote='',comment.char='')
japan = read.table('../prnp_onset_private/data/japan.tsv',sep='\t',header=T,quote='',comment.char='')
france = read.table('../prnp_onset_private/data/france.tsv',sep='\t',header=T,quote='',comment.char='')
spain = read.table('../prnp_onset_private/data/spain.tsv',sep='\t',header=T,quote='',comment.char='')


data = rbind(ucsf, german, mrc, e200k, bologna, veneto, japan, france, spain)

# some final data cleaning
data$questionable[is.na(data$questionable)] = 0 # NA or blank is equivalent to "0", not questionable
data$sex[data$sex=='1'] = 'M'
data$sex[data$sex=='2'] = 'F'
data$sex[data$sex=='-1'] = '?'

data$cis129 = gsub(' ','',data$cis129)
data$trans129 = gsub(' ','',data$trans129)

# check on predictive tests
sum(data$predictive==1, na.rm=T)
sum(data$predictive==1 & data$prion_disease_status == 1, na.rm=T)
sum(data$predictive==1 & data$study=='MRC', na.rm=T)
sum(data$predictive==1 & data$study=='MRC' & data$prion_disease_status > 1, na.rm=T)
data[data$predictive==1 & data$study=='MRC' & data$prion_disease_status > 1 & !is.na(data$predictive),]

# first of all, how many people in this dataset have a positive predictive testing result
# for a high penetrance mutation and are currently asymptomatic?
sum(data$predictive==1 & data$mutant_alleles==1 & data$prion_disease_status==1 & data$family_mutation %in% mutprev$variant, na.rm=T)
sum(data$predictive==1 & data$mutant_alleles==1 & data$prion_disease_status==1 & data$family_mutation %in% mutprev$variant & data$surv_age >= 40, na.rm=T)

# remove any cases where some data are uncertain
data = data[data$questionable==0,]
# subset to mutations of interest
data = data[data$family_mutation %in% c(top_hp_muts,supp_hp_muts), ]
# remove people without mutations from the dataset
data = data[data$mutant_alleles==1, ]

# choose survival ages -- onset or death for rapid mutations, onset only for slowly progressive mutations
data$surv_age[data$family_mutation %in% rapid_muts] = coalesce(data$onset_surv_age[data$family_mutation %in% rapid_muts], data$death_surv_age[data$family_mutation %in% rapid_muts])
data$surv_status[data$family_mutation %in% rapid_muts] = pmax(data$onset_surv_status[data$family_mutation %in% rapid_muts], data$death_surv_status[data$family_mutation %in% rapid_muts], na.rm=TRUE)
data$surv_age[data$family_mutation %in% slow_muts] = data$onset_surv_age[data$family_mutation %in% slow_muts]
data$surv_status[data$family_mutation %in% slow_muts] = data$onset_surv_status[data$family_mutation %in% slow_muts]

# write out life tables -- these will be the input to public_dataset_analysis.R and will also
# be released as supplementary tables.
lt_list = {}
regenerate_life_tables = FALSE
if (regenerate_life_tables) {
  
  # generate life tables from raw data
  for (mut in c(top_hp_muts, supp_hp_muts)) {
    cat(paste("Generating life table for mutation ",mut,"...\n",sep=""))
    flush.console()
    lt = create_life_table(subset(data, family_mutation==mut))
    lt = bootstrap_confidence_intervals(lt, n_iter = 1000)
    lt_list[[mut]] = lt
    write.table(lt, paste('data/lt_',mut,'.tsv',sep=''), sep='\t', row.names=F, col.names=T, quote=F)
  } 
  
  # write out as XLS for the Supplementary Life Tables
  supp_life_tables = {}
  for (mut in c(top_hp_muts, supp_hp_muts)) {
    lt = lt_list[[mut]][,c('age','lt','dt','wt','qt','qt_smooth','pt','cm')]
    colnames(lt) = c('age','lives','onsets','withdrawals','raw_hazard','smoothed_hazard','cumulative_survival','conditional_median')
    supp_life_tables[[mut]] = format(lt, digits=2)
  }
  supp_life_tables[['definitions']] = data.frame(
    column_name = c('age','lives','onsets','withdrawals','raw_hazard','smoothed_hazard','cumulative_survival','conditional_median'),
    definition = c('Age in years',
                   'Number of healthy individuals at the beginning of this age',
                   'Number of individuals with onset or death at this age',
                   'Number of individuals lost to followup, alive and well, or dead of unrelated causes at this age',
                   'Probability of onset at this age',
                   'Smoothed probability of onset at this age',
                   'Proportion of individuals still alive and well at this age',
                   'Median expected age of onset conditional on being alive and well at this age')
  )
  WriteXLS(supp_life_tables, 'figures/supplementary_life_tables.xls', SheetNames=c(top_hp_muts, supp_hp_muts,'definitions'))
  # note - couldn't figure out how to set formatting in Excel programmatically
  # so will still need to open and manually set bold, centered, number formats etc
  
} else {
  
  # if not re-generating the life tables then read them in from disk
  for (mut in c(top_hp_muts, supp_hp_muts)) {
    lt = read.table(paste('data/lt_',mut,'.tsv',sep=''), sep='\t', header=T)
    lt_list[[mut]] = lt
  } 
  
}


regenerate_duration_tables = FALSE
if (regenerate_duration_tables) {
  
  # create the tables from the raw data
  supp_duration_tables = {}
  for (mut in c(top_hp_muts, supp_hp_muts)) {
    supp_duration_tables[[mut]] = create_duration_table(subset(data, family_mutation==mut))
    write.table(supp_duration_tables[[mut]], paste('data/duration_',mut,'.tsv',sep=''), sep='\t', row.names=F, col.names=T, quote=F)
  }
  
  # write out to Excel
  supp_duration_tables[['definitions']] = data.frame(
    column_name = c('months','alive','dead','censored'),
    definition = c('Months since first symptom',
                   'Number of individuals symptomatic but alive at the beginning of the time period',
                   'Number of individuals who died this number of months after symptoms',
                   'Number of individuals who were symptomatic but alive at last followup at this number of months')
  )
  WriteXLS(supp_duration_tables, 'figures/supplementary_duration_tables.xls', SheetNames=c(top_hp_muts, supp_hp_muts,'definitions'))
  
}


desc_stats_by_study = sqldf("
                          select   study,
                                   sum(case when surv_age is not null then 1 else 0 end) age_known
                          from     data
                          where    mutant_alleles == 1
                          and      family_mutation in (select mut from mutparms)
                          group by 1
                          order by 2 desc
                          ;")
desc_stats_by_study

desc_stats_by_ascertainment = sqldf("
                                    select   ascertainment,
                                             sum(case when surv_age is not null then 1 else 0 end) age_known
                                    from     data
                                    where    mutant_alleles == 1
                                    and      family_mutation in (select mut from mutparms)
                                    group by 1
                                    order by 2 desc
                                    ;")
desc_stats_by_ascertainment

desc_stats_by_vital_status = sqldf("
                                    select   prion_disease_status,
                                    sum(case when surv_age is not null then 1 else 0 end) age_known
                                    from     data
                                    where    mutant_alleles == 1
                                    and      family_mutation in (select mut from mutparms)
                                    group by 1
                                    order by 1 asc
                                    ;")
desc_stats_by_vital_status

table_s2_colnames = c('source','n')
table_s2 = rbind(setNames(desc_stats_by_study, table_s2_colnames), setNames(data.frame(cbind("total",sum(desc_stats_by_study$age_known))),table_s2_colnames),
                 setNames(desc_stats_by_ascertainment, table_s2_colnames), setNames(data.frame(cbind("total",sum(desc_stats_by_ascertainment$age_known))),table_s2_colnames),
                 setNames(desc_stats_by_vital_status, table_s2_colnames), setNames(data.frame(cbind("total",sum(desc_stats_by_vital_status$age_known))),table_s2_colnames))

write.table(table_s2, 'figures/table_s2.tsv', sep='\t', row.names=F, col.names=T, quote=F)


# descriptive stats
ao_stats = sqldf("
                 select   family_mutation variant,
                 sum(case when surv_age is not null then 1 else 0 end) age_known
                 from     data
                 where    mutant_alleles == 1
                 and      family_mutation in (select mut from mutparms)
                 group by 1
                 order by 2 asc
                 ;")

ao_stats$mean = 0.0
ao_stats$sd = 0.0
ao_stats$median = 0
ao_stats$tile25 = 0
ao_stats$tile75 = 0

for (i in 1:nrow(ao_stats)) {
  ao_stats$mean[i] = mean(data$surv_age[data$surv_status==1 & data$family_mutation==ao_stats$variant[i]], na.rm=T)
  ao_stats$sd[i] = sd(data$surv_age[data$surv_status==1 & data$family_mutation==ao_stats$variant[i]], na.rm=T)
  ao_stats$onset_n[i] = sum(!is.na(data$surv_age[data$surv_status==1 & data$family_mutation==ao_stats$variant[i]]))
  quartiles = find_quartiles(lt_list[[ao_stats$variant[i]]])
  ao_stats$median[i] = quartiles[2]
  ao_stats$tile25[i] = quartiles[1]
  ao_stats$tile75[i] = quartiles[3]
  ao_stats$surv_n[i] = max(lt_list[[ao_stats$variant[i]]]$lt)
  ao_stats$range_min[i] = min(data$surv_age[data$surv_status==1 & data$family_mutation==ao_stats$variant[i]], na.rm=T)
  whichmax = which(data$family_mutation==ao_stats$variant[i] & data$surv_age==max(data$surv_age[data$family_mutation==ao_stats$variant[i]],na.rm=T) & !is.na(data$surv_age)) 
  ao_stats$range_max[i] = data$surv_age[whichmax]
  ao_stats$max_status[i] = data$surv_status[whichmax]
}

ao_disp = ao_stats
ao_disp$meansd = paste(formatC(ao_disp$mean,format='f',digits=1), '+-', formatC(ao_disp$sd,format='f',digits=1))
ao_disp$surv = paste(ao_disp$median, ' (', ao_disp$tile25, ' - ', ao_disp$tile75, ')', sep='')
ao_disp$star = ifelse(ao_disp$max_status==0, '*', '')
ao_disp$range = paste(ao_disp$range_min, ' - ', ao_disp$range_max, ao_disp$star, sep='')
ao_disp = ao_disp[,c('variant','meansd','onset_n','surv','range','surv_n')]

ao_disp

sum(ao_disp$surv_n[ao_disp$variant %in% top_hp_muts])




# --- figure S2. duration by top 3 mutations
pdf('figures/figure_s2.pdf',width=8,height=8)


par(mfrow=c(2,1),mar=c(4,4,2,4))
panel = 1

data$duration_mo = round(data$duration_d/30.44)
plot(NA, NA, xlim=c(0,3), ylim=c(0,1.05), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
axis(side=1, at=0:3, line=0.5, labels=0:3, lwd=0, lwd.ticks=1)
axis(side=1, at=(0:6)/2, line=0, labels=NA, lwd=1, lwd.ticks=1)
axis(side=2, at=(0:4)/4, labels=percent(0:4/4), las=2)
abline(h=0:4/4, lwd=.25, col='#777777')
abline(v=0:6/2, lwd=.25, col='#777777')
for (mutation in top_hp_muts) {
  sfobj = survfit(Surv(duration_mo,death_surv_status)~1,data=subset(data, surv_age > 1 & data$family_mutation==mutation))
  points(c(0,sfobj$time/12), c(1,sfobj$surv), type = 'l', lwd=3, col=mutparms$col[mutparms$mut==mutation])
}
mtext(side=1, line=2.5, text='years from first symptom')
mtext(side=2, line=3, text='proportion surviving')
#mtext(side=3, line=0, cex=1.2, font=2, text='disease duration in genetic prion disease')
mtext(side=4, line=0.5, at=c(.75,.05,0), text=mutparms[top_hp_muts,'mut'], col=mutparms[top_hp_muts,'col'], font=2, las=2)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)
panel = panel + 1



data$duration_mo = round(data$duration_d/30.44)
plot(NA, NA, xlim=c(0,30), ylim=c(0,1.05), xaxs='i', yaxs='i', axes=FALSE, ann=FALSE)
axis(side=1, at=(0:3)*10, line=0.5, labels=0:3*10, lwd=0, lwd.ticks=1)
axis(side=1, at=(0:6)*5, line=0, labels=NA, lwd=1, lwd.ticks=1)
axis(side=2, at=(0:4)/4, labels=percent(0:4/4), las=2)
abline(h=0:4/4, lwd=.25, col='#777777')
abline(v=0:6*5, lwd=.25, col='#777777')
for (mutation in c(top_hp_muts,supp_hp_muts)) {
  sfobj = survfit(Surv(duration_mo,death_surv_status)~1,data=subset(data, surv_age > 1 & data$family_mutation==mutation))
  par(new=TRUE)
  #plot(sfobj, lwd=3, col=mutparms$col[mutparms$mut==mutation], ann=FALSE, axes=FALSE, conf.int=FALSE, xlim=c(0,30*12), ylim=c(0,1.05), xaxs='i', yaxs='i')
  points(c(0,sfobj$time/12), c(1,sfobj$surv), type = 'l', lwd=3, col=mutparms$col[mutparms$mut==mutation])
  #text(x=max(sfobj$time)/12, y=sfobj$surv[length(sfobj$surv)], pos=3, labels=mutation, col=mutparms$col[mutparms$mut==mutation], font=2)
}
mtext(side=1, line=2.5, text='years from first symptom')
mtext(side=2, line=3, text='proportion surviving')
#mtext(side=3, line=0, cex=1.2, font=2, text='disease duration in genetic prion disease')
#legend('topright',mutparms$mut,col=mutparms$col,text.col=mutparms$col,text.font=2)

# hard-code custom text locations - no better way to do this it seems
hardcoded_order = c('E200K','D178N','A117V','P102L','6-OPRI','P105L')
mtext(side=1, at=c(1,4,7,12.5,21,25), text=hardcoded_order, col=mutparms$col[match(hardcoded_order,mutparms$mut)], font=2)
text(x=20, y=.6, pos=4, labels='5-OPRI', col=mutparms$col[mutparms$mut=='5-OPRI'], font=2)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)

dev.off()

# explain why the 5-OPRI curve looks weird
data[data$family_mutation=='5-OPRI' & !is.na(data$duration_mo), c('duration_mo','death_surv_status')]
summary(survfit(Surv(duration_mo,death_surv_status)~1,data=subset(data, surv_age > 1 & data$family_mutation=='5-OPRI')))


data$duration_mo = round(data$duration_d/30.44)


ntests = 22
modifiers = data.frame(variable=character(ntests), 
                       mutation=character(ntests), 
                       comparison=character(ntests), 
                       n=character(ntests), 
                       test=character(ntests), 
                       p=character(ntests), 
                       bc=character(ntests),
                       literature=character(ntests))

# semi-automatedly, mostly manually, fill in this table
i = 1
# age of onset ~ codon 129 diplotype
for (mut in c('P102L','D178N','E200K')) {
  possibleError = tryCatch({
    diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mut & surv_age > 1))
    p_raw = surv_pval(diffobj)
    p_string = formatC(p_raw,format='g',digits=2)
    compared_groups = paste(gsub(', ','/',gsub("(cis|trans)129=","",names(diffobj$n))), collapse=' vs. ')
    n_vals = paste(diffobj$n, collapse=' vs. ')
    modifiers$variable[i] = 'onset'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = compared_groups 
    modifiers$n[i] = paste(diffobj$n, collapse=' vs. ')
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = p_string
    modifiers$literature[i] = ''
  }, warning = function(w) { w
  }, error = function (e) { e
  }, finally = { 
  })
  if (inherits(possibleError, 'error')) {
    modifiers$variable[i] = 'onset'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = 'diplotypes'
    #n_m = sum(data$cis129[data$family_mutation==mut & data$surv_age > 1 & !is.na(data$surv_age)]=='M',na.rm=T)
    modifiers$n[i] = 'error, fill in manually'
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = 'N/A'
    modifiers$literature[i] = ''
  }
  i = i + 1
  
  # age of onset ~ codon 129 genotype
  possibleError = tryCatch({
    diffobj = survdiff(Surv(surv_age,surv_status)~gt129,data=subset(data, family_mutation==mut & surv_age > 1))
    p_raw = surv_pval(diffobj)
    p_string = formatC(p_raw,format='g',digits=2)
    compared_groups = paste(gsub("gt129=","",names(diffobj$n)), collapse=' vs. ')
    n_vals = paste(diffobj$n, collapse=' vs. ')
    modifiers$variable[i] = 'onset'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = compared_groups 
    modifiers$n[i] = paste(diffobj$n, collapse=' vs. ')
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = p_string
    modifiers$literature[i] = ''
  }, warning = function(w) { w
  }, error = function (e) { e
  }, finally = { 
  })
  if (inherits(possibleError, 'error')) {
    modifiers$variable[i] = 'onset'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = 'genotypes'
    #n_m = sum(data$cis129[data$family_mutation==mut & data$surv_age > 1 & !is.na(data$surv_age)]=='M',na.rm=T)
    modifiers$n[i] = 'error, fill in manually'
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = 'N/A'
    modifiers$literature[i] = ''
  }
  i = i + 1
}




# duration ~ codon 129 diplotype
for (mut in c('P102L','D178N','E200K')) {
  possibleError = tryCatch({
    diffobj = survdiff(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mut & surv_age > 1))
    p_raw = surv_pval(diffobj)
    p_string = formatC(p_raw,format='g',digits=2, flag='#')
    compared_groups = paste(gsub(', ','/',gsub("(cis|trans)129=","",names(diffobj$n))), collapse=' vs. ')
    n_vals = paste(diffobj$n, collapse=' vs. ')
    modifiers$variable[i] = 'duration'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = compared_groups 
    modifiers$n[i] = paste(diffobj$n, collapse=' vs. ')
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = p_string
    modifiers$literature[i] = ''
  }, warning = function(w) { w
  }, error = function (e) { e
  }, finally = { 
  })
  if (inherits(possibleError, 'error')) {
    modifiers$variable[i] = 'duration'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = 'diplotypes'
    #n_m = sum(data$cis129[data$family_mutation==mut & data$surv_age > 1 & !is.na(data$surv_age)]=='M',na.rm=T)
    modifiers$n[i] = 'error, fill in manually'
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = 'N/A'
    modifiers$literature[i] = ''
  }
  i = i + 1
  
  # duration ~ codon 129 genotype
  possibleError = tryCatch({
    diffobj = survdiff(Surv(duration_mo,death_surv_status)~gt129, data=subset(data, family_mutation==mut & surv_age > 1))
    p_raw = surv_pval(diffobj)
    p_string = formatC(p_raw,format='g',digits=2, flag='#')
    compared_groups = paste(gsub("gt129=","",names(diffobj$n)), collapse=' vs. ')
    n_vals = paste(diffobj$n, collapse=' vs. ')
    modifiers$variable[i] = 'duration'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = compared_groups 
    modifiers$n[i] = paste(diffobj$n, collapse=' vs. ')
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = p_string
    modifiers$literature[i] = ''
  }, warning = function(w) { w
  }, error = function (e) { e
  }, finally = { 
  })
  if (inherits(possibleError, 'error')) {
    modifiers$variable[i] = 'duration'
    modifiers$mutation[i] = mut
    modifiers$comparison[i] = 'genotypes'
    #n_m = sum(data$cis129[data$family_mutation==mut & data$surv_age > 1 & !is.na(data$surv_age)]=='M',na.rm=T)
    modifiers$n[i] = 'error, fill in manually'
    modifiers$test[i] = 'log-rank'
    modifiers$p[i] = 'N/A'
    modifiers$literature[i] = ''
  }
  i = i + 1
}


# parent-child correlation
for (mut in c('P102L','D178N','E200K')) {
  pc = sqldf(paste("
                   select   parent.yob pyob, child.yob cyob, parent.ad pad, child.ad cad
                   from     data parent, data child
                   where    (parent.iid = child.father or parent.iid = child.mother)
                   and      child.yob is not null and parent.ad is not null and child.ad is not null
                   and      parent.family_mutation = '",mut,"'
                   ;",sep=''))
  m = lm(cad ~ pad + cyob, data=pc)
  p = summary(m)$coefficients["pad","Pr(>|t|)"]
  n = nrow(pc)
  modifiers$variable[i] = 'onset'
  modifiers$mutation[i] = mut
  modifiers$comparison[i] = 'parent vs. child'
  modifiers$n[i] = n
  modifiers$test[i] = 'linear regression'
  modifiers$p[i] = formatC(p, format='fg', digits=2, flag='#')
  modifiers$literature[i] = ''
  i = i + 1
}


# age of onset vs. sex, controlling for mutation
coxobj = coxph(Surv(surv_age,surv_status)~sex+family_mutation,data=subset(data, surv_age > 1 & sex %in% c('M','F') & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'onset'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'men vs. women'
n_m = sum(data$sex[!is.na(data$surv_age) & data$surv_age > 1]=='M', na.rm=T)
n_f = sum(data$sex[!is.na(data$surv_age) & data$surv_age > 1]=='F', na.rm=T)
modifiers$n[i] = paste(n_m, 'vs.', n_f)
modifiers$test[i] = 'Cox'
modifiers$p[i] = formatC(summary(coxobj)$coefficients['sexM','Pr(>|z|)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i + 1

# duration vs. sex, controlling for mutation
m = lm(duration_d ~ sex + family_mutation, data = subset(data, sex %in% c('M','F') & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'duration'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'men vs. women'
n_m = sum(data$sex[!is.na(data$duration_d) & data$surv_age > 1]=='M', na.rm=T)
n_f = sum(data$sex[!is.na(data$duration_d) & data$surv_age > 1]=='F', na.rm=T)
modifiers$n[i] = paste(n_m, 'vs.', n_f)
modifiers$test[i] = 'linear regression'
modifiers$p[i] = formatC(summary(m)$coefficients['sexM','Pr(>|t|)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i + 1

# onset vs. ascertainment, controlling for mutation
coxobj = coxph(Surv(surv_age,surv_status)~as.factor(ascertainment)+family_mutation,data=subset(data, surv_age > 1 & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'onset'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'direct vs. indirect ascertainment'
n_d = sum(data$ascertainment[!is.na(data$surv_age) & data$surv_age > 1]==1, na.rm=T)
n_i = sum(data$ascertainment[!is.na(data$surv_age) & data$surv_age > 1]==2, na.rm=T)
modifiers$n[i] = paste(n_d, 'vs.', n_i)
modifiers$test[i] = 'Cox'
modifiers$p[i] = formatC(summary(coxobj)$coefficients['as.factor(ascertainment)2','Pr(>|z|)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i + 1

# duration vs. ascertainment
m = lm(duration_d ~ as.factor(ascertainment) + family_mutation, data=subset(data, family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'duration'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'direct vs. indirect ascertainment'
n_d = sum(data$ascertainment[!is.na(data$duration_d)]==1, na.rm=T)
n_i = sum(data$ascertainment[!is.na(data$duration_d)]==2, na.rm=T)
modifiers$n[i] = paste(n_d, 'vs.', n_i)
modifiers$test[i] = 'linear regression'
modifiers$p[i] = formatC(summary(m)$coefficients['as.factor(ascertainment)2','Pr(>|t|)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i+1 


# onset vs. year of birth
m = lm(surv_age ~ family_mutation + yob, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'onset'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'year of birth'
n = length(summary(m)$residuals)
modifiers$n[i] = paste(n)
modifiers$test[i] = 'linear regression'
modifiers$p[i] = formatC(summary(m)$coefficients['yob','Pr(>|t|)'], format='g', digits=2, flag='#')
modifiers$literature[i] = ''
i = i+1 


# onset vs. study
anova = aov(surv_age ~ family_mutation + study, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'onset'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'study centers'
n_tot = nrow(subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
n_centers = length(unique(data$study))
modifiers$n[i] = paste(n_tot)
modifiers$test[i] = 'two-way ANOVA'
modifiers$p[i] = formatC(summary(anova)[[1]]['study','Pr(>F)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i+1 


# year of onset or death -- use the surv_age variable which is onset for slow mutations and coalesce(onset,death) for rapid mutations
data$yood = data$yob + data$surv_age
# onset vs. year of onset
m = lm(surv_age ~ family_mutation + study + yood, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
modifiers$variable[i] = 'onset'
modifiers$mutation[i] = 'top three'
modifiers$comparison[i] = 'year of onset'
n = length(summary(m)$residuals)
modifiers$n[i] = paste(n)
modifiers$test[i] = 'linear regression'
modifiers$p[i] = formatC(summary(m)$coefficients['yood','Pr(>|t|)'], format='fg', digits=2, flag='#')
modifiers$literature[i] = ''
i = i+1 




modifiers$bc = formatC(pmin(1.0,as.numeric(modifiers$p) * nrow(modifiers)), format='g', digits=2, flag='#')
modifiers$bc[modifiers$p=='N/A'] = 'N/A'

modifiers
write.table(modifiers,'figures/table_s7.tsv',sep='\t',row.names=F,col.names=T,quote=F)
WriteXLS(modifiers,'figures/table_s7.xls')


three_colors = c('#66c2a5','#fc8d62','#8da0cb')
four_colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c')


# ---- Figure S4. onset and codon 129
pdf('figures/figure_s4.pdf',width=8,height=8)
par(mfrow=c(3,2), mar=c(4,5,3,1))
panel = 1
for (mutation in top_hp_muts) {
  
  diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1))
  p_raw = surv_pval(diplotype_diffobj)
  p_string = formatC(p_raw,format='g',digits=2)
  
  data_subset = data$family_mutation==mutation & data$surv_age > 1 & !is.na(data$cis129) & !is.na(data$trans129) & !is.na(data$surv_age) & !is.na(data$family_mutation)
  groups = unique(paste(data$cis129[data_subset], data$trans129[data_subset], sep='/'))
  n = length(groups)
  plot(survfit(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1)),
       xlim=c(0,100),ylim=c(0,1.05),lwd=3,col=four_colors,
       yaxt='n',xaxs='i',yaxs='i',xlab='',ylab='')
  axis(side=2,at=(0:4)/4,labels=percent((0:4/4)),las=2)
  legend('topright',groups,col=four_colors,lwd=3,title='codon 129 (cis/trans)')
  mtext(side=1,line=2.0,text='age',cex=0.7)
  mtext(side=2,line=3.5,text='proportion alive and well',cex=0.7)
  mtext(side=3,line=0.5,cex=0.7,text=paste('P = ',p_string,', log-rank test',sep=''))
  mtext(side=3,line=1.5,cex=0.8,text=paste(mutation,'age of onset by codon 129 diplotype'),font=2)
  
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)
  panel = panel + 1
  
  genotype_diffobj = survdiff(Surv(surv_age,surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1))
  p_raw = surv_pval(genotype_diffobj)
  p_string = formatC(p_raw,format='g',digits=2)
  
  data_subset = data$family_mutation==mutation & data$surv_age > 1 & !is.na(data$gt129) & !is.na(data$surv_age) & !is.na(data$family_mutation)
  groups = unique(paste(data$gt129[data_subset], sep='/'))
  n = length(groups)
  plot(survfit(Surv(surv_age,surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1)),
       xlim=c(0,100),ylim=c(0,1.05),lwd=3,col=three_colors,
       yaxt='n',xaxs='i',yaxs='i',xlab='',ylab='')
  axis(side=2,at=(0:4)/4,labels=percent((0:4/4)),las=2)
  legend('topright',groups,col=three_colors,lwd=3,title='codon 129 (ignoring phase)')
  mtext(side=1,line=2.0,text='age',cex=0.7)
  mtext(side=2,line=3.5,text='proportion alive and well',cex=0.7)
  mtext(side=3,line=0.5,cex=0.7,text=paste('P = ',p_string,', log-rank test',sep=''))
  mtext(side=3,line=1.5,cex=0.8,text=paste(mutation,'age of onset by codon 129 genotype'),font=2)
  
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)
  panel = panel + 1
}
dev.off()






# ---- Figure S5. duration and codon 129
pdf('figures/figure_s5.pdf',width=8,height=8)
par(mfrow=c(3,2), mar=c(4,5,3,1))
panel = 1
for (mutation in top_hp_muts) {
  
  # for x-axis scale, round 90th percentile of duration up to the nearest year
  xmax = ceiling(quantile(data$duration_mo[data$family_mutation==mutation & data$surv_age > 1],.9,na.rm=TRUE)/12)*12
  
  diplotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1))
  p_raw = surv_pval(diplotype_diffobj)
  p_string = formatC(p_raw,format='g',digits=2)
  
  data_subset = data$family_mutation==mutation & data$surv_age > 1 & !is.na(data$cis129) & !is.na(data$trans129) & !is.na(data$surv_age) & !is.na(data$family_mutation)
  groups = unique(paste(data$cis129[data_subset], data$trans129[data_subset], sep='/'))
  n = length(groups)
  plot(survfit(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1)),
       xlim=c(0,xmax),ylim=c(0,1.05),lwd=3,col=four_colors,
       axes=FALSE,xaxs='i',yaxs='i',xlab='',ylab='')
  axis(side=1,at=seq(0,xmax,by=12))
  axis(side=1,at=seq(0,xmax,by=6),labels=NA,cex.axis=0.6)
  axis(side=2,at=(0:4)/4,labels=percent((0:4/4)),las=2)
  legend('topright',groups,col=four_colors,lwd=3,title='codon 129 (cis/trans)')
  mtext(side=1,line=2.0,text='months since first symptom',cex=0.7)
  mtext(side=2,line=3.5,text='proportion surviving',cex=0.7)
  mtext(side=3,line=0.5,cex=0.7,text=paste('P = ',p_string,', log-rank test',sep=''))
  mtext(side=3,line=1.5,cex=0.8,text=paste(mutation,'disease duration by codon 129 diplotype'),font=2)
  
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)
  panel = panel + 1
  
  genotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1))
  p_raw = surv_pval(genotype_diffobj)
  p_string = formatC(p_raw,format='g',digits=2)
  
  data_subset = data$family_mutation==mutation & data$surv_age > 1 & !is.na(data$gt129) & !is.na(data$surv_age) & !is.na(data$family_mutation)
  groups = unique(paste(data$gt129[data_subset], sep='/'))
  n = length(groups)
  plot(survfit(Surv(duration_mo,death_surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1)),
       xlim=c(0,xmax),ylim=c(0,1.05),lwd=3,col=three_colors,
       axes=FALSE,xaxs='i',yaxs='i',xlab='',ylab='')
  axis(side=1,at=seq(0,xmax,by=12))
  axis(side=1,at=seq(0,xmax,by=6),labels=NA,cex.axis=0.6)
  axis(side=2,at=(0:4)/4,labels=percent((0:4/4)),las=2)
  legend('topright',groups,col=three_colors,lwd=3,title='codon 129 (ignoring phase)')
  mtext(side=1,line=2.0,text='months since first symptom',cex=0.7)
  mtext(side=2,line=3.5,text='proportion surviving',cex=0.7)
  mtext(side=3,line=0.5,cex=0.7,text=paste('P = ',p_string,', log-rank test',sep=''))
  mtext(side=3,line=1.5,cex=0.8,text=paste(mutation,'disease duration by codon 129 genotype'),font=2)
  
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.3)
  panel = panel + 1
}
dev.off()



# ---- Additional stats for supplementary discussion


# onset in D178N
mutation = 'D178N'
diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','VV')))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MV','VV')))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV')))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV') & cis129=='M'))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(surv_age,surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV') & trans129=='M'))
diplotype_diffobj

# duration in E200K
mutation = 'E200K'
diplotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','VV')))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV')))
diplotype_diffobj
diplotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~cis129+trans129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV') & cis129=='M'))
diplotype_diffobj
genotype_diffobj = survdiff(Surv(duration_mo,death_surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV')))
genotype_diffobj

mutation = "P102L"
genotype_diffobj = survdiff(Surv(surv_age,surv_status)~gt129,data=subset(data, family_mutation==mutation & surv_age > 1 & gt129 %in% c('MM','MV')))
genotype_diffobj


# effect of diplotype and other factors on variance explained  - for supplementary discussion
m = lm(surv_age ~ family_mutation, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
summary(m)

m = lm(surv_age ~ family_mutation + cis129 + trans129, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
summary(m)
round(summary(m)$adj.r.squared,2)

m = lm(surv_age ~ family_mutation + study, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
summary(m)
round(summary(m)$adj.r.squared,2)

m = lm(surv_age ~ family_mutation + study + yood, data=subset(data, surv_status==1 & family_mutation %in% top_hp_muts))
summary(m)
round(summary(m)$adj.r.squared,2)

# some dissection of the natural history issue using Cox and left-truncation

# first, do our (very limited) prospective data accord with our retrospective data?
retrospective = data[is.na(data$ascertainment_age) & !is.na(data$surv_age),c('iid','family_mutation','surv_age','surv_status')]
retrospective$ascertainment_age = 0 # suppose that 'retrospective' is equivalent to ascertainment from birth
prospective = data[!is.na(data$ascertainment_age) & data$ascertainment_age < data$surv_age & !is.na(data$surv_age),c('iid','family_mutation','surv_age','surv_status','ascertainment_age')]

nrow(prospective)
sum(prospective$surv_age - prospective$ascertainment_age) / nrow(prospective)

retrospective$asc = 'retro'
prospective$asc = 'pro'
prore = rbind(retrospective, prospective)

coxph(Surv(time=ascertainment_age,time2=surv_age,event=surv_status,type='counting')~asc+family_mutation,data=prore)
# we see no indication that hazard is any different for people followed prospectively,
# however, we only have 24 individuals and only 5 events, so not much power


