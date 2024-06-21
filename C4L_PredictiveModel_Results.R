### C4-L Final Results ###


library(ROCR)
library(ggplot2)
library(gridExtra)
library(survminer)
library(survival)

# Set up lists of discovery and validation samples
data_total = read.csv("SuppTable2.csv")

#67 discovery patients
data_D = data_total[1:67,]
#55 validation patients
data_V = data_total[68:122,]

data_TP1D=list()   
  
for(i in 1:ncol(data_D)) {       
  data_TP1D[[i]] <- data_D[ , i]     
} 
  
names(data_TP1D)=colnames(data_D)

data_TP1D= as_tibble(data_TP1D)


data_TP1V=list()   
  
for(i in 1:ncol(data_V)) {       
  data_TP1V[[i]] <- data_V[ , i]     
} 
  
names(data_TP1V)=colnames(data_V)

data_TP1V= as_tibble(data_TP1V)

# Univariate model results for UC-CaRE, CNA score, PSC for discovery data

log_reg_CNA <- glm(ProgressiontoAN~CNA_MSI_score, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_CNA), confint((log_reg_CNA)))),pval=coef(summary(log_reg_CNA))[,4])

log_reg_unresec <- glm(ProgressiontoAN~UCCaRETP1unresected, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_unresec), confint((log_reg_unresec)))),pval=coef(summary(log_reg_unresec))[,4])

log_reg_multif <- glm(ProgressiontoAN~UCCaRETP1multifocal, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_multif), confint((log_reg_multif)))),pval=coef(summary(log_reg_multif))[,4])

log_reg_inflam <- glm(ProgressiontoAN~UCCaRETP1inflammation, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_inflam), confint((log_reg_inflam)))),pval=coef(summary(log_reg_inflam))[,4])

log_reg_size <- glm(ProgressiontoAN~UCCaRETP1size, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_size), confint((log_reg_size)))),pval=coef(summary(log_reg_size))[,4])

log_reg_PSC <- glm(ProgressiontoAN~PSC, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_PSC), confint((log_reg_PSC)))),pval=coef(summary(log_reg_PSC))[,4])



# Calculate PPV and NPV of CNA score over time in validation data.

month6 = 365.25/2
month12 = 365.25
month24 = 2*365.25
month36 = 3*365.25
month48 = 4*365.25
month60 = 5*365.25

validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ CNA_MSI_score, data = data_TP1V)

validfit$strata_group = c(rep(1,validfit$strata[1]),rep(2,validfit$strata[2]))

# NPV for lowest group with 0 risk factors (strata_group ==1)

month6ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month6])
npvmonth6V = validfit$surv[validfit$strata_group==1][month6ind]
month12ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month12])
npvmonth12V = validfit$surv[validfit$strata_group==1][month12ind]
year2ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month24])
npv2V = validfit$surv[validfit$strata_group==1][year2ind]
year3ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month36])
npv3V = validfit$surv[validfit$strata_group==1][year3ind]
year4ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month48])
npv4V = validfit$surv[validfit$strata_group==1][year4ind]
year5ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month60])
npv5V = validfit$surv[validfit$strata_group==1][year5ind]

print(paste('npv 6 month VALIDATION = ',npvmonth6V ))
print(paste('npv 12 month VALIDATION = ',npvmonth12V ))
print(paste('npv year 2 VALIDATION = ',npv2V  ))
print(paste('npv year 3 VALIDATION = ',npv3V  ))
print(paste('npv year 4 VALIDATION = ',npv4V  ))
print(paste('npv year 5 VALIDATION = ',npv5V  ))

month6ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month6])
ppvmonth6V = 1-validfit$surv[validfit$strata_group==2][month6ind]
month12ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month12])
ppvmonth12V = 1-validfit$surv[validfit$strata_group==2][month12ind]
year2ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month24])
ppv2V = 1-validfit$surv[validfit$strata_group==2][year2ind]
year3ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month36])
ppv3V = 1-validfit$surv[validfit$strata_group==2][year3ind]
year4ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month48])
ppv4V = 1-validfit$surv[validfit$strata_group==2][year4ind]
year5ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month60])
ppv5V = 1-validfit$surv[validfit$strata_group==2][year5ind]

print(paste('ppv 6 month VALIDATION = ',ppvmonth6V ))
print(paste('ppv 12 month VALIDATION = ',ppvmonth12V ))
print(paste('ppv year 2 VALIDATION = ',ppv2V  ))
print(paste('ppv year 3 VALIDATION = ',ppv3V  ))
print(paste('ppv year 4 VALIDATION = ',ppv4V  ))
print(paste('ppv year 5 VALIDATION = ',ppv5V  ))


ppv_all=c(ppvmonth6V,ppvmonth12V,ppv2V ,ppv3V,ppv4V,ppv5V)
npv_all=c(npvmonth6V,npvmonth12V,npv2V ,npv3V,npv4V,npv5V)

### BOOTSTRAPPING by leave one out in validation data
v_total=length(data_TP1V$ProgressiontoAN)
ppv_boot = matrix(rep(0,6*1000),ncol=1000)
npv_boot = matrix(rep(0,6*1000),ncol=1000)
for (i in 1:1000){
	ind_to_remove = sample(x=1:v_total,size=v_total/4,prob=rep(1/v_total,v_total))
	validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ CNA_MSI_score, data = data_TP1V[-ind_to_remove,])
	validfit$strata_group = c(rep(1,validfit$strata[1]),rep(2,validfit$strata[2]))
	# NPV for low risk group (strata_group ==1)
	if (min(validfit$time[validfit$strata_group==1])>month6){
		npvmonth6V=1
	}
	else{
		month6ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month6])
		npvmonth6V = validfit$surv[validfit$strata_group==1][month6ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month12){
		npvmonth12V=1
	}
	else{
		month12ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month12])
		npvmonth12V = validfit$surv[validfit$strata_group==1][month12ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month24){
		npv2V=1
	}
	else {
		year2ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month24])
		npv2V = validfit$surv[validfit$strata_group==1][year2ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month36){
		npv3V=1
	}
	else {
		year3ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month36])
		npv3V = validfit$surv[validfit$strata_group==1][year3ind]	
	}
	if (min(validfit$time[validfit$strata_group==1])>month48){
		npv4V=1
	}
	else{
		year4ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month48])
		npv4V = validfit$surv[validfit$strata_group==1][year4ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month60){
		npv5V=1
	}
	else{
		year5ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month60])
		npv5V = validfit$surv[validfit$strata_group==1][year5ind]
	}

	npv_boot[,i] = c(npvmonth6V,npvmonth12V,npv2V ,npv3V,npv4V,npv5V)
	# PPV for high risk group (strata_group ==2)
	if (min(validfit$time[validfit$strata_group==2])>month6){
		ppvmonth6V=0
	}
	else{
		month6ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month6])
		ppvmonth6V = 1-validfit$surv[validfit$strata_group==2][month6ind]
	}
	if (min(validfit$time[validfit$strata_group==2])>month12){
		ppvmonth12V=0
	}
	else {
		month12ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month12])
		ppvmonth12V = 1-validfit$surv[validfit$strata_group==2][month12ind]
	}
	year2ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month24])
	ppv2V = 1-validfit$surv[validfit$strata_group==2][year2ind]
	year3ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month36])
	ppv3V = 1-validfit$surv[validfit$strata_group==2][year3ind]
	year4ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month48])
	ppv4V = 1-validfit$surv[validfit$strata_group==2][year4ind]
	year5ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month60])
	ppv5V = 1-validfit$surv[validfit$strata_group==2][year5ind]

	ppv_boot[,i] = c(ppvmonth6V,ppvmonth12V,ppv2V ,ppv3V,ppv4V,ppv5V)

}

lowerCI <- function(x){
	quantile(x,prob=0.025)
}

upperCI <- function(x){
	quantile(x,prob=0.975)
}

####################
### Figure 4B. #####
####################
## Plot NPV and PPV over time in Validation data with confidence intervals
dev.new()
plot(c(.5,1:5),ppv_all, ylim=c(0,1), type='o',col="orange",xaxt = "n",xlab="time since LGD diagnosis",cex=2,lwd=3,ylab="Predictive value in validation",pch=16,cex.axis=1.5,cex.lab=1.5)
lines(c(.5,1:5),npv_all, type='o', col='darkgreen',lwd=3,cex=2,pch=16)
axis(1,at=c(.5,1:5),labels=c('6 mos', '1 yr', '2 yr', '3 yr', '4 yr', '5 yr'), cex.lab=1.5, cex.axis=1.5)
legend("bottomright", legend = c("NPV", "PPV"), col = c("darkgreen", "orange"), lty = 1, cex = 2, bty='n',pch=16,lwd=3)
arrows(x0=c(.5,1:5), y0=apply(ppv_boot,1,FUN = lowerCI), x1=c(.5,1:5), y1=apply(ppv_boot,1,FUN = upperCI), code=3, angle=90, length=0.1, col="orange", lwd=2)
arrows(x0=c(.5,1:5), y0=apply(npv_boot,1,FUN = lowerCI), x1=c(.5,1:5), y1=apply(npv_boot,1,FUN = upperCI), code=3, angle=90, length=0.1, col="darkgreen", lwd=2)

# Output for paper
# ppv_all
#[1] 0.2727273 0.4545455 0.6363636 0.8181818 0.8181818 0.9090909
#npv_all
#[1] 1.0000000 0.9767442 0.9534884 0.9534884 0.9033048 0.9033048



## KM plots

# univariate CNA score
####################
### Figure 4C. #####
####################
dev.new()
discovfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ CNA_MSI_score, data = data_TP1D)
ggsurvplot(discovfit, data = data_TP1D,risk.table = TRUE ,pval=T,legend.labs=c('low-risk','high-risk'), ,pval.coord = c(4000,0.25), legend.title =  "Discovery (n=67)",  ylab = "HGD/CRC progression-free survival", xlab= "Follow-up (days)", palette='npg', size=2, censor.size=6, xlim=c(0,6000),break.x.by= 2000)

####################
### Figure 4D. #####
####################
dev.new()
validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ CNA_MSI_score, data = data_TP1V)
ggsurvplot(validfit, data = data_TP1V,risk.table = TRUE ,pval=T,legend.labs=c('low-risk','high-risk'), ,pval.coord = c(4000,0.25), legend.title =  "Validation (n=55)",  ylab = "HGD/CRC progression-free survival", xlab= "Follow-up (days)", palette='npg', size=2, censor.size=6)



## Multivariate model selection

log_reg_all <- glm(ProgressiontoAN~CNA_MSI_score+UCCaRETP1size+UCCaRETP1inflammation+UCCaRETP1multifocal+UCCaRETP1unresected, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_all), confint((log_reg_all)))),pval=coef(summary(log_reg_all))[,4])

step(log_reg_all)

# multivariate CNA + unresected selected

log_reg_CNA_unresec <- glm(ProgressiontoAN~CNA_MSI_score+UCCaRETP1unresected, data=data_TP1D, family="binomial")
cbind(exp(cbind(OR = coef(log_reg_CNA_unresec), confint((log_reg_CNA_unresec)))),pval=coef(summary(log_reg_CNA_unresec))[,4])



discov_data_glm_predictions = predict(log_reg_CNA_unresec,type='response')
temp= rep(0,length(discov_data_glm_predictions))
temp[discov_data_glm_predictions>0.5]=1

#########################
### same as plot 4C #####
#########################
dev.new()
discovfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ temp, data = data_TP1D)
ggsurvplot(discovfit, data = data_TP1D,risk.table = TRUE ,pval=T,legend.labs=c('low-risk','high-risk'), ,pval.coord = c(4000,0.25), legend.title =  "Discovery (n=67)",  ylab = "HGD/CRC progression-free survival", xlab= "Follow-up (days)", palette='npg', size=2, censor.size=6)



valid_data = data_TP1V[!is.na(data_TP1V$UCCaRETP1unresected),]

valid_data_glm_predictions = predict(log_reg_CNA_unresec,newdata=valid_data,type='response')
temp= rep(0,length(valid_data_glm_predictions))
temp[valid_data_glm_predictions>0.5]=1

####################
### Figure 4G. #####
####################
dev.new()
validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ temp, data = valid_data)
ggsurvplot(validfit, data = valid_data,risk.table = TRUE ,pval=T,legend.labs=c('low-risk','high-risk'), ,pval.coord = c(4000,0.25), legend.title =  "Validation (n=37)",  ylab = "HGD/CRC progression-free survival", xlab= "Follow-up (days)", palette='npg', size=2, censor.size=6)

## Calculate HR with 95% confidence intervals
mysurv=Surv(valid_data$CensorTime,valid_data$ProgressiontoAN) ~ temp
sdf=survdiff(mysurv)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

cbind("HR =" , HR, "95% CI:" , low95, up95)




validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ temp, data = valid_data)
validfit$strata_group = c(rep(1,validfit$strata[1]),rep(2,validfit$strata[2]))

# NPV for lowest group with 0 risk factors (strata_group ==1)

month6ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month6])
npvmonth6V = validfit$surv[validfit$strata_group==1][month6ind]
month12ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month12])
npvmonth12V = validfit$surv[validfit$strata_group==1][month12ind]
year2ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month24])
npv2V = validfit$surv[validfit$strata_group==1][year2ind]
year3ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month36])
npv3V = validfit$surv[validfit$strata_group==1][year3ind]
year4ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month48])
npv4V = validfit$surv[validfit$strata_group==1][year4ind]
year5ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month60])
npv5V = validfit$surv[validfit$strata_group==1][year5ind]

print(paste('npv 6 month VALIDATION = ',npvmonth6V ))
print(paste('npv 12 month VALIDATION = ',npvmonth12V ))
print(paste('npv year 2 VALIDATION = ',npv2V  ))
print(paste('npv year 3 VALIDATION = ',npv3V  ))
print(paste('npv year 4 VALIDATION = ',npv4V  ))
print(paste('npv year 5 VALIDATION = ',npv5V  ))

month6ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month6])
ppvmonth6V = 1-validfit$surv[validfit$strata_group==2][month6ind]
month12ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month12])
ppvmonth12V = 1-validfit$surv[validfit$strata_group==2][month12ind]
year2ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month24])
ppv2V = 1-validfit$surv[validfit$strata_group==2][year2ind]
year3ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month36])
ppv3V = 1-validfit$surv[validfit$strata_group==2][year3ind]
year4ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month48])
ppv4V = 1-validfit$surv[validfit$strata_group==2][year4ind]
year5ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month60])
ppv5V = 1-validfit$surv[validfit$strata_group==2][year5ind]

print(paste('ppv 6 month VALIDATION = ',ppvmonth6V ))
print(paste('ppv 12 month VALIDATION = ',ppvmonth12V ))
print(paste('ppv year 2 VALIDATION = ',ppv2V  ))
print(paste('ppv year 3 VALIDATION = ',ppv3V  ))
print(paste('ppv year 4 VALIDATION = ',ppv4V  ))
print(paste('ppv year 5 VALIDATION = ',ppv5V  ))


ppv_all=c(ppvmonth6V,ppvmonth12V,ppv2V ,ppv3V,ppv4V,ppv5V)
npv_all=c(npvmonth6V,npvmonth12V,npv2V ,npv3V,npv4V,npv5V)
### BOOTSTRAPPING by leave one out in validation data
v_total=length(valid_data$ProgressiontoAN)
ppv_boot = matrix(rep(0,6*1000),ncol=1000)
npv_boot = matrix(rep(0,6*1000),ncol=1000)
for (i in 1:1000){
	ind_to_remove = sample(x=1:v_total,size=v_total/4,prob=rep(1/v_total,v_total))
	#validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ CNA_MSI_score, data = valid_data[-ind_to_remove,])
	valid_data_glm_predictions = predict(log_reg_CNA_unresec,newdata=valid_data[-ind_to_remove,],type='response')
	temp= rep(0,length(valid_data_glm_predictions))
	temp[valid_data_glm_predictions>0.5]=1
	validfit = survfit(Surv(CensorTime,ProgressiontoAN) ~ temp, data = valid_data[-ind_to_remove,])

	validfit$strata_group = c(rep(1,validfit$strata[1]),rep(2,validfit$strata[2]))
	# NPV for low risk group (strata_group ==1)
	if (min(validfit$time[validfit$strata_group==1])>month6){
		npvmonth6V=1
	}
	else{
		month6ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month6])
		npvmonth6V = validfit$surv[validfit$strata_group==1][month6ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month12){
		npvmonth12V=1
	}
	else{
		month12ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month12])
		npvmonth12V = validfit$surv[validfit$strata_group==1][month12ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month24){
		npv2V=1
	}
	else {
		year2ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month24])
		npv2V = validfit$surv[validfit$strata_group==1][year2ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month36){
		npv3V=1
	}
	else {
		year3ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month36])
		npv3V = validfit$surv[validfit$strata_group==1][year3ind]	
	}
	if (min(validfit$time[validfit$strata_group==1])>month48){
		npv4V=1
	}
	else{
		year4ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month48])
		npv4V = validfit$surv[validfit$strata_group==1][year4ind]
	}
	if (min(validfit$time[validfit$strata_group==1])>month60){
		npv5V=1
	}
	else{
		year5ind = which.max(validfit$time[validfit$strata_group==1][validfit$time[validfit$strata_group==1]<=month60])
		npv5V = validfit$surv[validfit$strata_group==1][year5ind]
	}

	npv_boot[,i] = c(npvmonth6V,npvmonth12V,npv2V ,npv3V,npv4V,npv5V)
	# PPV for high risk group (strata_group ==2)
	if (min(validfit$time[validfit$strata_group==2])>month6){
		ppvmonth6V=0
	}
	else{
		month6ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month6])
		ppvmonth6V = 1-validfit$surv[validfit$strata_group==2][month6ind]
	}
	if (min(validfit$time[validfit$strata_group==2])>month12){
		ppvmonth12V=0
	}
	else {
		month12ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month12])
		ppvmonth12V = 1-validfit$surv[validfit$strata_group==2][month12ind]
	}
	if (min(validfit$time[validfit$strata_group==2])>month24){
		ppv2V=0
	}
	else {
		month24ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month24])
		ppv2V = 1-validfit$surv[validfit$strata_group==2][month24ind]
	}
	# year2ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month24])
	# ppv2V = 1-validfit$surv[validfit$strata_group==2][year2ind]
	year3ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month36])
	ppv3V = 1-validfit$surv[validfit$strata_group==2][year3ind]
	year4ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month48])
	ppv4V = 1-validfit$surv[validfit$strata_group==2][year4ind]
	year5ind = which.max(validfit$time[validfit$strata_group==2][validfit$time[validfit$strata_group==2]<=month60])
	ppv5V = 1-validfit$surv[validfit$strata_group==2][year5ind]

	ppv_boot[,i] = c(ppvmonth6V,ppvmonth12V,ppv2V ,ppv3V,ppv4V,ppv5V)

}

lowerCI <- function(x){
	quantile(x,prob=0.025)
}

upperCI <- function(x){
	quantile(x,prob=0.975)
}


####################
### Figure 4F. #####
####################


## Plot NPV and PPV over time in Validation data with confidence intervals
dev.new()
plot(c(.5,1:5),ppv_all, ylim=c(0,1), type='o',col="orange",xaxt = "n",xlab="time since LGD diagnosis",cex=2,lwd=3,ylab="Predictive value in validation",pch=16,cex.axis=1.5,cex.lab=1.5)
lines(c(.5,1:5),npv_all, type='o', col='darkgreen',lwd=3,cex=2,pch=16)
axis(1,at=c(.5,1:5),labels=c('6 mos', '1 yr', '2 yr', '3 yr', '4 yr', '5 yr'), cex.lab=1.5, cex.axis=1.5)
legend("bottomright", legend = c("NPV", "PPV"), col = c("darkgreen", "orange"), lty = 1, cex = 2, bty='n',pch=16,lwd=3)
arrows(x0=c(.5,1:5), y0=apply(ppv_boot,1,FUN = lowerCI), x1=c(.5,1:5), y1=apply(ppv_boot,1,FUN = upperCI), code=3, angle=90, length=0.1, col="orange", lwd=2)
arrows(x0=c(.5,1:5), y0=apply(npv_boot,1,FUN = lowerCI), x1=c(.5,1:5), y1=apply(npv_boot,1,FUN = upperCI), code=3, angle=90, length=0.1, col="darkgreen", lwd=2)


# Output for paper
# ppv_all
# [1] 0.4 0.8 0.8 1.0 1.0 1.0

# npv_all
# [1] 1.0000000 1.0000000 0.9677419 0.9677419 0.8960573 0.8960573



# Calculate AUCs


valid_data = data_TP1V
#### Remove validation patients who either progressed after 5 years or were censored before 5 years (13 patients) to compare with discovery data more closely ####

late_P = which(data_TP1V$CensorTime>1826 & data_TP1V$ProgressiontoAN==1)
early_NP = which(data_TP1V$CensorTime<1826 & data_TP1V$ProgressiontoAN==0)
	
data_TP1_valid = data_TP1V[-c(late_P, early_NP),]


valid_data = data_TP1_valid

# Discov
pred1 <- prediction(predict(log_reg_CNA), data_TP1D$ProgressiontoAN)
perf_AUC=performance(pred1,"auc")

print(perf_AUC@y.values[[1]])
#[1] 0.8535354

# Valid
pred1 <- prediction(predict(log_reg_CNA, newdata=valid_data), valid_data$ProgressiontoAN)
perf_AUC=performance(pred1,"auc")

print(perf_AUC@y.values[[1]])
#[1] 0.8571429


# CNA score + unresected only

#### Remove validation patients who do not have resection completion information####

valid_data = data_TP1_valid[!is.na(data_TP1_valid$UCCaRETP1unresected),]

# Discov
pred1 <- prediction(predict(log_reg_CNA_unresec), data_TP1D$ProgressiontoAN)
perf_AUC=performance(pred1,"auc")

print(perf_AUC@y.values[[1]])
#[1] 0.8974747

# Valid
pred1 <- prediction(predict(log_reg_CNA_unresec, newdata=valid_data), valid_data$ProgressiontoAN)
perf_AUC=performance(pred1,"auc")

print(perf_AUC@y.values[[1]])
#[1] 0.9464286
