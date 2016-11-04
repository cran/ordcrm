fiterations <-
function(dosehat, cohortsize=3, pseudotox, old.d, old.y,
	dosetox, truedosetoxmodeltype, design, targetDLT=0.30,
	pseudoweights=NA, discrete=FALSE, discretedoses=NA, numberdltrule=NA,
	lowerlimitrule=NA, upperlimitrule=NA, dltrule=NA,
	increaserule=NA, minimum=NA, maximum=NA,
	combine01=FALSE, iteration, rounddown=FALSE)
	{

#WARNINGS#
if ((discrete==FALSE) & rounddown==TRUE)	{
	stop('Rounding down only applied to discrete dose levels')
}
if ((discrete==TRUE) & is.na(discretedoses[1])==TRUE)	{
	stop('If using discrete doses you must specify dose levels to run this function')
}
if ((is.na(discretedoses[1])==FALSE) & discrete==FALSE)	{
	stop('If specifying discrete dose levels, discrete must equal TRUE')
}
if (truedosetoxmodeltype=='CRM' & (design=='POM' | design=='CR'))	{
	stop('If using an underlying model that is a 2 parameter binary logistic model (CRM), your design
	must also be a binary outcome CRM (design=CRM)')
}
if (truedosetoxmodeltype!='POM' & truedosetoxmodeltype!='CR' & truedosetoxmodeltype!='CRM')	{
	stop('You must specify if the underlying model is a POM (proportional odds model or violation),
	a CR (continuation ratio model), or a 2 parameter binary logistic model (CRM)')
}
if (design!='POM' & design!='CR' & design!='CRM')	{
	stop('You must specify the dose finding design you wish to use.
	Choices are POM (proportional odds model ordinal design), CRM (original binary CRM),
	or CR (continuation ratio model ordinal design)')
}

if (is.na(pseudoweights)==TRUE)	{
	pseudoweights<-cohortsize
}

lastdose<-old.d[length(old.d)]
d<-rep(dosehat, cohortsize)

if (truedosetoxmodeltype=="CRM")	{
	prob1<- dosetox[1] + dosetox[2]*dosehat
	prob1<-1/(1+exp(-prob1))
	prob0<-1-prob1

	phat<- prob1
	y<-stats::rbinom(cohortsize, 1, phat)
	numberdlty<-sum(y)
}

if (truedosetoxmodeltype=="POM")	{
	phat1ormore<- dosetox[1] + dosetox[2]*dosehat
	phat1ormore<-1/(1+exp(-phat1ormore))
	phat2ormore<- dosetox[3] + dosetox[4]*dosehat
	phat2ormore<-1/(1+exp(-phat2ormore))
	phat3ormore<- dosetox[5] + dosetox[6]*dosehat
	phat3ormore<-1/(1+exp(-phat3ormore))
	phat4ormore<- dosetox[7] + dosetox[8]*dosehat
	phat4ormore<-1/(1+exp(-phat4ormore))
	phat4<-phat4ormore
	phat3<-phat3ormore-phat4ormore
	phat2<-phat2ormore-phat3ormore
	phat1<-phat1ormore-phat2ormore
	phat0<-1-phat1ormore
	phat01<-phat0+phat1

	if (combine01==FALSE & design!='CRM')	{
		gen<-stats::rmultinom(n=1, size=cohortsize,
			prob=c(phat0, phat1, phat2, phat3, phat4))
		y <- c(rep(0,gen[1]), rep(1,gen[2]), rep(2, gen[3]), rep(3, gen[4]), rep(4, gen[5]))
		numberdlty<-sum(ifelse(y==3 | y==4,1,0))
	}
	if (combine01==TRUE & design!='CRM')	{
		gen<-stats::rmultinom(n=1, size=cohortsize,
			prob=c(phat01, phat2, phat3, phat4))
		y <- c(rep(0,gen[1]), rep(1,gen[2]), rep(2, gen[3]), rep(3, gen[4]))
		numberdlty<-sum(ifelse(y==2 | y==3,1,0))
	}
	if (design=='CRM')	{
		phat<- phat3ormore
		y<-stats::rbinom(cohortsize, 1, phat)
		numberdlty<-sum(y)}
}

if (truedosetoxmodeltype=="CR")	{
	crequal0<-1/(1+(exp(-(dosetox[1]+(dosetox[2]*dosehat)))))
	crequal1<-1/(1+(exp(-(dosetox[1]+dosetox[3]+(dosetox[4]*dosehat)))))
	crequal2<-1/(1+(exp(-(dosetox[1]+dosetox[5]+(dosetox[6]*dosehat)))))
	crequal3<-1/(1+(exp(-(dosetox[1]+dosetox[7]+(dosetox[8]*dosehat)))))
	phatmorethan0<-1-crequal0
	phatmorethan1<-(1-crequal1)*phatmorethan0
	phatmorethan2<-(1-crequal2)*phatmorethan1
	phatmorethan3<-(1-crequal3)*phatmorethan2
	phat4<-phatmorethan3
	phat3<-phatmorethan2-phatmorethan3
	phat2<-phatmorethan1-phatmorethan2
	phat1<-phatmorethan0-phatmorethan1
	phat0<-1-phatmorethan0
	phat01<-phat0+phat1

	if (combine01==FALSE & design!='CRM')	{
		gen<-stats::rmultinom(n=1, size=cohortsize,
		prob=c(phat0, phat1, phat2, phat3, phat4))
		y <- c(rep(0,gen[1]), rep(1,gen[2]), rep(2, gen[3]), rep(3, gen[4]), rep(4, gen[5]))
		numberdlty<-sum(ifelse(y==3 | y==4,1,0))
	}
	if (combine01==TRUE & design!='CRM')	{
		gen<-stats::rmultinom(n=1, size=cohortsize, prob=c(phat01, phat2, phat3, phat4))
		y <- c(rep(0,gen[1]), rep(1,gen[2]), rep(2, gen[3]), rep(3, gen[4]))
		numberdlty<-sum(ifelse(y==2 | y==3,1,0))
	}
	if (design=='CRM')	{
		phat<- phatmorethan2
		y<-stats::rbinom(cohortsize, 1, phat)
		numberdlty<-sum(y)
	}
}

yall<-append(old.y, y, after=length(old.y))
dall<-append(old.d, d, after=length(old.d))
lengthdata<-iteration*cohortsize
weight<-c(rep(pseudoweights/length(pseudotox), length(pseudotox)),
	rep(1, (iteration*cohortsize)))
weightpseudo<-round(((sum(weight[1:length(pseudotox)]))/(lengthdata+pseudoweights))*100,2)
weightdata<-round((sum(weight[(length(pseudotox)+1):length(weight)])/(lengthdata+pseudoweights))*100,2)

################
###BINARY CRM###
################
if (design=='CRM')	{
	reg<-suppressWarnings(rms::lrm(yall~dall, weights=weight))
	dose1<-round((log(targetDLT/(1-targetDLT))-reg$coef['Intercept'])/reg$coef['dall'])
	actualdltcohort<-numberdlty
	a<-constraints(dose1, lastdose, actualdltcohort, numberdltrule,
		lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum)
	rulesused<-a$rulesused
	dose<-a$dosehat.new
if (discrete==TRUE)	{
	discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
	esttoxdiscdose<-diffdose<-absdiffdose<-c()
	for (i in 1:length(discretedoses))	{
		esttoxdiscdose[i]<-1/(1+exp(-1*(reg$coef['Intercept']+ (reg$coef['dall']*discretedoses[i]))))
			diffdose[i]<-(dose-discretedoses[i])
			absdiffdose[i]<-abs(dose-discretedoses[i])
	}
	tt<-order(absdiffdose, decreasing=FALSE)
	ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
	if (rounddown==TRUE)	{
	discdose<-ifelse(ordernextdose[1,2]>dose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
	}
	if (rounddown==FALSE)	{
	discdose<-ordernextdose[1,2]
	}
}
}
########################
###POM ORDINAL DESIGN###
########################
if (design=='POM')	{
	reg<-suppressWarnings(rms::lrm(yall~dall, weights=weight))
	dose1<- ifelse(combine01==FALSE, round((log(targetDLT/(1-targetDLT)) - reg$coef[3])/reg$coef['dall']),
		round((log(targetDLT/(1-targetDLT)) - reg$coef[2])/reg$coef['dall']))
	actualdltcohort<-numberdlty
	a<-constraints(dose1, lastdose, actualdltcohort, numberdltrule,
		lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum)
	rulesused<-a$rulesused
	dose<-a$dosehat.new
if (discrete==TRUE)	{
	discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
	esttoxdiscdose<-diffdose<-absdiffdose<-c()
	for (i in 1:length(discretedoses))	{
	esttoxdiscdose[i]<-ifelse(combine01==FALSE, 1/(1+exp(-1*(reg$coef[3]+ (reg$coef['dall']*discretedoses[i])))),
		1/(1+exp(-1*(reg$coef[2]+ (reg$coef['dall']*discretedoses[i])))))
			diffdose[i]<-(dose-discretedoses[i])
			absdiffdose[i]<-abs(dose-discretedoses[i])
	}
	tt<-order(absdiffdose, decreasing=FALSE)
	ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
	if (rounddown==TRUE)	{
	discdose<-ifelse(ordernextdose[1,2]>dose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
	}
	if (rounddown==FALSE)	{
	discdose<-ordernextdose[1,2]
	}
}
}

#############################
###CR MODEL ORDINAL DESIGN###
#############################
if (design=='CR')	{
	y1<-rms::cr.setup(yall)
	newweights<-weight[y1$subs]
	dose<-dall[y1$subs]
	y<-y1$y
	cohort<-y1$cohort
	reg<-suppressWarnings(rms::lrm(y ~cohort+dose, weights=newweights))

	if (combine01==FALSE)	{
		d<-(-((1-targetDLT)/(targetDLT)))
		c<-(exp(reg$coef['Intercept'])+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1'])+
			exp(reg$coef['Intercept']+reg$coef['cohort=yall>=2']))
		b<-(exp((2*reg$coef['Intercept'])+reg$coef['cohort=yall>=1'])+exp((2*reg$coef['Intercept'])+reg$coef['cohort=yall>=2'])+
			exp((2*reg$coef['Intercept'])+reg$coef['cohort=yall>=1']+reg$coef['cohort=yall>=2']))
		a<-(exp((3*reg$coef['Intercept'])+reg$coef['cohort=yall>=1']+reg$coef['cohort=yall>=2']))
		m<-((2*(b^3))-(9*a*b*c)+(27*(a^2)*d))
		n2<-(m^2)-(4*(((b^2)-(3*a*c))^3))
		n<-sqrt(abs(n2))
		n<-ifelse(n2>=0, n, as.complex(n*(1i)))
		xx<-ifelse(n2>=0, ifelse((.5*(m+n))>0,(.5*(m+n))^(1/3),-(-(.5*(m+n)))^(1/3)),
		as.complex((.5*(m+n))^(1/3)))		# s #
		yy<-ifelse(n2>=0, ifelse((.5*(m-n))>0,(.5*(m-n))^(1/3),-(-(.5*(m-n)))^(1/3)),
		as.complex((.5*(m-n))^(1/3)))       # t #
		firstterm<-(-1*b)/(3*a)
		coef1<-as.complex((1+(1i*sqrt(3)))/(6*a))
		coef2<-as.complex((1-(1i*sqrt(3)))/(6*a))
		xi1<-as.complex(firstterm-(xx/(3*a))-(yy/(3*a))) # w #
		xi2<-as.complex(firstterm+(coef1*xx)+(coef2*yy)) # w #
		xi3<-as.complex(firstterm+(coef2*xx)+(coef1*yy)) # w #
		xi<-c(xi1, xi2, xi3)					#Gives the 3 roots#
		blah<-(log(xi)/reg$coef['dose'])
		blah<-suppressWarnings(as.double(blah))
		dcheck<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*blah)))*
			(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*blah)))*
			(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*blah))))
		possibledoses<-ifelse((dcheck>=targetDLT-0.01 & dcheck<=targetDLT+0.01),blah,NA)
		d.new<-blah[which(is.na(possibledoses)==FALSE)]
		dose1<-ifelse(length(d.new)==1, round(d.new,0), NA)
		actualdltcohort<-numberdlty
		a<-constraints(dose1, lastdose, actualdltcohort, numberdltrule,
			lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum)
		rulesused<-a$rulesused
		dose<-a$dosehat.new
if (discrete==TRUE)	{
	discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
	esttoxdiscdose<-diffdose<-absdiffdose<-c()
	for (i in 1:length(discretedoses))	{
	esttoxdiscdose[i]<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*discretedoses[i])))*
		(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*discretedoses[i])))*
		(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*discretedoses[i]))))
		diffdose[i]<-(dose-discretedoses[i])
		absdiffdose[i]<-abs(dose-discretedoses[i])
	}
	tt<-order(absdiffdose, decreasing=FALSE)
	ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
	if (rounddown==TRUE)	{
	discdose<-ifelse(ordernextdose[1,2]>dose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
	}
	if (rounddown==FALSE)	{
	discdose<-ordernextdose[1,2]
	}
}
	}

	if (combine01==TRUE)	{
		c<-(-((1-targetDLT)/(targetDLT)))
		b<-(exp(reg$coef['Intercept'])+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']))
		a<-(exp((2*reg$coef['Intercept'])+reg$coef['cohort=yall>=1']))
		z<-c(c,b,a)
		possibledoses<-polyroot(z)
		blah<-(log(possibledoses)/reg$coef['dose'])
		blah<-suppressWarnings(as.double(blah))
		dcheck<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*blah)))*
			(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*blah))))
		possibledoses<-ifelse((dcheck>=targetDLT-0.01 & dcheck<=targetDLT+0.01),blah,NA)
		d.new<-blah[which(is.na(possibledoses)==FALSE)]
		dose1<-(ifelse(length(d.new)==1, round(d.new,0), NA))
		actualdltcohort<-numberdlty
		a<-constraints(dose1, lastdose, actualdltcohort, numberdltrule,
			lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum)
		rulesused<-a$rulesused
		dose<-a$dosehat.new
if (discrete==TRUE)	{
	discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
	esttoxdiscdose<-diffdose<-absdiffdose<-c()
	for (i in 1:length(discretedoses))	{
	esttoxdiscdose[i]<-dcheck<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*discretedoses[i])))*
		(1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*discretedoses[i]))))
		diffdose[i]<-(dose-discretedoses[i])
		absdiffdose[i]<-abs(dose-discretedoses[i])
	}
	tt<-order(absdiffdose, decreasing=FALSE)
	ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
	if (rounddown==TRUE)	{
	discdose<-ifelse(ordernextdose[1,2]>dose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
	}
	if (rounddown==FALSE)	{
	discdose<-ordernextdose[1,2]
	}
}
	}
}
#Counting number of people at each dose#
if (discrete==TRUE)	{
dosestested<-dall[-c(1:length(pseudotox))]
numalreadytested<-c()
for (i in 1:length(discretedoses))	{
numalreadytested[i]<-length(which(dosestested==discretedoses[i]))}
}

if (discrete==FALSE)	{
dosestested<-NA
numalreadytested<-NA
esttoxdiscdose<-NA
	discdose<-NA}

	return(list(dose=dose, discdose=discdose, y=yall, d=dall,
		reg=reg, weightpseudo=weightpseudo, weightdata=weightdata,
		numberdlty=numberdlty, rulesused=rulesused,
		numalreadytested=numalreadytested,discretedoses=discretedoses))
}

