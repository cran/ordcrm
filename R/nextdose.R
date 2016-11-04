#' Estimates the updated model and next dose to be tested for the Proportional Odds Model (POM),
#' Continuation Ratio (CR) Model CRM, or binary 2-parameter logistic model CRM.
#'
#' This function incorporates pseudodata created in \code{\link{pseudodata}} with the collected cohorts
#' of patient information and updates the fitted model selected (choice of POM, CR, or binary 2-parameter
#' logistic model). Using this newly updated model, it then re-estimates the next dose to be tested for
#' a specified target dose limiting toxicity (DLT) rate. It has the flexibility to add in various
#' safety constraints and can accommodate continuous or discrete doses. This function also has to option to
#' combine toxicity grades 0 and 1 as specified by CTCAEv4.0 into one category should clinical investigators
#' wish to do so.
#'
#' @param design Specifies which dose finding design you are running simulations on.
#' Choices are: POM, CR, or CRM.
#' @param pseudotox A vector of toxicity grades specified by the pseudodata.
#' @param pseudodose A vector of dose levels specified by the pseudodata.
#' @param collectedtox A vector of treated patient toxicity outcomes.
#' For j patients treated up to this point in the trial, collecteddose is in the form:
#' c(toxicity grade patient 1, ..., toxicity grade patient j). \cr
#' If combine01 is false, toxicity grades range from 0, 1, 2, 3, and 4.
#' If combine01 is true, toxicity grades are coded 0 (grades 0 and 1), 1 (grade 2), 2 (grade 3),
#' and 3 (grade 4). Toxicity grades are specified according to CTCAEv4.0.
#' @param collecteddose A vector of treated patient dose levels. For j patients treated up to this point
#' in the trial, collecteddose is in the form: c(dose patient 1, ..., dose patient j).
#' @param cohortsize Number of patients treated in a cohort. Defaults to 3.
#' @param targetDLT Target dose limiting toxicity (DLT) rate pre-specified by clinical investigators prior to the start of the trial.
#' Must be specified between 0 and 1. Defaults to 0.30.
#' @param pseudoweights Pseudoweights determine the amount of influence the pseudodata has on estimating the updated model.
#' Please refer to Van Meter et al (2010) for more information. We suggest setting the pseudoweights equal to the cohortsize.
#' When pseudoweights = X, the total pseudodata represent X individuals in the updated model. If not specified,
#' pseudoweights defaults to cohortsize.
#' Must be specified between 0 and 1. Defaults to 0.30.
#' @param discrete True/False. If discrete = TRUE this allows for discrete dose levels
#' to be specified prior to the start of the trial. Defaults to FALSE.
#' @param discretedoses Specified discrete dose levels if desired. They must be specified if discrete is equal to TRUE.
#' It is written for j dose levels as: c(d1, d2,..., dj).
#' @param rounddown True/False. Only applicable when using discrete dose levels.\cr
#' If rounddown = TRUE, the estimate dose from specified model will round down to the more conservative discrete dose level.\cr
#' If rounddown = FALSE, it will select the discrete dose closest to the estimated model selection. Defaults to FALSE.
#' @param numberdltrule Number of dose limiting toxicities (DLTs) that would be considered excessively unsafe if observed in a cohort of patients.
#' This would prompt the DLT constraint if there is one specified. Defaults to NA.
#' @param lowerlimitrule If lowerlimitrule is specified, this constraint will stop the trial if the dose
#' is estimated to be under a certain level as a safety precaution. Investigators may feel any dose under
#' this lowerlimitrule would suggest excessive toxicity with this investigational drug and no true MTD may exist.\cr
#' This can be specified as a dose level or percentage.\cr
#' If lowerlimitrule <1 (percentage), it takes the range of possible data (maximum - minimum) and makes the lower bound
#' when the trial would stop early as: minimum + (maximum - minimum)*lowerlimitrule.\cr
#' If lowerlimitrule >1 (dose level) then the lower bound when the trial would stop early is the lowerlimitrule.
#' @param upperlimitrule If upperlimitrule is specified, this constraint will stop the trial if the dose is estimated to be higher
#' than a certain level. \cr
#' If upperlimitrule <1 (percentage), it takes the range of possible data (maximum - minimum) and makes the upper bound
#' when the trial would stop early as: maximum - (maximum - minimum)*upperlimitrule. \cr
#' If upperlimitrule >1 (dose level) then the upper bound
#' when the trial would stop early is the upperlimitrule.
#' @param dltrule If the numberdltrule (i.e. 2 DLTs) occurs in the last cohort of patients tested,
#' the next estimated dose must decrease by an amount specified by the dltrule.
#' Defaults to NA.  \cr
#' If 0 < dltrule < 1 (percentage), with excessive DLTs in the last cohort,
#' the next estimated dose must decrease by the dltrule (percentage) of
#' the last tested dose level, i.e. (1-dltrule)*lastdose.\cr
#' If dltrule >= 1 (dose level), with excessive DLTs in the last cohort
#' the next estimated dose must decrease by dltrule dose level, i.e. lastdose-dltrule.
#' @param increaserule If increaserule is specified, then the next estimated dose
#' can only increase by an amount specified by the increaserule between tested cohorts. Defaults to NA.\cr
#' If 0 < increaserule < 1 (percentage), the next dose can only increase to a maximum of lastdose*(1+increaserule).\cr
#' If increaserule >= 1 (dose level), the next dose can only increase to a maximum of lastdose + increaserule.
#' @param minimum Minimum dose that will be considered to test patients. Must be specified if lowerlimitrule<1.
#' @param maximum Maximum dose that will be considered to test patients. Must be specified if upperlimitrule<1.
#' @param combine01 True/False. If combine01 = TRUE, toxicity grades 0 and 1 are combined into 1 category.
#' Therefore all toxicities must be coded: 0 (grades 0 and 1), 1 (grade2), 2 (grade 3),
#' and 3 (grade 4) according to CTCAEv4.0. Defaults to FALSE.
#' @param plotit True/False. If true, returns a plot of the updated model with the next dose to test in a cohort identified
#' as well as the current accrued patient toxicity outcomes. Defaults to TRUE.
#'
#' @details
#' If using a POM CRM, this function assumes a proportional odds model as described by Harrell (Harrell, 2001).
#' For combine01=FALSE, y equal to toxicity grade outcomes j in c(0, 1, 2, 3, 4) as specified by CTCAEv4.0, and x equal to the dose,
#' this is written in the form:\cr
#' \deqn{P(y>=j|x)=1/(1+exp(-(\alpha_j + \beta*x))),for j = 1, 2, 3, 4}\cr
#' For combine01=Ture, toxicity grades are now 0/1, 2, 3, and 4, and y is recoded as 0 = grades 0 and 1, 1 = grade 2,
#' 2 = grade 3, and 3 = grade 4.\cr
#'
#' If using a CR model design, this function assumes a continuation ratio model as described by Harrell (Harrell, 2001).
#' For combine01=FALSE, y equal to toxicity grade outcomes j in c(0, 1, 2, 3, 4) as specified by CTCAEv4.0, and x equal to the dose,
#' this is written in the form:\cr
#' \deqn{P(y=j|y>=j,x)=1/(1+exp(-(\alpha + \theta_j + \gamma*x))), j = 0, 1, 2, 3}\cr
#' For combine01=TURE, toxicity grades are now 0/1, 2, 3, and 4, and y is recoded as 0 = grades 0 and 1, 1 = grade 2,
#' 2 = grade 3, and 3 = grade 4.\cr
#'
#' If using a binary CRM, this assumes a standard 2-parameter logistic model.
#' For y equal to toxicity outcome (1 is a DLT, 0 is not a DLT) and x equal to the dose, it is written in the form:\cr
#' \deqn{P(y=1|x)=1/(1+exp(-(\alpha + \beta*x)))}\cr
#'
#' Estimated next dose is calculated by Pr(3 or 4 toxicity) <= targetDLT for the model utilized.
#'
#' @return \item{Next Dose with No Used Constraints}{Returns estimated dose with no safety constraints
#' implemented for a specified target DLT rate.}
#' @return \item{Next Dose with Constraints if Applicable}{Returns estimated dose with safety constraints used if
#' applicable for a specified target DLT rate.}
#' @return \item{Next Dose if Using Discrete Dose Levels}{Returns the next dose estimated if
#' discrete dose levels are requested.}
#' @return \item{Estimated Probability of a Grade 0 or 1 Toxicity}{Using the newly re-estimated model
#' based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 0 or 1 toxicity as specified by CTCAEv4.0.
#' Only applicable if combine01 = TRUE.}
#' @return \item{Estimated Probability of a Grade 0 Toxicity}{Using the newly re-estimated model based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 0 toxicity as specified by CTCAEv4.0.
#' Only applicable if combine01 = FALSE.}
#' @return \item{Estimated Probability of a Grade 1 Toxicity}{Using the newly re-estimated model based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 1 toxicity as specified by CTCAEv4.0.
#' Only applicable if combine01 = FALSE.}
#' @return \item{Estimated Probability of a Grade 2 Toxicity}{Using the newly re-estimated model based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 2 toxicity as specified by CTCAEv4.0.}
#' @return \item{Estimated Probability of a Grade 3 Toxicity}{Using the newly re-estimated model based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 3 toxicity as specified by CTCAEv4.0.}
#' @return \item{Estimated Probability of a Grade 4 Toxicity}{Using the newly re-estimated model based on the combination of pseudodata and updated patient information,
#' this returns the probability of experiencing a grade 4 toxicity as specified by CTCAEv4.0.}
#' @return \item{Estimated Regression Model}{Returns parameter estimates for the newly re-estimated model based on the combination of pseudodata and updated patient information.}
#' @return \item{Influence from Pseudodata}{Returns the percentage of influence the pseudodata has on estimating the new model.
#' Please refer to Van Meter et al (2010) for more information.}
#' @return \item{Influence from Collected Data}{Returns the percentage of influence the treated patient data has on estimating the new model.
#' Please refer to Van Meter et al (2010) for more information.}
#' @return \item{Constraints Used}{Identifies which safety constraint was implemented if any. Returns: DLT Rule Used, Increase Rule Used,
#'  Lower Bound Stop Rule Used, Upper Bound Stop Rule Used, or None Used.}
#'
#' @author Emily V. Dressler, PhD \cr
#' Markey Cancer Center\cr
#' Division of Cancer Biostatistics\cr
#' University of Kentucky\cr
#' \email{EmilyVDressler@@gmail.com}
#' @references 1. Van Meter EM, Garrett-Mayer E, Bandyopadhyay D. Proportional odds model for dose finding clinical trial designs with ordinal toxicity grading. Statistics in Medicine 2011; 30: 2070-2080. \cr
#'2. Van Meter EM, Garrett-Mayer E, Bandyopadhyay. Dose finding clinical trial design for ordinal toxicity grades using the continuation ratio model: an extension of the continual reassessment method. Clinical Trials 2012; 9(3): 303-313. \cr
#'3. Garrett-Mayer E. The continual reassessment method for dose-finding studies: a tutorial. Clinical Trials 2006; 3: 57-71. \cr
#'4. Piantadosi S, Fisher JD, Grossman S. Practical implementation of a modified continual reassessment method for dose-finding trials. Cancer Chemother Pharmacol 1998; 41: 429-436. \cr
#'5. Harrell FE, Jr. Regression Modeling Strategies with Application to Linear Models, Logistic Regression, and Survival Analysis. Springer: New York, NY, 2001. \cr
#'6. McCullagh P. Regression Models for Ordinal Data. Journal of the Royal Statistical Society. Series B (Methodological). 1980; 42: 109-142. \cr
#'7. CTCAE. Cancer Therapy Evaluation Program, Common Terminology Criteria for Adverse Events, Version 4.0, DCTD, NCI, NIH, DHHS (http://ctep/cancer.gov). In Cancer Therapy Evaluation Program, 2010.
#' @export
#' @seealso \code{\link{pseudodata}}
#' @examples
#' #Pseudodata toxicity grades#
#' initial.cr.y1 <- c(rep(0, 45), rep(1, 36), rep(2, 9), rep(3, 8), rep(4, 2),
#'                    rep(0, 24), rep(1, 31), rep(2, 15), rep(3, 26), rep(4, 4),
#'                    rep(0, 14), rep(1, 23), rep(2, 13), rep(3, 40), rep(4, 10),
#'                    rep(0, 1), rep(1, 4), rep(2, 5), rep(3, 35), rep(4, 55))
#'
#' #Pseudodata dose levels#
#' initial.cr.d1 <- c(rep(200, 100), rep(934, 100), rep(1467, 100), rep(3000, 100))
#'
#' #Pseudodata toxicity grades if combining grades 0 and 1 into 1 category#
#' combine.cr.y1 <- c(rep(0, 81), rep(1, 9), rep(2, 8), rep(3, 2),
#'                    rep(0, 55), rep(1, 15), rep(2, 26), rep(3, 4),
#'                    rep(0, 37), rep(1, 13), rep(2, 40), rep(3, 10),
#'                    rep(0, 5), rep(1, 5), rep(2, 35), rep(3, 55))
#'
#' #6 patients already treated at doses 1060 and 800 respectively
#' #CR model assumption
#'
#' #Example 1
#' nextdose(design='CR', pseudotox = initial.cr.y1, pseudodose = initial.cr.d1, cohortsize =,
#'          collectedtox = c(1, 4, 2, 0, 0, 1), collecteddose = c(1060, 1060, 1060, 800, 800, 800),
#'          targetDLT = 0.30, pseudoweights = 3, discrete = FALSE, discretedoses = NA,
#'          combine01 = FALSE, lowerlimitrule = 500)
#'
#' #Example 2
#' #Discrete doses and combining grades 0 and 1
#' nextdose(design='CR', pseudotox = combine.cr.y1, pseudodose = initial.cr.d1, cohortsize =,
#'         collectedtox = c(1, 3, 2, 0, 0, 1), collecteddose = c(1060, 1060, 1060, 800, 800, 800),
#'         targetDLT = 0.30, pseudoweights = 3, discrete = TRUE,
#'         discretedoses = c(200, 500, 100, 1200, 1500, 1800), combine01 = TRUE,
#'         lowerlimitrule = 500)







nextdose <-
function(design, pseudotox, pseudodose, collectedtox, collecteddose,
	cohortsize=3, targetDLT=0.30, pseudoweights=NA, discrete=FALSE, discretedoses=NA,
	rounddown=FALSE,
	numberdltrule=NA, lowerlimitrule=NA, upperlimitrule=NA, dltrule=NA,
	increaserule=NA,  minimum=NA, maximum=NA, combine01=FALSE, plotit=TRUE)
	{

#WARNINGS#
if ((design!='POM' & design!='CR' & design!='CRM'))	{
	stop('You must specify the dose finding design you wish to use.
	Choices are POM (proportional odds model ordinal design), CRM (original binary CRM),
	or CR (continuation ratio model ordinal design)')
}
if ((discrete==FALSE) & rounddown==TRUE)	{
	stop('Rounding down only applied to discrete dose levels')
}
if ((discrete==TRUE) & is.na(discretedoses[1])==TRUE)	{
	stop('If using discrete doses you must specify dose levels to run this function')
}
if ((is.na(discretedoses[1])==FALSE) & discrete==FALSE)	{
	stop('If specifying discrete dose levels, discrete must equal TRUE')
}
if (is.na(pseudoweights))	{
	pseudoweights<-cohortsize
}

###################################################
###Proportional Odds Model Design##################
###################################################
if (design=='POM')	{
lastdose<-collecteddose[length(collecteddose)]
yall<-append(pseudotox, collectedtox, after=length(pseudotox))
dall<-append(pseudodose, collecteddose, after=length(pseudodose))
weight<-c(rep(pseudoweights/length(pseudotox), length(pseudotox)),
	rep(1, (length(collectedtox))))
weightpseudo<-round(((sum(weight[1:length(pseudotox)]))/(length(collectedtox)+pseudoweights))*100,2)
weightdata<-round((sum(weight[(length(pseudotox)+1):length(weight)])/(length(collectedtox)+pseudoweights))*100,2)
reg<-suppressWarnings(rms::lrm(yall~dall, weights=weight))

if (combine01==FALSE)	{
		dose<- round((log(targetDLT/(1-targetDLT)) - reg$coef[3])/reg$coef['dall'])
		dosenoconstraint<-dose
		dltlastcohort<-collectedtox[((length(collectedtox)-cohortsize)+1):length(collectedtox)]
		dltlastcohort<-ifelse(dltlastcohort==3 | dltlastcohort==4, 1, 0)
		actualdltcohort<-sum(dltlastcohort)
		a<-constraints(dose, lastdose, actualdltcohort, numberdltrule,
			lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum)
		dose<-a$dosehat.new
		constraintused<-a$constraintused

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
		phat1ormore<- reg$coef[1]+ (reg$coef['dall']*discdose)
		phat1ormore<-1/(1+exp(-phat1ormore))
		phat2ormore<- reg$coef[2]+ (reg$coef['dall']*discdose)
		phat2ormore<-1/(1+exp(-phat2ormore))
		phat3ormore<- reg$coef[3]+ (reg$coef['dall']*discdose)
		phat3ormore<-1/(1+exp(-phat3ormore))
		phat4ormore<- reg$coef[4]+ (reg$coef['dall']*discdose)
		phat4ormore<-1/(1+exp(-phat4ormore))
}

if (discrete==FALSE)	{
	esttoxdiscdose<-NA
	discdose<-NA
		phat1ormore<- reg$coef[1]+ (reg$coef['dall']*dose)
		phat1ormore<-1/(1+exp(-phat1ormore))
		phat2ormore<- reg$coef[2]+ (reg$coef['dall']*dose)
		phat2ormore<-1/(1+exp(-phat2ormore))
		phat3ormore<- reg$coef[3]+ (reg$coef['dall']*dose)
		phat3ormore<-1/(1+exp(-phat3ormore))
		phat4ormore<- reg$coef[4]+ (reg$coef['dall']*dose)
		phat4ormore<-1/(1+exp(-phat4ormore))
	}

	phat4<-round(phat4ormore,4)*100
	phat3<-round(phat3ormore-phat4ormore,4)*100
	phat2<-round(phat2ormore-phat3ormore,4)*100
	phat1<-round(phat1ormore-phat2ormore,4)*100
	phat0<-round(1-phat1ormore,4)*100
	phat0or1<-NA
	probDLT<-phat3+phat4
	probnonDLT<-ifelse(combine01==TRUE, phat0or1+phat2, phat0+phat1+phat2)
}

if (combine01==TRUE)	{
		dose<- round((log(targetDLT/(1-targetDLT)) - reg$coef[2])/reg$coef['dall'])
		dosenoconstraint<-dose
		dltlastcohort<-collectedtox[((length(collectedtox)-cohortsize)+1):length(collectedtox)]
		dltlastcohort<-ifelse(dltlastcohort==2 | dltlastcohort==3, 1, 0)
		actualdltcohort<-sum(dltlastcohort)
		a<-constraints(dose, lastdose, actualdltcohort, numberdltrule,
			lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum)
		dose<-a$dosehat.new
		constraintused<-a$constraintused
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
		phat2ormore<- reg$coef[1]+ (reg$coef['dall']*discdose)
		phat2ormore<-1/(1+exp(-phat2ormore))
		phat3ormore<- reg$coef[2]+ (reg$coef['dall']*discdose)
		phat3ormore<-1/(1+exp(-phat3ormore))
		phat4ormore<- reg$coef[3]+ (reg$coef['dall']*discdose)
		phat4ormore<-1/(1+exp(-phat4ormore))
}

if (discrete==FALSE)	{
		esttoxdiscdose<-NA
		discdose<-NA
		phat2ormore<- reg$coef[1]+ (reg$coef['dall']*dose)
		phat2ormore<-1/(1+exp(-phat2ormore))
		phat3ormore<- reg$coef[2]+ (reg$coef['dall']*dose)
		phat3ormore<-1/(1+exp(-phat3ormore))
		phat4ormore<- reg$coef[3]+ (reg$coef['dall']*dose)
		phat4ormore<-1/(1+exp(-phat4ormore))
}

	phat4<-round(phat4ormore,4)*100
	phat3<-round(phat3ormore-phat4ormore,4)*100
	phat2<-round(phat2ormore-phat3ormore,4)*100
	phat0or1<-round(1-phat2ormore,4)*100
	phat0<-NA
	phat1<-NA
	probDLT<-phat3+phat4
	probnonDLT<-ifelse(combine01==TRUE, phat0or1+phat2, phat0+phat1+phat2)
}

if (plotit==TRUE)	{
if (discrete==FALSE)	{
	if (combine01==FALSE)	{
		graphics::plot(c(0,(max(dall)*1.05)),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey85",colors)
			colors<-ifelse(collectedtox==2,"grey20",colors)
			colors<-ifelse(collectedtox==3,"grey60",colors)
			colors<-ifelse(collectedtox==4,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			type<-ifelse(collectedtox==4,19,type)
			area<-ifelse(collectedtox==3 | collectedtox==4, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp1<- reg$coef[1]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		tmp2<- reg$coef[2]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		tmp3<- reg$coef[3]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		tmp4<- reg$coef[4]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp1)),lty=1, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		graphics::points(dose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x=(dose*1.05), y=(targetDLT+.03), labels=dose)
		graphics::legend(0,.95,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Next Dose"),
			lty=c(1,2,3,4,-1),pch=c(-1,-1,-1,-1,15),lwd=c(2,2,2,2,1),
			col=c(1,1,1,1,"red"),bg="white")
		graphics::legend("bottom",c("Grade 0","Grade 1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey82","grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Proportional Odds Model CRM")
		}
	if (combine01==TRUE)	{
		graphics::plot(c(0,(max(dall)*1.05)),c(-0.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey60",colors)
			colors<-ifelse(collectedtox==2,"grey85",colors)
			colors<-ifelse(collectedtox==3,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			area<-ifelse(collectedtox==2 | collectedtox==3, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp2<- reg$coef[1]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		tmp3<- reg$coef[2]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		tmp4<- reg$coef[3]+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		graphics::points(dose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x=(dose*1.05), y=(targetDLT+.03), labels=dose)
		graphics::legend(0,.95,c("p(y>=2)","p(y>=3)","p(y>=4)","Next Dose"),
			lty=c(2,3,4,-1),pch=c(-1,-1,-1,15),lwd=c(2,2,2,1),
			col=c(1,1,1,"red"),bg="white")
		graphics::legend("bottom",c("Grade 0/1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Proportional Odds Model CRM")
		}
	}
if (discrete==TRUE)	{
	if (combine01==FALSE)	{
		lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
		graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0)))),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey20",colors)
			colors<-ifelse(collectedtox==2,"grey60",colors)
			colors<-ifelse(collectedtox==3,"grey85",colors)
			colors<-ifelse(collectedtox==4,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			type<-ifelse(collectedtox==4,19,type)
			area<-ifelse(collectedtox==3 | collectedtox==4, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp1<- reg$coef[1]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp2<- reg$coef[2]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp3<- reg$coef[3]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp4<- reg$coef[4]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp1)),lty=1, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		for (i in 1:length(discretedoses))	{
		graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
		graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
		graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
		graphics::legend(0,.95,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Next Discrete Dose","Estimated Toxicity"),
			lty=c(1,2,3,4,-1,-1),pch=c(-1,-1,-1,-1,15,15),lwd=c(2,2,2,2,1,1),
			col=c(1,1,1,1,"red","gray90"),bg="white")
		graphics::legend("bottom",c("Grade 0","Grade 1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey82","grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Proportional Odds Model CRM")

		}
	if (combine01==TRUE)	{
		lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
		graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0)))),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey60",colors)
			colors<-ifelse(collectedtox==2,"grey85",colors)
			colors<-ifelse(collectedtox==3,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			area<-ifelse(collectedtox==2 | collectedtox==3, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp2<- reg$coef[1]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp3<- reg$coef[2]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp4<- reg$coef[3]+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		for (i in 1:length(discretedoses))	{
		graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
		graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
		graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
		graphics::legend(0,.95,c("p(y>=2)","p(y>=3)","p(y>=4)","Next Discrete Dose","Estimated Toxicity"),
			lty=c(2,3,4,-1,-1),pch=c(-1,-1,-1,15,15),lwd=c(2,2,2,1,1),
			col=c(1,1,1,"red","gray90"),bg="white")
		graphics::legend("bottom",c("Grade 0/1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Proportional Odds Model CRM")
		}
	}
}
}

######################################################
###Continuation Ratio Design##########################
######################################################
if (design=='CR')	{

lastdose<-collecteddose[length(collecteddose)]
yall<-append(pseudotox, collectedtox, after=length(pseudotox))
dall<-append(pseudodose, collecteddose, after=length(pseudodose))
weight<-c(rep((pseudoweights)/length(pseudotox), length(pseudotox)),
	rep(1, (length(collectedtox))))
weightpseudo<-round(((sum(weight[1:length(pseudotox)]))/(length(collectedtox)+pseudoweights))*100,2)
weightdata<-round((sum(weight[(length(pseudotox)+1):length(weight)])/(length(collectedtox)+pseudoweights))*100,2)
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
	dose<-ifelse(length(d.new)==1, round(d.new,0), NA)
	dosenoconstraint<-dose
	dltlastcohort<-collectedtox[((length(collectedtox)-cohortsize)+1):length(collectedtox)]
	dltlastcohort<-ifelse(dltlastcohort==3 | dltlastcohort==4, 1, 0)
	actualdltcohort<-sum(dltlastcohort)
	a<-constraints(dose, lastdose, actualdltcohort, numberdltrule,
		lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum)
	constraintused<-a$constraintused
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
		crequal0<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*discdose)))))
		crequal1<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*discdose)))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*discdose)))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*discdose)))))
}
if (discrete==FALSE)	{
	esttoxdiscdose<-NA
	discdose<-NA
		crequal0<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*dose)))))
		crequal1<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*dose)))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*dose)))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*dose)))))
}
	phatmorethan0<-1-crequal0
	phatmorethan1<-(1-crequal1)*phatmorethan0
	phatmorethan2<-(1-crequal2)*phatmorethan1
	phatmorethan3<-(1-crequal3)*phatmorethan2
	phat4<-round(phatmorethan3,4)*100
	phat3<-round(phatmorethan2-phatmorethan3,4)*100
	phat2<-round(phatmorethan1-phatmorethan2,4)*100
	phat1<-round(phatmorethan0-phatmorethan1,4)*100
	phat0<-round(1-phatmorethan0,4)*100
	phat0or1<-NA
	probDLT<-phat3+phat4
	probnonDLT<-ifelse(combine01==TRUE, phat0or1+phat2, phat0+phat1+phat2)
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
	dose<-(ifelse(length(d.new)==1, round(d.new,0), NA))
	dosenoconstraint<-dose
	dltlastcohort<-collectedtox[((length(collectedtox)-cohortsize)+1):length(collectedtox)]
	dltlastcohort<-ifelse(dltlastcohort==2 | dltlastcohort==3, 1, 0)
	actualdltcohort<-sum(dltlastcohort)
	a<-constraints(dose, lastdose, actualdltcohort, numberdltrule,
		lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum)
	constraintused<-a$constraintused
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
		crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*discdose)))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*discdose)))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*discdose)))))
}
if (discrete==FALSE)	{
	esttoxdiscdose<-NA
	discdose<-NA
		crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*dose)))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*dose)))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*dose)))))
}
	phatmorethan01<-1-crequal01
	phatmorethan2<-(1-crequal2)*phatmorethan01
	phatmorethan3<-(1-crequal3)*phatmorethan2
	phat4<-round(phatmorethan3,4)*100
	phat3<-round(phatmorethan2-phatmorethan3,4)*100
	phat2<-round(phatmorethan01-phatmorethan2,4)*100
	phat0or1<-round(1-phatmorethan01,4)*100
	phat0<-NA
	phat1<-NA
	probDLT<-phat3+phat4
	probnonDLT<-ifelse(combine01==TRUE, phat0or1+phat2, phat0+phat1+phat2)
}

########################################
#Plotting New Model#####################
########################################
if (plotit==TRUE)	{
if (discrete==FALSE)	{
	if (combine01==FALSE)	{
		graphics::plot(c(0,(max(dall)*1.05)),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey85",colors)
			colors<-ifelse(collectedtox==2,"grey20",colors)
			colors<-ifelse(collectedtox==3,"grey60",colors)
			colors<-ifelse(collectedtox==4,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			type<-ifelse(collectedtox==4,19,type)
			area<-ifelse(collectedtox==3 | collectedtox==4, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp0<- reg$coef['Intercept']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))
		tmp1<- reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))
		tmp2<- reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))
		tmp3<- reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))
		graphics::lines(seq(0,(max(dall)*1.05),1),(1/(1+exp(tmp0))),lty=1, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),(1/(1+exp(tmp3)))*(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		graphics::points(dose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x=(dose*1.05), y=(targetDLT+.03), labels=dose)
		graphics::legend(0,.95,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Next Dose"),
			lty=c(1,2,3,4,-1),pch=c(-1,-1,-1,-1,15),lwd=c(2,2,2,2,1),
			col=c(1,1,1,1,"red"),bg="white")
		graphics::legend("bottom",c("Grade 0","Grade 1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey82","grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Continuation Ratio Model CRM")
		}
	if (combine01==TRUE)	{
		graphics::plot(c(0,(max(dall)*1.05)),c(-0.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey60",colors)
			colors<-ifelse(collectedtox==2,"grey85",colors)
			colors<-ifelse(collectedtox==3,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			area<-ifelse(collectedtox==2 | collectedtox==3, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,(max(dall)*1.05),1))))))
		graphics::lines(seq(0,(max(dall)*1.05),1),1-crequal01,lty=2, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),(1-crequal01)*(1-crequal2),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,(max(dall)*1.05),1),(1-crequal01)*(1-crequal2)*(1-crequal3),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		graphics::points(dose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x=(dose*1.05), y=(targetDLT+.03), labels=dose)
		graphics::legend(0,.95,c("p(y>=2)","p(y>=3)","p(y>=4)","Next Dose"),
			lty=c(2,3,4,-1),pch=c(-1,-1,-1,15),lwd=c(2,2,2,1),
			col=c(1,1,1,"red"),bg="white")
		graphics::legend("bottom",c("Grade 0/1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Continuation Ratio Model CRM")
		}
	}
if (discrete==TRUE)	{
	if (combine01==FALSE)	{
		lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
		graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0)))),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey20",colors)
			colors<-ifelse(collectedtox==2,"grey60",colors)
			colors<-ifelse(collectedtox==3,"grey85",colors)
			colors<-ifelse(collectedtox==4,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			type<-ifelse(collectedtox==4,19,type)
			area<-ifelse(collectedtox==3 | collectedtox==4, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		tmp0<- reg$coef['Intercept']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp1<- reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp2<- reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		tmp3<- reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1/(1+exp(tmp0))),lty=1, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=2, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1/(1+exp(tmp3)))*(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		for (i in 1:length(discretedoses))	{
		graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
		graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
		graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
		graphics::legend(0,.95,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Next Discrete Dose","Estimated Toxicity"),
			lty=c(1,2,3,4,-1,-1),pch=c(-1,-1,-1,-1,15,15),lwd=c(2,2,2,2,1,1),
			col=c(1,1,1,1,"red","gray90"),bg="white")
		graphics::legend("bottom",c("Grade 0","Grade 1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey82","grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Continuation Ratio Model CRM")
		}
	if (combine01==TRUE)	{
		lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
		graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0)))),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
			colors<-ifelse(collectedtox==0,1,NA)
			colors<-ifelse(collectedtox==1,"grey60",colors)
			colors<-ifelse(collectedtox==2,"grey85",colors)
			colors<-ifelse(collectedtox==3,1,colors)
			type<-ifelse(collectedtox==0,1,NA)
			type<-ifelse(collectedtox==1,19,type)
			type<-ifelse(collectedtox==2,19,type)
			type<-ifelse(collectedtox==3,19,type)
			area<-ifelse(collectedtox==2 | collectedtox==3, 1, 0)
			graphics::points(jitter(collecteddose),area,pch=type,col=colors)
		crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))))))
		crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))))))
		crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))))))
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1-crequal01,lty=2, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1-crequal01)*(1-crequal2),lty=3, col=1, lwd=2)
		graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),(1-crequal01)*(1-crequal2)*(1-crequal3),lty=4, col=1, lwd=2)
		graphics::abline(h=targetDLT,lty=3)
		graphics::abline(v=dosenoconstraint,lty=3)
		graphics::abline(h=0,lty=1)
		for (i in 1:length(discretedoses))	{
			graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
			graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
		graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
		graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
		graphics::legend(0,.95,c("p(y>=2)","p(y>=3)","p(y>=4)","Next Discrete Dose","Estimated Toxicity"),
			lty=c(2,3,4,-1,-1),pch=c(-1,-1,-1,15,15),lwd=c(2,2,2,1,1),
			col=c(1,1,1,"red","gray90"),bg="white")
		graphics::legend("bottom",c("Grade 0/1","Grade 2","Grade 3","Grade 4"),
			pch=c(1,19,19,19,19),col=c(1,"grey20","grey60",1),bg="white", horiz=TRUE)
		graphics::title("Updated Continuation Ratio Model CRM")
		}
	}
}
}

###################################################
###Binary CRM Design ##############################
###################################################
if (design=='CRM')	{

lastdose<-collecteddose[length(collecteddose)]
yall<-append(pseudotox, collectedtox, after=length(pseudotox))
dall<-append(pseudodose, collecteddose, after=length(pseudodose))
weight<-c(rep(pseudoweights/length(pseudotox), length(pseudotox)),
	rep(1, (length(collectedtox))))
weightpseudo<-round(((sum(weight[1:length(pseudotox)]))/(length(collectedtox)+pseudoweights))*100,2)
weightdata<-round((sum(weight[(length(pseudotox)+1):length(weight)])/(length(collectedtox)+pseudoweights))*100,2)
reg<-suppressWarnings(rms::lrm(yall~dall, weights=weight))

	dose<-round((log(targetDLT/(1-targetDLT))-reg$coef['Intercept'])/reg$coef['dall'])
	dosenoconstraint<-dose
	dltlastcohort<-collectedtox[((length(collectedtox)-cohortsize)+1):length(collectedtox)]
	actualdltcohort<-sum(dltlastcohort)
	a<-constraints(dose, lastdose, actualdltcohort, numberdltrule,
		lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum)
	dose<-a$dosehat.new
	constraintused<-a$constraintused
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
if (discrete==FALSE)	{
	esttoxdiscdose<-NA
	discdose<-NA
}

	prob1<- ifelse(discrete==FALSE, reg$coef['Intercept']+ (reg$coef['dall']*dose),
		reg$coef['Intercept']+ (reg$coef['dall']*discdose))
	prob1<-1/(1+exp(-prob1))
	prob0<-1-prob1
	probDLT<-round(prob1*100,2)
	probnonDLT<-round(prob0*100,2)
	phat0or1<-phat0<-phat1<-phat2<-phat3<-phat4<-NA

if (plotit==TRUE)	{
if (discrete==FALSE)	{
	graphics::plot(c(0,(max(dall)*1.05)),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
	colors<-ifelse(collectedtox==0,1,1)
	type<-ifelse(collectedtox==0, 1, 19)
	graphics::points(jitter(collecteddose),collectedtox,pch=type,col=colors)
	crm<- reg$coef['Intercept']+(reg$coef['dall']*seq(0,(max(dall)*1.05),1))
	graphics::lines(seq(0,(max(dall)*1.05),1),1/(1+exp(-crm)),lty=3, col=1, lwd=2)
	graphics::abline(h=targetDLT,lty=3)
	graphics::abline(v=dosenoconstraint,lty=3)
	graphics::abline(h=0,lty=1)
	graphics::points(dose,targetDLT,pch=15,cex=1.3,col="red",bg="red")
	graphics::text(x=(dose*1.05), y=(targetDLT+.03), labels=dose)
	graphics::legend(0,.95,c("Probability of a DLT","Next Dose","Patient DLTs", "Patient Non-DLTs"),
		lty=c(3,-1,-1,-1),pch=c(-1,15,19,1),lwd=c(2,1,1,1), col=c(1,"red",1,1),bg="white")
	graphics::title("Updated Model for Binary CRM")
	}
if (discrete==TRUE)	{
	lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
	graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0)))),c(-.1,1), type="n",ylab="Probability",xlab="Dose")
	colors<-ifelse(collectedtox==0,1,1)
	type<-ifelse(collectedtox==0, 1, 19)
	graphics::points(jitter(collecteddose),collectedtox,pch=type,col=colors)
	crm<- reg$coef['Intercept']+(reg$coef['dall']*seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1))
	graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(max(collecteddose)*1.05,0))),1),1/(1+exp(-crm)),lty=3, col=1, lwd=2)
	graphics::abline(h=targetDLT,lty=3)
	graphics::abline(v=dosenoconstraint,lty=3)
	graphics::abline(h=0,lty=1)
	for (i in 1:length(discretedoses))	{
		graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
		graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
	graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
	graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
	graphics::legend(0,.95,c("Probability of a DLT","Next Dose","Patient DLTs", "Patient Non-DLTs","Estimated Toxicity"),
		lty=c(3,-1,-1,-1,-1),pch=c(-1,15,19,1,15),lwd=c(2,1,1,1,1), col=c(1,"red",1,1,"gray90"),bg="white")
	graphics::title("Updated Model for Binary CRM")
	}
}
}
if (discrete==TRUE)	{
esttoxdiscdose2<-c()
	for (k in 1:length(discretedoses))	{
	esttoxdiscdose2[k]<-round(esttoxdiscdose[k]*100,2)}
	}
if (discrete==FALSE)	{
	esttoxdiscdose2<-NA}

	return(list("Next Dose with No Used Constraints"=dosenoconstraint,
	"Next Dose with Constraints if Applicable"=dose,
	"Next Dose if Using Discrete Dose Levels"=discdose,
	"Estimated Probability of a non-DLT"=probnonDLT,
	"Estimated Probability of a DLT"=probDLT,
	"Estimated Probability of a Grade 0 or 1 Toxicity"=phat0or1,
	"Estimated Probability of a Grade 0 Toxicity"=phat0,
	"Estimated Probability of a Grade 1 Toxicity"=phat1,
	"Estimated Probability of a Grade 2 Toxicity"=phat2,
	"Estimated Probability of a Grade 3 Toxicity"=phat3,
	"Estimated Probability of a Grade 4 Toxicity"=phat4,
	"Estimated Regression Model"=reg,
	"Estimated Toxicity at Discrete Doses"=esttoxdiscdose2,
	"Discrete Doses"=discretedoses,
	"Influence from Pseudodata"=weightpseudo,
	"Influence from Collected Data"=weightdata,
	"Constraints Used"=constraintused))
}

