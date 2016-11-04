#' Simulations to Assess operating characteristics of proportional odds model (POM), continuation ratio (CR) Model, and binary CRM Designs
#'
#' This function performs simulations to assess design performance characteristics of the proportional odds model (POM),
#' continuation ratio (CR) Model, and binary CRM likelihood-based dose finding designs.
#' There are many different options in this function
#' to vary sample size, cohort sizes, target dose-limiting toxicity (DLT) rates, true underlying dose-toxicity models,
#' discrete or continuous dose levels, combining ordianl grades 0 and 1 into one category, and the design you wish to assess.
#' Returns results collected over all simulations run for these specified safety and dose estimation performance criteria.
#' @param startdose Starting dose for a trial for a specified target DLT rate.
#' @param numbersims Number of simulations you wish to run (We used 2000 simulations in our studies).
#' @param cohortsize Number of patients treated in each cohort. Defaults to 3.
#' @param samplesize Total sample size of each trial. Defaults to 30.
#' @param pseudoweights Pseudoweights determine the amount of influence the pseudodata has on estimating the updated model.
#' Please refer to Van Meter et al (2010) for more information. We suggest setting the pseudoweights equal to the cohortsize.
#' When pseudoweights = X, the total pseudodata represent X individuals in the updated model. If not specified,
#' pseudoweights defaults to cohortsize.
#' @param stopearly True/False. Allows early stopping if a certain number of patients are treated at the same dose level and
#' the model estimates the same dose to treat the next cohort.
#' Must specify stopearlynumber when stopearly = TRUE. Defaults to FALSE.
#' @param pseudotox A vector of toxicity grades specified by the pseudodata.
#' @param pseudodose A vector of dose levels specified by the pseudodata.
#' @param dosetox Specifies the true underlying dose-toxicity relationship model. Please refer to
#' Van Meter et al (2010) for specification of the POM and CR model.
#' Assumes 5 toxicity categories as specified by CTCAEv4.0, i.e. grades 0, 1, 2,
#' 3, and 4. Grade 5 toxicities (death directly related to an adverse event are not considered when running simulations).\cr
#' If using a POM, it must be in the form:
#'  c(\eqn{\alpha1, \beta1, \alpha2, \beta2, \alpha3, \beta3, \alpha4, \beta4}).\cr
#' If using a CR model it must be in the form: c(\eqn{\alpha, \gamma0, \theta1,
#' \gamma1, \theta2, \gamma2, \theta3, \gamma3}).
#' \strong{IMPORTANT!} If combine01 is true, you still must specify the dosetox model in
#' terms of grades 0, 1, 2, 3, and 4 for the simulations to run correctly!!!
#' @param truedosetoxmodeltype Must specify what type of true underlying dose-toxicity model you are using.
#' Choices are: POM, CR, or CRM.
#' @param design Specifies which dose finding design you are running simulations on. Choices are: POM, CR, or CRM.
#' @param targetDLT Target dose limiting toxicity (DLT) rate pre-specified by clinical investigators prior to the start of the trial.
#' Must be specified between 0 and 1. Defaults to 0.30.
#' @param discrete True/False. If discrete = TRUE this allows for discrete dose levels
#' to be specified prior to the start of the trial. Defaults to FALSE
#' @param discretedoses Specified discrete dose levels if desired. They must be specified if discrete is equal to TRUE.
#' It is written for j dose levels as: c(d1, d2,..., dj).
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
#' @param minimum Minimum dose that will be considered to test patients. Must be specified if lowerlimitrule<1.
#' @param maximum Maximum dose that will be considered to test patients. Must be specified if upperlimitrule<1.
#' @param combine01 True/False. If combine01 = TRUE, toxicity grades 0 and 1 are combined into 1 category.
#' Therefore all toxicities must be coded: 0 (grades 0 and 1), 1 (grade2), 2 (grade 3),
#' and 3 (grade 4) according to CTCAEv4.0. Defaults to FALSE.
#' @param stopearlynumber Total number of patients treated at the same dose level before stopping.
#' Must set stopearly = TRUE when a number is specified.
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
#' @param rounddown True/False. Only applicable when using discrete dose levels.\cr
#' If rounddown = TRUE, the estimate dose from specified model will round down to the more conservative discrete dose level.\cr
#' If rounddown = FALSE, it will select the discrete dose closest to the estimated model selection. Defaults to FALSE.
#'
#' @details
#' If using a POM CRM, this function assumes a proportional odds model as described by Harrell (Harrell, 2001).
#' For combine01=FALSE, y equal to toxicity grade outcomes j in c(0, 1, 2, 3, 4) as specified by CTCAEv4.0, and x equal to the dose,
#' this is written in the form:\cr
#' \deqn{P(y>=j|x)=1/(1+exp(-(\alpha_j + \beta*x))), j = 1, 2, 3, 4}\cr
#' For combine01=TRUE, toxicity grades are now 0/1, 2, 3, and 4, and y is recoded as 0 = grades 0 and 1, 1 = grade 2,
#' 2 = grade 3, and 3 = grade 4.\cr
#'
#' If using a CR model design, this function assumes a continuation ratio model as described by Harrell (Harrell, 2001).
#' For combine01=FALSE, y equal to toxicity grade outcomes j in c(0, 1, 2, 3, 4) as specified by CTCAEv4.0, and x equal to the dose,
#' this is written in the form:\cr
#' \deqn{P(y=j|y>=j,x)=1/(1+exp(-(\alpha + \theta_j + \gamma*x))), j = 0, 1, 2, 3}\cr
#' For combine01=TRUE, toxicity grades are now 0/1, 2, 3, and 4, and y is recoded as 0 = grades 0 and 1, 1 = grade 2,
#' 2 = grade 3, and 3 = grade 4.\cr
#'
#' If using a binary CRM, this assumes a standard 2-parameter logistic model.
#' For y equal to toxicity outcome (1 is a DLT, 0 is not a DLT) and x equal to the dose, it is written in the form:\cr
#' \deqn{P(y=1|x)=1/(1+exp(-(\alpha + \beta*x)))}\cr
#'
#' Estimated next dose is calculated by Pr(3 or 4 toxicity) <= targetDLT for the model utilized.
#'
#' @return Returns various safety and dose estimation performance statistics collected over all simulations.
#' Please see the list below for details.
#' @return \item{`Median Total Sample Size of Trials Not Stopped Early Due to Safety Concerns`}{Median
#' total sample size of trial for all simulations not stopped early due to ANY safety constraints.}
#' @return \item{`Patients per Cohort`}{Pre-defined cohortsize}
#' @return \item{`Proportion of Trials Stopped Early Due to Safety Concerns`}{Proportion of trials
#' stopped early due to ANY safety constraints among all simulations.}
#' @return \item{`Proportion of Trials Stopped Early Used Constraint`}{Proportion of trials stopped
#' early due to any constraints while estimating the final dose for all simulations.}
#' @return \item{`Proportion of Trials Stopped Early Used the DLT Constraint`}{Proportion of trials
#' stopped early only due to the DLT constraint while estimating the final dose for all simulations.}
#' @return \item{`Proportion of Trials Stopped Early Used the Increase Constraint`}{Proportion of trials
#' stopped early only due to the increase constraint while estimating the final dose for all simulations.}
#' @return \item{`25\% Quantile Sample Size`}{25\% quantile total Sample Size of trial for all simulations.}
#' @return \item{`75\% Quantile Sample Size`}{75\% quantile total Sample Size of trial for all simulations.}
#' @return \item{`Median Dose`}{Median dose estimated from all simulations.}
#' @return \item{`25\% Quantile Dose`}{25\% quantile dose estimated from all simulations.}
#' @return \item{`75\% Quantile Dose`}{75\% quantile dose estimated from all simulations.}
#' @return \item{`Percentage of Times Each Discrete Dose Level Was Selected`}{This value is only returned when "discrete" is TRUE and
#' "discretedoses" are pre-specified. Each discrete dose level has a percentage indicating the chance that the dose level
#' been selected as the final estimated dose for all simulations.}
#' @return \item{`Median Percent Difference Between Estimated Dose and MTD`}{Median difference in percentage between the final
#' estimated dose and the true maximum tolerated dose (MTD) for all simulations.}
#' @return \item{`Median Expected DLT for the Final Estimated Dose`}{Median expected DLT rate
#' for the final estimated dose for all simulations.}
#' @return \item{`Median Sample Size For All Trials`}{Median total sample size of trial for all simulations.}
#' @return \item{`Median \% of patients in all trials treated at doses with >40\% DLT rate`}{Median percentage of patients
#' in the trial treated at doses with a DLT rate larger than 40\% for all simulations.}
#' @return \item{`Median \% of patients in all trials treated at doses with >50\% DLT rate`}{Median percentage of patients
#' in the trial treated at doses with a DLT rate larger than 50\% for all simulations.}
#' @return \item{`Median \% of patients in all trials treated at doses with <20\% DLT rate`}{Median percentage of patients
#' in the trial treated at doses with a DLT rate less than 20\% for all simulations.}
#' @return \item{`Median \% of patients in all trials treated at doses with <10\% DLT rate`}{Median percentage of patients
#' in the trial treated at doses with a DLT rate less than 10\% for all simulations.}
#' @return \item{`Median Number of patients in all trials with a DLT (grade 3 or 4)`}{Median number of patients
#' in trial with a grade 3 or grade 4 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with a DLT (grade 3 or 4)`}{Median percentage of patients
#' in trial with a grade 3 or grade 4 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with grade 0 and 1`}{Median percentage of patients
#' in trial with a grade 1 or grade 2 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with a non-DLT (grade 1 or 2)`}{Median percentage of patients
#' in trial with a grade 1 or grade 2 non-DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with no toxicity`}{Median percentage of patients
#' in trial without DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with grade 1`}{Median percentage of patients
#' in trial with grade 1 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with grade 2`}{Median percentage of patients
#' in trial with grade 2 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with grade 3`}{Median percentage of patients
#' in trial with grade 3 DLT for all simulations.}
#' @return \item{`Median \% of patients in all trials with grade 4`}{Median percentage of patients
#' in trial with grade 4 DLT for all simulations.}
#' @return \item{`Acutal MTD`}{True MTD specified from dosetox.}
#' @return \item{`True MTD Discrete`}{This value is only returned when "discrete" is TRUE and
#' "discretedoses" are pre-specified. Rounds the "Acutal MTD" to the nearest pre-specified discrete dose.}
#' @return \item{Dose}{A vector of final estimated dose from each simulation.}
#' @return \item{Listit}{A matrix of estimated doses after each cohort for each simulation.
#' For example, if "numbersims" is 5 and there are 10 cohorts in each trial, this will result in a 5*10 matrix.}

#' @author Emily V. Dressler, PhD \cr
#' Markey Cancer Center \cr
#' Division of Cancer Biostatistics \cr
#' University of Kentucky \cr
#' \email{EmilyVDressler@@gmail.com}
#' @references 1. Van Meter EM, Garrett-Mayer E, Bandyopadhyay D. Proportional odds model for dose finding clinical trial designs with ordinal toxicity grading. Statistics in Medicine 2011; 30: 2070-2080. \cr
#'2. Van Meter EM, Garrett-Mayer E, Bandyopadhyay. Dose finding clinical trial design for ordinal toxicity grades using the continuation ratio model: an extension of the continual reassessment method. Clinical Trials 2012; 9(3): 303-313. \cr
#'3. Garrett-Mayer E. The continual reassessment method for dose-finding studies: a tutorial. Clinical Trials 2006; 3: 57-71. \cr
#'4. Piantadosi S, Fisher JD, Grossman S. Practical implementation of a modified continual reassessment method for dose-finding trials. Cancer Chemother Pharmacol 1998; 41: 429-436. \cr
#'5. Harrell FE, Jr. Regression Modeling Strategies with Application to Linear Models, Logistic Regression, and Survival Analysis. Springer: New York, NY, 2001. \cr
#'6. McCullagh P. Regression Models for Ordinal Data. Journal of the Royal Statistical Society. Series B (Methodological). 1980; 42: 109-142. \cr
#'7. CTCAE. Cancer Therapy Evaluation Program, Common Terminology Criteria for Adverse Events, Version 4.0, DCTD, NCI, NIH, DHHS (http://ctep/cancer.gov). In Cancer Therapy Evaluation Program, 2010.
#' @export
#' @examples
#' #Example 1
#' #Underlying True POM dosetox
#' truePO1<-c(-0.4, 0.0011, -1.3, 0.0011, -2.8, 0.0011, -3.9, 0.0011)
#' #POM Pseudodata
#' initial.pom.y1<-c(rep(0,45),rep(1,36),rep(2,9),rep(3,8),rep(4,2),rep(0,18),
#' rep(1,35),rep(2,17),rep(3,24),rep(4,6),rep(0,8),rep(1,24),rep(2,18),rep(3,38),
#' rep(4,12),rep(0,1),rep(1,4),rep(2,5),rep(3,35),rep(4,55))
#' initial.pom.d1<-c(rep(200,100),rep(1067,100),rep(1613,100),rep(3000,100))
#' #POM CRM Design Simulation run 10 times, with a true dosetox model equal to a POM
#' crmsimulations(startdose = 1060, numbersims = 10, cohortsize = 3, samplesize = 30,
#' pseudoweights = 3,pseudotox = initial.pom.y1, pseudodose = initial.pom.d1, dosetox = truePO1,
#' truedosetoxmodeltype = 'POM', design = 'POM', targetDLT = 0.30, discrete = FALSE,
#' discretedoses = NA, numberdltrule = NA, lowerlimitrule = NA, upperlimitrule = NA,
#' dltrule = NA, increaserule= NA, minimum = NA, maximum = NA, combine01 = FALSE)
#'
#' #Example 2
#'#Underlying true POM dosetox
#'truePO1<-c(-0.4, 0.0011, -1.3, 0.0011, -2.8, 0.0011, -3.9, 0.0011) #POM 1#
#'#CR Model Pseudodata with Grades 0 and 1 combined into one category
#'combine.cr.y1<-c(rep(0,81),rep(1,9),rep(2,8),rep(3,2),
#'rep(0,55),rep(1,15),rep(2,26),rep(3,4),
#'rep(0,37),rep(1,13),rep(2,40),rep(3,10),
#'rep(0,5),rep(1,5),rep(2,35),rep(3,55))
#'initial.cr.d1<-c(rep(200,100),rep(934,100),rep(1467,100),rep(3000,100))
#'#CR Model CRM Design Simulation run 10 times, with a true dosetox model equal to a POM
#'crmsimulations(startdose = 1060, numbersims = 10, cohortsize = 3,
#'samplesize = 30, pseudoweights = 3,pseudotox = combine.cr.y1,
#'pseudodose = initial.cr.d1, dosetox = truePO1, truedosetoxmodeltype = 'POM',
#'design = 'CR', targetDLT = 0.30, discrete = FALSE, discretedoses = NA,
#'numberdltrule = NA, lowerlimitrule = NA, upperlimitrule = NA, dltrule = NA,
#'increaserule = NA, minimum = NA, maximum = NA, combine01 = TRUE)
#'
#'#Example 3
#'#Underlying True CR dosetox
#'crmodel1<-c(1.58,-0.0018,-0.45,-0.0018,1.20,-0.0018,3.14,-0.0018) #CR Model 1#
#'#CRM Pseudodata
#'initial.crm.y1<-c(rep(0,90),rep(1,10),rep(0,70),rep(1,30),
#'rep(0,50),rep(1,50),rep(0,10),rep(1,90))
#'initial.crm.d1<-c(rep(200,100),rep(1060,100),rep(1600,100),rep(3000,100))
#'#Binary CRM Design Simulation run 10 times, with a true dosetox model equal to a CR
#Sample size 20, cohortsize 2, and discrete doses
#'crmsimulations(startdose = 1060, numbersims = 10, cohortsize = 2,
#'samplesize = 20, pseudoweights = 3,pseudotox = initial.crm.y1,
#'pseudodose = initial.crm.d1, dosetox = crmodel1, truedosetoxmodeltype = 'CR',
#'design = 'CRM', targetDLT = 0.30, discrete=TRUE,
#'discretedoses=c(200,400,800,1200,1800,2500),
#'numberdltrule = NA, lowerlimitrule = NA, upperlimitrule = NA, dltrule = NA,
#'increaserule = NA, minimum = NA, maximum = NA, combine01 = FALSE)
#'


crmsimulations <-
function(startdose, numbersims, cohortsize=3, samplesize=30, pseudoweights=NA, stopearly=FALSE, stopearlynumber=NA,
	pseudotox, pseudodose, dosetox, truedosetoxmodeltype, design, targetDLT=0.30,
	discrete=FALSE, discretedoses=NA, numberdltrule=NA, lowerlimitrule=NA, upperlimitrule=NA,
	dltrule=NA, increaserule=NA, minimum=NA, maximum=NA, combine01=FALSE, rounddown=FALSE)
	{

if (is.na(pseudoweights))	{
	pseudoweights<-cohortsize
}

#WARNINGS#
if ((discrete==FALSE) & stopearly==TRUE)	{
	stop('Stopping early currently only applied to discrete dose levels')
}
if (is.na(stopearlynumber) & stopearly==TRUE)	{
	stop('If stopping early, you must specify how many patients treated at a dose level before stopping')
}
if ((is.na(stopearlynumber)==FALSE) & stopearly==FALSE)	{
	stop('If specifying a number treated at a dose level to stop early, you must set stopearly=TRUE')
}
if ((discrete==FALSE) & rounddown==TRUE)	{
	stop('Rounding down only applied to discrete dose levels')
}
if ((is.na(lowerlimitrule)==FALSE & is.na(upperlimitrule)==TRUE) | (is.na(upperlimitrule)==FALSE & is.na(lowerlimitrule)==TRUE))	{
	stop('For these simulations, you must specify both a lower and upper bound rule')}
if ((is.na(lowerlimitrule)==FALSE & is.na(upperlimitrule)==FALSE) & (((lowerlimitrule<1 & lowerlimitrule>0)
	& (upperlimitrule>=1 | upperlimitrule==0)) |
	((lowerlimitrule==0 | lowerlimitrule>=1) & (upperlimitrule>0 & upperlimitrule<1))))	{
	stop("For these simulations, lowerlimitrules and upperlimitrules must either be both in terms
	of percentages or a specifed dose level")
}
if (design=='CRM' & combine01==TRUE)	{
	stop("If running simulations on the binary CRM, combine01 must be FALSE")
}
if (discrete==TRUE)	{
numalreadytested2<-matrix(nrow=numbersims, ncol=length(discretedoses))}
if (discrete==FALSE)	{
numalreadytested2<-NA}
listit<-matrix(nrow=numbersims,ncol=samplesize/cohortsize)
sim<-under20.n<-under20.pct<-over40.n<-
over40.pct<-over50.n<-over50.pct<-dlts.n<-dlts.pct<-under10.pct<-under10.n<-ts1.n<-
ts1.pct<-ts2.n<-ts2.pct<-ts3.n<-ts3.pct<-ts4.n<-ts4.pct<-
allunder20.n<-allunder20.pct<-allover40.n<-allover40.pct<-allover50.n<-allover50.pct<-
alldlts.n<-alldlts.pct<-allunder10.pct<-allunder10.n<-allts1.n<-allts1.pct<-
allts2.n<-allts2.pct<-allts3.n<-allts3.pct<-allts4.n<-allts4.pct<-N<-Nna<-
if10percent<-if20percent<-if40<-if50<-if10<-
if20<-constraintused<-dltruleused<-incrruleused<-lowruleused<-
upruleused<-numberdlty<-expectedDLT<-
allts01.n<-allts01.pct<-ts01.n<-ts01.pct<-rep(0,numbersims)

DLT<-targetDLT
dosehat<-startdose
initial.y<-pseudotox
initial.d<-pseudodose

#Identifying true doses for different DLT rates#
if (truedosetoxmodeltype=='CRM')	{
truedose<- round((log(targetDLT/(1-targetDLT))-dosetox[1])/dosetox[2],0)
truedlt40<-round((log(.4/(1-.4))-dosetox[1])/dosetox[2],0)
truedlt50<-round((log(.5/(1-.5))-dosetox[1])/dosetox[2],0)
truedlt10<-round((log(.1/(1-.1))-dosetox[1])/dosetox[2],0)
truedlt20<-round((log(.2/(1-.2))-dosetox[1])/dosetox[2],0)
}

if (truedosetoxmodeltype=='POM')	{
truedose<- round((log(targetDLT/(1-targetDLT))-dosetox[5])/dosetox[6],0)
truedlt40<-round((log(.4/(1-.4))-dosetox[5])/dosetox[6],0)
truedlt50<-round((log(.5/(1-.5))-dosetox[5])/dosetox[6],0)
truedlt10<-round((log(.1/(1-.1))-dosetox[5])/dosetox[6],0)
truedlt20<-round((log(.2/(1-.2))-dosetox[5])/dosetox[6],0)
}

if (truedosetoxmodeltype=='CR')	{
		truedltlevel<-c(targetDLT,.4,.5,.1,.2)
		d<-(-((1-truedltlevel)/(truedltlevel)))
		c<-(exp(dosetox[1])+exp(dosetox[1]+dosetox[3])+
			exp(dosetox[1]+dosetox[5]))
		b<-(exp((2*dosetox[1])+dosetox[3])+
			exp((2*dosetox[1])+dosetox[5])+
			exp((2*dosetox[1])+dosetox[3]+dosetox[5]))
		a<-(exp((3*dosetox[1])+dosetox[3]+dosetox[5]))
		m<-((2*(b^3))-(9*a*b*c)+(27*(a^2)*d))
		n2<-(m^2)-(4*(((b^2)-(3*a*c))^3))
		n<-sqrt(abs(n2))
		n<-ifelse(n2>=0, n, as.complex(n*(1i)))
		xx<-ifelse(n2>=0, ifelse((.5*(m+n))>0,(.5*(m+n))^(1/3),-(-(.5*(m+n)))^(1/3)),
			as.complex((.5*(m+n))^(1/3)))
		yy<-ifelse(n2>=0, ifelse((.5*(m-n))>0,(.5*(m-n))^(1/3),-(-(.5*(m-n)))^(1/3)),
			as.complex((.5*(m-n))^(1/3)))
		firstterm<-(-1*b)/(3*a)
		coef1<-as.complex((1+(1i*sqrt(3)))/(6*a))
		coef2<-as.complex((1-(1i*sqrt(3)))/(6*a))
		xi1<-as.complex(firstterm-(xx/(3*a))-(yy/(3*a)))
		xi2<-as.complex(firstterm+(coef1*xx)+(coef2*yy))
		xi3<-as.complex(firstterm+(coef2*xx)+(coef1*yy))

	xi<-blah1<-blah<-dcheck<-possibledoses<-d.new<-truedose1<-vector("list",5)
	for (i in 1:5)	{
		xi[[i]]<-c(xi1[i], xi2[i], xi3[i])
		blah1[[i]]<-(log(xi[[i]])/dosetox[4])
		blah[[i]]<-suppressWarnings(as.double(blah1[[i]]))
		dcheck[[i]]<-1/((1+exp(dosetox[1]+(dosetox[2]*blah[[i]])))*
			(1+exp(dosetox[1]+dosetox[3]+(dosetox[4]*blah[[i]])))*
			(1+exp(dosetox[1]+dosetox[5]+(dosetox[6]*blah[[i]]))))
		possibledoses[[i]]<-ifelse((dcheck[[i]]>=truedltlevel[i]-0.01 & dcheck[[i]]<=truedltlevel[i]+0.01),blah[[i]],NA)
		d.new[[i]]<-blah[[i]][which(is.na(possibledoses[[i]])==FALSE)]
		truedose1[[i]]<-ifelse(length(d.new[[i]])==1, round(d.new[[i]],0), NA)
	}
	truedose<- truedose1[[1]]
	truedlt40<-truedose1[[2]]
	truedlt50<-truedose1[[3]]
	truedlt10<-truedose1[[4]]
	truedlt20<-truedose1[[5]]
}

numberiterations<-samplesize/cohortsize

#####################
#RUNNING SIMULATIONS#
#####################
i<-0
while((i<-i+1)<=numbersims) {

if (stopearly==FALSE)	{
###Upper and Lower Rules are in terms of Percentages###
if (is.na(lowerlimitrule)==FALSE)	{
if ((lowerlimitrule>0 & lowerlimitrule<1) | (upperlimitrule>0 & upperlimitrule<1)) 	{
stop<-0
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

for (j in 2:numberiterations)	{
if (stop==0)	{
	it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=j, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	if ((it$dose==(minimum+((maximum-minimum)*lowerlimitrule)) | it$dose==(maximum-((maximum-minimum)*upperlimitrule)))
		& it$dose==it$d[length(it$d)])	{
				it$dose<-it$bse<-it$var.dosehat<-NA
				stop<-1
	}
	listit[i,j]<-it$dose
	if (j==numberiterations)	{
		stop<-1}
	}
}
}

###Upper and Lower Limit Rules are numbers###
if ((lowerlimitrule>=1 | lowerlimitrule==0) | (upperlimitrule==0 | upperlimitrule>=1))	{
	stop<-0
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

	for (j in 2:numberiterations)	{
	if (stop==0)	{
		it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
			dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
			numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum, combine01, iteration=j, rounddown)
		it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	if ((it$dose==lowerlimitrule | it$dose==upperlimitrule) & it$dose==it$d[length(it$d)])	{
				it$dose<-it$bse<-it$var.dosehat<-NA
				stop<-1
	}
	listit[i,j]<-it$dose
	if (j==numberiterations)	{
		stop<-1}
	}
}
}
}

###No upper or lower limits###
if (is.na(lowerlimitrule)==TRUE & is.na(upperlimitrule)==TRUE)	{
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

	for (j in 2:numberiterations)	{
		it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
			dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
			numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum, combine01, iteration=j, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,j]<-it$dose
	}
}
}

if (stopearly==TRUE)	{
###Upper and Lower Rules are in terms of Percentages###
if (is.na(lowerlimitrule)==FALSE)	{
if ((lowerlimitrule>0 & lowerlimitrule<1) | (upperlimitrule>0 & upperlimitrule<1)) 	{
stop<-0
stopearly1<-0
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

for (j in 2:numberiterations)	{
if (stop==0 & stopearly1==0)	{
	it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=j, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	if ((it$dose==(minimum+((maximum-minimum)*lowerlimitrule)) | it$dose==(maximum-((maximum-minimum)*upperlimitrule)))
		& it$dose==it$d[length(it$d)])	{
				it$dose<-it$bse<-it$var.dosehat<-NA
				stop<-1
	}
	if  (it$numalreadytested[which(it$discretedoses==it$discdose)]>=stopearlynumber)	{
	stopearly1<-1}
	listit[i,j]<-it$dose
	if (j==numberiterations)	{
		stop<-1}
	}
}
}

###Upper and Lower Limit Rules are numbers###
if ((lowerlimitrule>=1 | lowerlimitrule==0) | (upperlimitrule==0 | upperlimitrule>=1))	{
	stop<-0
	stopearly1<-0
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

	for (j in 2:numberiterations)	{
	if (stop==0 & stopearly1==0)	{
		it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
			dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
			numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum, combine01, iteration=j, rounddown)
		it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	if ((it$dose==lowerlimitrule | it$dose==upperlimitrule) & it$dose==it$d[length(it$d)])	{
				it$dose<-it$bse<-it$var.dosehat<-NA
				stop<-1
	}
	if  (it$numalreadytested[which(it$discretedoses==it$discdose)]>=stopearlynumber)	{
	stopearly1<-1}
	listit[i,j]<-it$dose
	if (j==numberiterations)	{
		stop<-1}
	}
}
}
}

###No upper or lower limits###
if (is.na(lowerlimitrule)==TRUE & is.na(upperlimitrule)==TRUE)	{
stopearly1<-0
	it<-fiterations(dosehat, cohortsize, pseudotox, pseudodose, pseudotox,
		dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
		numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
		increaserule, minimum, maximum, combine01, iteration=1, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	listit[i,1]<-it$dose

	for (j in 2:numberiterations)	{
	if (stopearly1==0)	{
		it<-fiterations(it$dose, cohortsize, pseudotox, it$d, it$y,
			dosetox, truedosetoxmodeltype, design, targetDLT, pseudoweights, discrete, discretedoses,
			numberdltrule, lowerlimitrule, upperlimitrule, dltrule,
			increaserule, minimum, maximum, combine01, iteration=j, rounddown)
	it$dose<-ifelse(discrete==TRUE,it$discdose,it$dose)
	if  (it$numalreadytested[which(it$discretedoses==it$discdose)]>=stopearlynumber)	{
	stopearly1<-1}
	listit[i,j]<-it$dose
	if (j==numberiterations)	{
		stop<-1}
	}
}
}}
#########
#RESULTS#
#########

#Gives y and d in terms of just the cohorts no pseudodata#
#Statistics with y and d consider all trials regardless of whether they stopped early or not#
	y<-yna<-it$y[-c(1:length(initial.y))]
	d<-dna<-it$d[-c(1:length(initial.y))]

#If the trial stopped early, assigns a NA to the yna and dna#
#So any statistics with a yna or dna only considers trials that reached total sample size#
#So y=yna when the trial did not stop early and same with d#
	if (is.na(it$dose))	{
	yna<-dna<-NA	}

	Nna[i]<-ifelse(is.na(it$dose),NA,length(y))
	N[i]<-length(y)
if (discrete==TRUE)	{
	numalreadytested2[i,]<- it$numalreadytested	}

#######################################################
#These results are collected inside of the simulations#

	###All statistics for all trials regardless of stopping early or not###
	allover40.n[i]<-sum(ifelse(d>truedlt40,1,0))
	allover40.pct[i]<-round(100*mean(ifelse(d>truedlt40,1,0)),digits=2)
	allover50.n[i]<-sum(ifelse(d>truedlt50,1,0))
	allover50.pct[i]<-round(100*mean(ifelse(d>truedlt50,1,0)),digits=2)
	allunder20.pct[i]<-round(100*mean(ifelse(d<truedlt20,1,0)),digits=2)
	allunder20.n[i]<-sum(ifelse(d<truedlt20,1,0))
	allunder10.pct[i]<-round(100*mean(ifelse(d<truedlt10,1,0)),digits=2)
	allunder10.n[i]<-sum(ifelse(d<truedlt10,1,0))

if (combine01==FALSE & design!='CRM')	{
	allts1.n[i]<-sum(ifelse(y==1,1,0))
	allts1.pct[i]<-round(100*mean(ifelse(y==1,1,0)),digits=2)
	allts2.n[i]<-sum(ifelse(y==2,1,0))
	allts2.pct[i]<-round(100*mean(ifelse(y==2,1,0)),digits=2)
	allts3.n[i]<-sum(ifelse(y==3,1,0))
	allts3.pct[i]<-round(100*mean(ifelse(y==3,1,0)),digits=2)
	allts4.n[i]<-sum(ifelse(y==4,1,0))
	allts4.pct[i]<-round(100*mean(ifelse(y==4,1,0)),digits=2)
	alldlts.n[i]<-sum(ifelse(y==3 | y==4, 1, 0))
	alldlts.pct[i]<-round(100*mean(ifelse(y==3 | y==4, 1, 0)),digits=2)
	allts01.n[i]<-NA
	allts01.pct[i]<-NA
	ts1.n[i]<-ifelse(is.na(yna[1]),NA, sum(ifelse(y==1,1,0)))
	ts1.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==1,1,0)),digits=2))
	ts2.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==2,1,0)))
	ts2.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==2,1,0)),digits=2))
	ts3.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==3,1,0)))
	ts3.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==3,1,0)),digits=2))
	ts4.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==4,1,0)))
	ts4.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==4,1,0)),digits=2))
	dlts.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==3 | y==4, 1, 0)))
	dlts.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==3 | y==4, 1, 0)),digits=2))
	ts01.n[i]<-NA
	ts01.pct[i]<-NA
}
if (combine01==TRUE & design!='CRM')	{
	allts1.n[i]<-NA
	allts1.pct[i]<-NA
	allts2.n[i]<-sum(ifelse(y==1,1,0))
	allts2.pct[i]<-round(100*mean(ifelse(y==1,1,0)),digits=2)
	allts3.n[i]<-sum(ifelse(y==2,1,0))
	allts3.pct[i]<-round(100*mean(ifelse(y==2,1,0)),digits=2)
	allts4.n[i]<-sum(ifelse(y==3,1,0))
	allts4.pct[i]<-round(100*mean(ifelse(y==2,1,0)),digits=2)
	alldlts.n[i]<-sum(ifelse(y==2 | y==3, 1, 0))
	alldlts.pct[i]<-round(100*mean(ifelse(y==2 | y==3, 1, 0)),digits=2)
	allts01.n[i]<-sum(ifelse(y==0,1,0))
	allts01.pct[i]<-round(100*mean(ifelse(y==0,1,0)),digits=2)
	ts1.n[i]<-NA
	ts1.pct[i]<-NA
	ts2.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==1,1,0)))
	ts2.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==1,1,0)),digits=2))
	ts3.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==2,1,0)))
	ts3.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==2,1,0)),digits=2))
	ts4.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==3,1,0)))
	ts4.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==3,1,0)),digits=2))
	dlts.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==2 | y==3, 1, 0)))
	dlts.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==2 | y==3, 1, 0)),digits=2))
	ts01.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==0,1,0)))
	ts01.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==0,1,0)),digits=2))
}
if (design=='CRM')	{
	allts1.n[i]<-NA
	allts1.pct[i]<-NA
	allts2.n[i]<-NA
	allts2.pct[i]<-NA
	allts3.n[i]<-NA
	allts3.pct[i]<-NA
	allts4.n[i]<-NA
	allts4.pct[i]<-NA
	alldlts.n[i]<-sum(ifelse(y==1, 1, 0))
	alldlts.pct[i]<-round(100*mean(ifelse(y==1, 1, 0)),digits=2)
	allts01.n[i]<-NA
	allts01.pct[i]<-NA
	ts1.n[i]<-NA
	ts1.pct[i]<-NA
	ts2.n[i]<-NA
	ts2.pct[i]<-NA
	ts3.n[i]<-NA
	ts3.pct[i]<-NA
	ts4.n[i]<-NA
	ts4.pct[i]<-NA
	dlts.n[i]<-ifelse(is.na(yna[1]),NA,sum(ifelse(y==1, 1, 0)))
	dlts.pct[i]<-ifelse(is.na(yna[1]),NA,round(100*mean(ifelse(y==1, 1, 0)),digits=2))
	ts01.n[i]<-NA
	ts01.pct[i]<-NA
}

	###Results only for trials that did not stop early###
	over40.n[i]<-ifelse(is.na(dna[1]),NA,sum(ifelse(d>truedlt40,1,0)))
	over40.pct[i]<-ifelse(is.na(dna[1]),NA,round(100*mean(ifelse(d>truedlt40,1,0)),digits=2))
	over50.n[i]<-ifelse(is.na(dna[1]),NA,sum(ifelse(d>truedlt50,1,0)))
	over50.pct[i]<-ifelse(is.na(dna[1]),NA,round(100*mean(ifelse(d>truedlt50,1,0)),digits=2))
	under20.pct[i]<-ifelse(is.na(dna[1]),NA,round(100*mean(ifelse(d<truedlt20,1,0)),digits=2))
	under20.n[i]<-ifelse(is.na(dna[1]),NA,sum(ifelse(d<truedlt20,1,0)))
	under10.pct[i]<-ifelse(is.na(dna[1]),NA,round(100*mean(ifelse(d<truedlt10,1,0)),digits=2))
	under10.n[i]<-ifelse(is.na(dna[1]),NA,sum(ifelse(d<truedlt10,1,0)))

	sim[i]<-ifelse(discrete==FALSE, it$dose, it$discdose)

if (truedosetoxmodeltype=='CRM')	{
		expectedDLT[i]<-100*(1/(1+exp(-(dosetox[1] + dosetox[2]*sim[i]))))
}
if (truedosetoxmodeltype=='POM')	{
		expectedDLT[i]<-100*(1/(1+exp(-(dosetox[5] + dosetox[6]*sim[i]))))
}
if (truedosetoxmodeltype=='CR')	{
		expectedDLT[i]<-100*(1/((1+exp(dosetox[1]+(dosetox[2]*sim[i])))*
        		(1+exp(dosetox[1]+dosetox[3]+(dosetox[4]*sim[i])))*
        		(1+exp(dosetox[1]+dosetox[5]+(dosetox[6]*sim[i])))))
}
	constraintused[i]<-ifelse(is.na(it$dose),NA,sum(it$rulesused))
	numberdlty[i]<-ifelse(is.na(it$dose),NA,it$numberdlty)
	dltruleused[i]<-it$rulesused[1]
	incrruleused[i]<-it$rulesused[2]
	lowruleused[i]<-it$rulesused[3]
	upruleused[i]<-it$rulesused[4]

#within 10 or 20%:
if10percent[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]<=(truedose+(truedose*.1)) & sim[i]>=(truedose-(truedose*.1)),1,0))
if20percent[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]<=(truedose+(truedose*.2)) & sim[i]>=(truedose-(truedose*.2)),1,0))

# % >40% DLT and >50%DLT
if40[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]>truedlt40,1,0))
if50[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]>truedlt50,1,0))

# % <10% DLT and <20%DLT
if10[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]<truedlt10,1,0))
if20[i]<-ifelse(is.na(it$dose),NA,ifelse(sim[i]<truedlt20,1,0))

print(i)
	}

##############################
#RESULTS OVER ALL SIMULATIONS#
##############################

#IF DISCRETE...PERCENTAGE OF TIMES SIMULATIONS PICKED EACH DOSE LEVEL#
if (discrete==TRUE)	{
percentagedose<-vector("list",length(discretedoses))
percentagepickeddiscretedose<-c()
	for (i in 1:length(discretedoses))	{
	percentagedose[[i]]<-which(sim==discretedoses[i])
	percentagepickeddiscretedose[i]<-round(((length(percentagedose[[i]]))/numbersims)*100,2)
	}
absdifftruedose<-c()
	for (i in 1:length(discretedoses))	{
			absdifftruedose[i]<-abs(truedose-discretedoses[i])
	}
	tt<-order(absdifftruedose, decreasing=FALSE)
	ordernextdose<-cbind(absdifftruedose,discretedoses)[tt,]
#	if (rounddown==TRUE)	{
#	truediscdose<-ifelse(ordernextdose[1,2]>truedose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
#	}
#	if (rounddown==FALSE)	{
#	truediscdose<-ordernextdose[1,2]
#	}
	truediscdose<-ordernextdose[1,2]
}

if (discrete==FALSE)	{
percentagepickeddiscretedose<-NA
truediscdose<-NA}

missing<-which(is.na(sim))
nummissing<-length(missing)
	propmissing<-round(100*(length(missing)/numbersims),digits=2)
totalobs<-numbersims-length(missing)
nuconstraintused<-length(which(constraintused==1))
	propconstraintused<-round(100*(nuconstraintused/totalobs),2)
nudltruleused<-length(which(dltruleused==1))
	propdltruleused<-round(100*(nudltruleused/totalobs),2)
nuincrruleused<-length(which(incrruleused==1))
	propincrruleused<-round(100*(nuincrruleused/totalobs),2)
nulowruleused<-length(which(lowruleused==1))
	proplowruleused<-round(100*(nulowruleused/totalobs),2)
nuupruleused<-length(which(upruleused==1))
	propupruleused<-round(100*(nuupruleused/totalobs),2)
###Median of Trials that Completed to Total SampleSize#
mediansize<-stats::median(Nna,na.rm=TRUE)

###Of Trials with an Estimated MTD###
trials10percent<-round(100*(sum(stats::na.omit(if10percent))/totalobs),digits=2)
trials20percent<-round(100*(sum(stats::na.omit(if20percent))/totalobs),digits=2)
trials40<-round(100*(sum(stats::na.omit(if40))/totalobs),digits=2)
trials50<-round(100*(sum(stats::na.omit(if50))/totalobs),digits=2)
trials20<-round(100*(sum(stats::na.omit(if20))/totalobs),digits=2)
trials10<-round(100*(sum(stats::na.omit(if10))/totalobs),digits=2)

###This is results for trials that did not stop early###
people40<-round(stats::median(over40.pct,na.rm=TRUE),digits=2)
people50<-round(stats::median(over50.pct,na.rm=TRUE),digits=2)
people20<-round(stats::median(under20.pct,na.rm=TRUE),digits=2)
people10<-round(stats::median(under10.pct,na.rm=TRUE),digits=2)
mediannumdlt<-round(stats::median(dlts.n, na.rm=TRUE),digits=2)
median3or4<-round(stats::median(dlts.pct,na.rm=TRUE),digits=2)
median1or2<-ifelse(combine01==TRUE,NA,round(stats::median(ts1.pct+ts2.pct,na.rm=TRUE),digits=2))
median0<-ifelse(combine01==TRUE,NA,stats::median(round(100-(ts1.pct+ts2.pct+
	ts3.pct+ts4.pct),digits=2),na.rm=TRUE))
median1<-ifelse(combine01==TRUE,NA,round(stats::median(ts1.pct,na.rm=TRUE),digits=2))
median2<-round(stats::median(ts2.pct,na.rm=TRUE),digits=2)
median3<-round(stats::median(ts3.pct,na.rm=TRUE),digits=2)
median4<-round(stats::median(ts4.pct,na.rm=TRUE),digits=2)
median01<-ifelse(combine01==TRUE,round(stats::median(ts01.pct,na.rm=TRUE),digits=2),NA)

###All Trials Even If They Stopped Early###
allmediansize<-round(stats::median(N),digits=2)
allpeople40<-round(stats::median(allover40.pct),digits=2)
allpeople50<-round(stats::median(allover50.pct),digits=2)
allpeople20<-round(stats::median(allunder20.pct),digits=2)
allpeople10<-round(stats::median(allunder10.pct),digits=2)
allmediannumdlt<-round(stats::median(alldlts.n),digits=2)
allmedian3or4<-round(stats::median(alldlts.pct),digits=2)
allmedian1or2<-ifelse(combine01==TRUE, NA, round(stats::median(allts1.pct+allts2.pct),digits=2))
allmedian0<-ifelse(combine01==TRUE,NA,stats::median(round(100-(allts1.pct+allts2.pct+
	allts3.pct+allts4.pct),digits=2)))
allmedian1<-ifelse(combine01==TRUE,NA,round(stats::median(allts1.pct),digits=2))
allmedian2<-round(stats::median(allts2.pct),digits=2)
allmedian3<-round(stats::median(allts3.pct),digits=2)
allmedian4<-round(stats::median(allts4.pct),digits=2)
allmedian01<-ifelse(combine01==TRUE,round(stats::median(allts01.pct),digits=2), NA)

stat<-summary(sim)
N=N
numalreadytested2=numalreadytested2
Nna=Nna
sim=sim
dosequantile<-round(stats::quantile(sim, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3),digits=0)
dltnquantile<-round(stats::quantile(dlts.n,prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3),digits=0)
dltquantile<-round(stats::quantile(dlts.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3),digits=2)
over40quantile<-round(stats::quantile(over40.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3),digits=2)
under20quantile<-round(stats::quantile(under20.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3),digits=2)

alldltnquantile<-stats::quantile(alldlts.n, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
alldltquantile<-stats::quantile(alldlts.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
allover40quantile<-stats::quantile(allover40.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
allunder20quantile<-stats::quantile(allunder20.pct, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
perchangefromMTD<-round((((sim-truedose)/truedose)*100),digits=2)
perchangefromMTDquantile<-stats::quantile(perchangefromMTD,prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
expectedDLTquantile<-stats::quantile(expectedDLT,prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)
samplesizequantile<-stats::quantile(N, prob=c(0.05,0.25,0.75,0.95),na.rm=TRUE,type=3)

	return(list('Median Total Sample Size of Trials Not Stopped Early Due to Safety Concerns'=mediansize,
		'Patients per Cohort'=cohortsize,
		'Proportion of Trials Stopped Early Due to Safety Concerns'=propmissing,
		'Proportion of Trials Stopped Early Used Constraint'=propconstraintused,
		'Proportion of Trials Stopped Early Used the DLT Constraint'=propdltruleused,
		'Proportion of Trials Stopped Early Used the Increase Constraint'=propincrruleused,
		'25% Quantile Sample Size'=samplesizequantile['25%'],
		'75% Quantile Sample Size'=samplesizequantile['75%'],
		'Median Dose'=stats::median(sim,na.rm=TRUE),
		'25% Quantile Dose'=dosequantile['25%'],
		'75% Quantile Dose'=dosequantile['75%'],
		'Percentage of Times Each Discrete Dose Level Was Selected'=percentagepickeddiscretedose,
		'Median Percent Difference Between Estimated Dose and MTD'=stats::median(perchangefromMTD,na.rm=TRUE),
		'Median Expected DLT for the Final Estimated Dose'=stats::median(expectedDLT,na.rm=TRUE),
		'Median Sample Size For All Trials'=allmediansize,
		'Median % of patients in all trials treated at doses with >40% DLT rate'=allpeople40,
		'Median % of patients in all trials treated at doses with >50% DLT rate'=allpeople50,
		'Median % of patients in all trials treated at doses with <20% DLT rate'=allpeople20,
		'Median % of patients in all trials treated at doses with <10% DLT rate'=allpeople10,
		'Median Number of patients in all trials with a DLT (grade 3 or 4)'=allmediannumdlt,
		'Median % of patients in all trials with a DLT (grade 3 or 4)'=allmedian3or4,
		'Median % of patients in all trials with grade 0 and 1'=allmedian01,
		'Median % of patients in all trials with a non-DLT (grade 1 or 2)'=allmedian1or2,
		'Median % of patients in all trials with no toxicity'=allmedian0,
		'Median % of patients in all trials with grade 1'=allmedian1,
		'Median % of patients in all trials with grade 2'=allmedian2,
		'Median % of patients in all trials with grade 3'=allmedian3,
		'Median % of patients in all trials with grade 4'=allmedian4,
		'Acutal MTD'=truedose,
		'True MTD Discrete'=truediscdose,
		'Dose'=sim,
		'Listit'=listit

))
}

