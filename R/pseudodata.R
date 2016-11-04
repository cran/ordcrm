#' Pseudodata for Likelihood-Based Continual Reassessment Method (CRM) Dose Finding Designs
#'
#' This function creates 2 pseudodata vectors (toxicity and dose levels) necessary to run any of the
#' 3 variations of the likelihood-based CRM design. If using the ordinal POM or CR model design,
#' ordinal toxicity grades are specified by CTCAEv4.0. This function also specifies the starting dose
#' for a trial with the choice of continuous or discrete dose levels for a specified target dose limiting
#' toxicity (DLT) rate.
#' @param design Specifies which dose finding design you are running simulations on.
#' Choices are: POM, CR, or CRM.
#' @param dose10 Hypothesized dose level for a 10 percent dose limiting toxicity (DLT) rate as specified
#' by clinical investigators prior to the start of the trial.
#' @param dose90 Hypothesized dose level for a 90 percent dose limiting toxicity (DLT) rate as specified
#' by clinical investigators prior to the start of the trial.
#' @param percentagegradedose10 Only applicable if using ordinal POM or CR model designs.
#' Expected percentage of toxicity grades at the specified dose10 in the form
#' c(percent grade 0, percent grade 1, grade 2, grade 3, grade 4) or
#' c(percent grades 0 and 1, grade 2, grade 3, grade 4) if combining grades 0 and 1 into one category.
#' The percentages must be specified between 0 and 100. If not specified,
#' it will default to c(45, 35, 10, 8, 2).
#' @param percentagegradedose90 Only applicable if using ordinal POM or CR model designs.
#' Expected percentage of toxicity grades at the specified dose90. This is also in the form:
#' c(percent grade 0, percent grade 1, grade 2, grade 3, grade 4) or
#' c(percent grades 0 and 1, grade 2, grade 3, grade 4) if combining grades 0 and 1 into one category.
#' The percentages must be specified between 0 and 100. If not specified,
#' it will default to c(2, 3, 5, 40, 50).
#' @param targetDLT Target dose limiting toxicity (DLT) rate pre-specified by clinical investigators
#' prior to the start of the trial. Must be specified between 0 and 1.
#' Defaults to 0.30.
#' @param stabilize True/False. If stabilize = TRUE, expected dose levels for 30 percent
#' and 50 percent DLT rates will be estimated from the initialized CR model and these levels will be
#' added into the pseudodata to help stabilize the model. This is particularly helpful when
#' collected data is initially sparse at the beginning of a trial. Defaults to TRUE.
#' @param discrete True/False. If discrete = TRUE, this allows for discrete dose levels to be
#' specified prior to the start of the trial. Defaults to FALSE.
#' @param discretedoses Specified discrete dose levels if desired. They must be specified if discrete is equal to TRUE.
#' It is written for j dose levels as: c(d1, d2,..., dj).
#' @param rounddown True/False. Only applicable when using discrete dose levels.\cr
#' If rounddown = TRUE, the estimate dose from specified model will round down to the more conservative discrete dose level.\cr
#' If rounddown = FALSE, it will select the discrete dose closest to the estimated model selection. Defaults to FALSE.
#' @param combine01 True/False. If combine01 = TRUE, toxicity grades 0 and 1 are combined into 1 category.
#' Therefore all toxicities must be coded: 0 (grades 0 and 1), 1 (grade2), 2 (grade 3),
#' and 3 (grade 4) according to CTCAEv4.0. Defaults to FALSE.
#' @param plotit True/False. If plotit = TRUE, returns a plot of the continuation ratio model estimated by
#' the pseudodata with the starting dose for a specified target DLT rate.
#' Also identifies dose10 and dose90 on the figure. Defaults to TRUE.
#'
#' @details
#' Fits a POM, CR model, or 2-parameter logistic model based on the dose where 10% and 90% DLTs are hypothesized to occur.
#' This outputs a vector of doses and toxicities to use as the starting model for this design (required to run
#' \code{\link{nextdose}} and \code{\link{crmsimulations}})
#'
#' @return \item{Pseudodata Toxicities}{A vector of toxicity grades as classified by CTCAEv4.0
#' created from the initialized model. If using a binary CRM this will be a vector of 0 (non-DLT) and
#' 1 (DLT) outcomes. This will be the pseudodata toxicity grades used for the rest of trial.}
#' @return \item{Pseudodata Doses}{A vector of dose levels created from the initialized model.
#' This will be the pseudodata dose levels used for the rest of the trial.}
#' @return \item{Starting Dose}{The starting dose for this trial given the specified target DLT rate assuming
#' continuous dose levels.}
#' @return \item{Starting Discrete Dose}{The starting dose for this trial given the specified
#' target DLT rate assuming discrete dose levels.}
#' @return \item{Estimated Toxicity at Discrete Doses}{Only applicable if using specified discrete dose levels.
#' Returns the expected probability of
#' a dose-limiting toxicity (DLT) at each discrete dose level given the initialized model.}
#' @return \item{discretedoses}{Returns the possible ordered discrete doses if applicable.}
#' @return \item{Regression Model}{Returns parameter estimates for the newly estimated model.}
#'
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
#' @seealso \code{\link{nextdose}}
#' @examples
#'
#' ######################
#' ###For POM Model######
#' ######################
#' #Creates pseudodata with no stabilization and continuous dose levels
#' pseudodata(design='POM', dose10 = 200, dose90 = 3600, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = FALSE, discretedoses = NA,
#'            percentagegradedose10 = c(45, 35, 10, 8, 2),
#'            percentagegradedose90 = c(2, 3, 5, 40, 50), combine01 = FALSE)
#' #Creates pseudodata with stabilize=T and discrete dose levels
#' pseudodata(design='POM', dose10 = 500, dose90 = 2000, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = TRUE,
#'            discretedoses = c(200, 500, 1000, 1200, 1500, 1800))
#' #Creates pseudodata when toxicity grades 0 and 1 are one category
#' pseudodata(design='POM', dose10 = 200, dose90 = 1000, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = FALSE, discretedoses = NA,
#'            percentagegradedose10 = c(80, 10, 8, 2),
#'            percentagegradedose90 = c(5, 5, 40, 50), combine01 = TRUE)
#'
#' ######################
#' ###For CR Model#######
#' ######################
#' #Creates pseudodata with no stabilization and continuous dose levels
#' pseudodata(design='CR', dose10 = 200, dose90 = 3600, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = FALSE, discretedoses = NA,
#'            percentagegradedose10 = c(45, 35, 10, 8, 2),
#'            percentagegradedose90 = c(2, 3, 5, 40, 50), combine01 = FALSE)
#' #Creates pseudodata with stabilize=T and discrete dose levels
#' pseudodata(design='CR', dose10 = 500, dose90 = 2000, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = TRUE,
#'            discretedoses = c(200, 500, 1000, 1200, 1500, 1800))
#' #Creates pseudodata when toxicity grades 0 and 1 are one category
#' pseudodata(design='CR', dose10 = 200, dose90 = 1000, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = FALSE, discretedoses = NA,
#'            percentagegradedose10 = c(80, 10, 8, 2),
#'            percentagegradedose90=c(5, 5, 40, 50), combine01 = TRUE)
#'
#' ######################
#' ###For Binary CRM#####
#' ######################
#' #Creates pseudodata with no stabilization and continuous dose levels
#' pseudodata(design='CRM', dose10 = 200, dose90 = 3600, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = FALSE, discretedoses = NA)
#' #Creates pseudodata with stabilize=T and discrete dose levels
#' pseudodata(design='CRM', dose10 = 200, dose90 = 3000, targetDLT = 0.3,
#'            stabilize = TRUE, discrete = TRUE,
#'            discretedoses = c(200, 500, 1000, 1200, 1500, 1800))




pseudodata <-
  function(design, dose10, dose90, percentagegradedose10=c(45,35,10,8,2),
           percentagegradedose90=c(2,3,5,40,50), targetDLT=0.30, stabilize=TRUE,
           discrete=FALSE, discretedoses=NA, rounddown=FALSE, combine01=FALSE, plotit=TRUE)	{

    #WARNINGS#
    if ((design!='POM' & design!='CR' & design!='CRM'))	{
      stop('You must specify the dose finding design you wish to use.
           Choices are POM (proportional odds model ordinal design), CRM (original binary CRM),
           or CR (continuation ratio model ordinal design)')
    }
    if ((discrete==FALSE) & rounddown==TRUE)	{
      stop('Rounding down only applied to discrete dose levels')
    }
    if ((combine01==TRUE) & ((length(percentagegradedose10)!=4) | length(percentagegradedose10)!=4)) {
      stop('If combining grades 0 and 1 into one category,
           percentagegradedose10 and percentagegradedose90
           must be in the form c(% grade 0/1, % grade 2, % grade 3, % grade 4)')
    }
    if ((discrete==TRUE) & is.na(discretedoses[1])==TRUE)	{
      stop('If using discrete doses you must specify dose levels to run this function')
    }
    if ((is.na(discretedoses[1])==FALSE) & discrete==FALSE)	{
      stop('If specifying discrete dose levels, discrete must equal TRUE')
    }
    if ((combine01==FALSE) & ((length(percentagegradedose10)==4) | length(percentagegradedose10)==4)) {
      stop('You specified percentagegradedose10 and percentagegradedose90
           in the form c(% grade 0/1, % grade 2, % grade 3, % grade 4).
           In this form, did you mean to combine grades 0 and 1 into one category?')
    }
    if (design=='CRM')	{
      percentagegradedose10<-percentagegradedose90<-NA}

    #########################################################
    ###POM Pseudodata #######################################
    #########################################################
    if (design=='POM')	{
      if (combine01==FALSE)	{
        d<-c(rep(dose10,100),rep(dose90,100))
        y<-c(rep(c(0,1,2,3,4),percentagegradedose10),rep(c(0,1,2,3,4),percentagegradedose90))
        reg<-rms::lrm(y~d)

        extraDLT<-c()
        extraDLT[1]<-ifelse(stabilize==TRUE,.30,NA)
        extraDLT[2]<-ifelse(stabilize==TRUE,.50,NA)
        if (is.na(extraDLT[1])==FALSE)	{
          extra<-vector("list",2)
          for (i in 1:2)	{
            dose<-round((log(extraDLT[i]/(1-extraDLT[i]))-reg$coef[3])/reg$coef['d'])
            phat1ormore<- reg$coef[1]+ (reg$coef['d']*dose)
            phat1ormore<-1/(1+exp(-phat1ormore))
            phat2ormore<- reg$coef[2]+ (reg$coef['d']*dose)
            phat2ormore<-1/(1+exp(-phat2ormore))
            phat3ormore<- reg$coef[3]+ (reg$coef['d']*dose)
            phat3ormore<-1/(1+exp(-phat3ormore))
            phat4ormore<- reg$coef[4]+ (reg$coef['d']*dose)
            phat4ormore<-1/(1+exp(-phat4ormore))
            phat4<-phat4ormore
            phat3<-phat3ormore-phat4ormore
            phat2<-phat2ormore-phat3ormore
            phat1<-phat1ormore-phat2ormore
            phat0<-1-phat1ormore

            extra[[i]][1]<-dose
            extra[[i]][2]<-phat0
            extra[[i]][3]<-phat1
            extra[[i]][4]<-phat2
            extra[[i]][5]<-phat3
            extra[[i]][6]<-phat4
          }

          dose30<-round(extra[[1]],2)
          dose50<-round(extra[[2]],2)
          diff30a<-.70-(dose30[2]+dose30[3]+dose30[4])
          dose30[4]<-dose30[4]+diff30a
          diff30b<-.30-(dose30[5]+dose30[6])
          dose30[5]<-dose30[5]+diff30b
          diff50a<-.50-(dose50[2]+dose50[3]+dose50[4])
          dose50[4]<-dose50[4]+diff50a
          diff50b<-.50-(dose50[5]+dose50[6])
          dose50[5]<-dose50[5]+diff50b
          percent30<-round(dose30[2:6]*100,0)
          percent50<-round(dose50[2:6]*100,0)

          pseudoy<-append(y,c(rep(c(0,1,2,3,4),percent30),rep(c(0,1,2,3,4),percent50)),after=length(y))
          pseudod<-append(d,c(rep(round(dose30[1],0),100),rep(round(dose50[1],0),100)),after=length(d))
        }
      }

      if (combine01==TRUE)	{
        d<-c(rep(dose10,100),rep(dose90,100))
        y<-c(rep(c(0,1,2,3),percentagegradedose10),rep(c(0,1,2,3),percentagegradedose90))
        reg<-rms::lrm(y~d)

        extraDLT<-c()
        extraDLT[1]<-ifelse(stabilize==TRUE,.30,NA)
        extraDLT[2]<-ifelse(stabilize==TRUE,.50,NA)
        if (is.na(extraDLT[1])==FALSE)	{
          extra<-vector("list",2)
          for (i in 1:2)	{
            dose<-round((log(extraDLT[i]/(1-extraDLT[i]))-reg$coef[2])/reg$coef['d'])
            phat2ormore<- reg$coef[1]+ (reg$coef['d']*dose)
            phat2ormore<-1/(1+exp(-phat2ormore))
            phat3ormore<- reg$coef[2]+ (reg$coef['d']*dose)
            phat3ormore<-1/(1+exp(-phat3ormore))
            phat4ormore<- reg$coef[3]+ (reg$coef['d']*dose)
            phat4ormore<-1/(1+exp(-phat4ormore))
            phat4<-phat4ormore
            phat3<-phat3ormore-phat4ormore
            phat2<-phat2ormore-phat3ormore
            phat0or1<-1-phat2ormore
            extra[[i]][1]<-dose
            extra[[i]][2]<-phat0or1 #prob 0 or 1#
            extra[[i]][3]<-phat2 #prob 2#
            extra[[i]][4]<-phat3 #prob 3#
            extra[[i]][5]<-phat4 #prob 4#
          }

          dose30<-round(extra[[1]],2)
          dose50<-round(extra[[2]],2)
          diff30a<-.70-(dose30[2]+dose30[3])
          dose30[3]<-dose30[3]+diff30a
          diff30b<-.30-(dose30[4]+dose30[5])
          dose30[4]<-dose30[4]+diff30b
          diff50a<-.50-(dose50[2]+dose50[3])
          dose50[3]<-dose50[3]+diff50a
          diff50b<-.50-(dose50[4]+dose50[5])
          dose50[4]<-dose50[4]+diff50b
          percent30<-round(dose30[2:5]*100,0)
          percent50<-round(dose50[2:5]*100,0)

          pseudoy<-append(y,c(rep(c(0,1,2,3),percent30),rep(c(0,1,2,3),percent50)),after=length(y))
          pseudod<-append(d,c(rep(round(dose30[1],0),100),rep(round(dose50[1],0),100)),after=length(d))
        }
      }

      if (is.na(extraDLT[1]))	{
        pseudoy<-y
        pseudod<-d
      }

      fullmodel<-rms::lrm(pseudoy~pseudod)
      startdose<-ifelse(combine01==FALSE,
                        round((log(targetDLT/(1-targetDLT))-fullmodel$coef[3])/fullmodel$coef['pseudod']),
                        round((log(targetDLT/(1-targetDLT))-fullmodel$coef[2])/fullmodel$coef['pseudod']))
      if (discrete==TRUE)	{
        discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
        esttoxdiscdose<-diffdose<-absdiffdose<-c()
        for (i in 1:length(discretedoses))	{
          esttoxdiscdose[i]<-ifelse(combine01==FALSE, 1/(1+exp(-1*(fullmodel$coef[3]+ (fullmodel$coef['pseudod']*discretedoses[i])))),
                                    1/(1+exp(-1*(fullmodel$coef[2]+ (fullmodel$coef['pseudod']*discretedoses[i])))))
          diffdose[i]<-(startdose-discretedoses[i])
          absdiffdose[i]<-abs(startdose-discretedoses[i])
        }
        tt<-order(absdiffdose, decreasing=FALSE)
        ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
        if (rounddown==TRUE)	{
          discdose<-ifelse(ordernextdose[1,2]>startdose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
        }
        if (rounddown==FALSE)	{
          discdose<-ordernextdose[1,2]
        }
      }
      if (discrete==FALSE)	{
        esttoxdiscdose<-NA
        discdose<-NA
      }

      if (plotit==TRUE)	{
        if (discrete==FALSE)	{
          if (combine01==FALSE)	{
            graphics::plot(c(0,round(dose90*1.05,0)),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp1<- fullmodel$coef[1]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            tmp2<- fullmodel$coef[2]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            tmp3<- fullmodel$coef[3]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            tmp4<- fullmodel$coef[4]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp1)),lty=1, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            graphics::points(startdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x=(startdose*1.05), y=(targetDLT+.03), labels=startdose)
            graphics::legend(0,1,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Starting Dose","Dose 10 & Dose 90"),
                   lty=c(1,2,3,4,-1,-1),pch=c(-1,-1,-1,-1,15,1),lwd=c(2,2,2,2,1,1),
                   col=c(1,1,1,1,"red",1),bg="white")
            graphics::title("Proportional Odds Model CRM Pseudodata")
          }
          if (combine01==TRUE)	{
            graphics::plot(c(0,round(dose90*1.05,0)),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp2<- fullmodel$coef[1]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            tmp3<- fullmodel$coef[2]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            tmp4<- fullmodel$coef[3]+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            graphics::points(startdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x=(startdose*1.05), y=(targetDLT+.03), labels=startdose)
            graphics::legend(0,1,c("p(y>=2)","p(y>=3)","p(y>=4)","Starting Dose","Dose 10 & Dose 90"),
                   lty=c(2,3,4,-1,-1),pch=c(-1,-1,-1,15,1),lwd=c(2,2,2,1,1),
                   col=c(1,1,1,"red",1),bg="white")
            graphics::title("Proportional Odds Model CRM Pseudodata")
          }
        }
        if (discrete==TRUE)	{
          if (combine01==FALSE)	{
            lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
            graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0)))),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp1<- fullmodel$coef[1]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp2<- fullmodel$coef[2]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp3<- fullmodel$coef[3]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp4<- fullmodel$coef[4]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp1)),lty=1, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            for (i in 1:length(discretedoses))	{
              graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
              graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
            graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
            graphics::legend(0,1,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Starting Discrete Dose","Dose 10 & Dose 90", "Estimated Toxicity"),
                   lty=c(1,2,3,4,-1,-1,-1),pch=c(-1,-1,-1,-1,15,1,15),lwd=c(2,2,2,2,1,1,1),
                   col=c(1,1,1,1,"red",1,"gray90"),bg="white")
            graphics::title("Proportional Odds Model CRM Pseudodata")
          }
          if (combine01==TRUE)	{
            lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
            graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0)))),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp2<- fullmodel$coef[1]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp3<- fullmodel$coef[2]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp4<- fullmodel$coef[3]+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp2)),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp3)),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-tmp4)),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            for (i in 1:length(discretedoses))	{
              graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
              graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
            graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
            graphics::legend(0,1,c("p(y>=2)","p(y>=3)","p(y>=4)","Starting Discrete Dose","Dose 10 & Dose 90", "Estimated Toxicity"),
                   lty=c(2,3,4,-1,-1,-1),pch=c(-1,-1,-1,15,1,15),lwd=c(2,2,2,1,1,1),
                   col=c(1,1,1,"red",1,"gray90"),bg="white")
            graphics::title("Proportional Odds Model CRM Pseudodata")
          }
        }
      }
    }

    #############################################
    ###Continuation Ratio Model #################
    #############################################
    if (design=='CR')	{
      if (combine01==FALSE)	{
        dall<-c(rep(dose10,100),rep(dose90,100))
        yall<-c(rep(c(0,1,2,3,4),percentagegradedose10),rep(c(0,1,2,3,4),percentagegradedose90))
        y1<-rms::cr.setup(yall)
        dose<-dall[y1$subs]
        y<-y1$y
        cohort<-y1$cohort
        reg<-rms::lrm(y ~cohort+dose)
        fullmodel<-reg

        extraDLT<-c()
        extraDLT[1]<-ifelse(stabilize==TRUE,.30,NA)
        extraDLT[2]<-ifelse(stabilize==TRUE,.50,NA)
        if (is.na(extraDLT[1])==FALSE)	{
          extra<-vector("list",2)
          for (i in 1:2)	{
            d<-(-((1-extraDLT[i])/(extraDLT[i])))
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
            possibledoses<-ifelse((dcheck>=extraDLT[i]-0.01 & dcheck<=extraDLT[i]+0.01),blah,NA)
            d.new<-blah[which(is.na(possibledoses)==FALSE)]
            dosehat<-ifelse(length(d.new)==1, round(d.new,0), NA)

            crequal0<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*dosehat)))))
            crequal1<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*dosehat)))))
            crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*dosehat)))))
            crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*dosehat)))))

            phatmorethan0<-1-crequal0
            phatmorethan1<-(1-crequal1)*phatmorethan0
            phatmorethan2<-(1-crequal2)*phatmorethan1
            phatmorethan3<-(1-crequal3)*phatmorethan2
            phat4<-phatmorethan3
            phat3<-phatmorethan2-phatmorethan3
            phat2<-phatmorethan1-phatmorethan2
            phat1<-phatmorethan0-phatmorethan1
            phat0<-1-phatmorethan0
            extra[[i]][1]<-dosehat
            extra[[i]][2]<-phat0
            extra[[i]][3]<-phat1
            extra[[i]][4]<-phat2
            extra[[i]][5]<-phat3
            extra[[i]][6]<-phat4
          }

          dose30<-round(extra[[1]],2)
          dose50<-round(extra[[2]],2)
          diff30a<-.70-(dose30[2]+dose30[3]+dose30[4])
          dose30[4]<-dose30[4]+diff30a
          diff30b<-.30-(dose30[5]+dose30[6])
          dose30[5]<-dose30[5]+diff30b
          diff50a<-.50-(dose50[2]+dose50[3]+dose50[4])
          dose50[4]<-dose50[4]+diff50a
          diff50b<-.50-(dose50[5]+dose50[6])
          dose50[5]<-dose50[5]+diff50b
          percent30<-round(dose30[2:6]*100,0)
          percent50<-round(dose50[2:6]*100,0)

          pseudoy<-append(yall,c(rep(c(0,1,2,3,4),percent30),rep(c(0,1,2,3,4),percent50)),after=length(yall))
          pseudod<-append(dall,c(rep(round(dose30[1],0),100),rep(round(dose50[1],0),100)),after=length(dall))
        }
      }

      if (combine01==TRUE)	{
        dall<-c(rep(dose10,100),rep(dose90,100))
        yall<-c(rep(c(0,1,2,3),percentagegradedose10),rep(c(0,1,2,3),percentagegradedose90))
        y1<-rms::cr.setup(yall)
        dose<-dall[y1$subs]
        y<-y1$y
        cohort<-y1$cohort
        reg<-rms::lrm(y ~cohort+dose)
        fullmodel<-reg

        extraDLT<-c()
        extraDLT[1]<-ifelse(stabilize==TRUE,.30,NA)
        extraDLT[2]<-ifelse(stabilize==TRUE,.50,NA)
        if (is.na(extraDLT[1])==FALSE)	{
          extra<-vector("list",2)
          for (i in 1:2)	{
            c<-(-((1-extraDLT[i])/(extraDLT[i])))
            b<-(exp(reg$coef['Intercept'])+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']))
            a<-(exp((2*reg$coef['Intercept'])+reg$coef['cohort=yall>=1']))
            z<-c(c,b,a)
            possibledoses<-polyroot(z)
            blah<-(log(possibledoses)/reg$coef['dose'])
            blah<-suppressWarnings(as.double(blah))
            dcheck<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*blah)))*
                         (1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*blah))))
            possibledoses<-ifelse((dcheck>=extraDLT[i]-0.01 & dcheck<=extraDLT[i]+0.01),blah,NA)
            d.new<-blah[which(is.na(possibledoses)==FALSE)]
            dosehat<-(ifelse(length(d.new)==1, round(d.new,0), NA))

            crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*dosehat)))))
            crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*dosehat)))))
            crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*dosehat)))))

            phatmorethan01<-1-crequal01
            phatmorethan2<-(1-crequal2)*phatmorethan01
            phatmorethan3<-(1-crequal3)*phatmorethan2
            phat4<-phatmorethan3
            phat3<-phatmorethan2-phatmorethan3
            phat2<-phatmorethan01-phatmorethan2
            phat01<-1-phatmorethan01
            extra[[i]][1]<-dosehat
            extra[[i]][2]<-phat01
            extra[[i]][3]<-phat2
            extra[[i]][4]<-phat3
            extra[[i]][5]<-phat4
          }

          dose30<-round(extra[[1]],2)
          dose50<-round(extra[[2]],2)
          diff30a<-.70-(dose30[2]+dose30[3])
          dose30[3]<-dose30[3]+diff30a
          diff30b<-.30-(dose30[4]+dose30[5])
          dose30[4]<-dose30[4]+diff30b
          diff50a<-.50-(dose50[2]+dose50[3])
          dose50[3]<-dose50[3]+diff50a
          diff50b<-.50-(dose50[4]+dose50[5])
          dose50[4]<-dose50[4]+diff50b
          percent30<-round(dose30[2:5]*100,0)
          percent50<-round(dose50[2:5]*100,0)

          pseudoy<-append(yall,c(rep(c(0,1,2,3),percent30),rep(c(0,1,2,3),percent50)),after=length(yall))
          pseudod<-append(dall,c(rep(round(dose30[1],0),100),rep(round(dose50[1],0),100)),after=length(dall))
        }
      }

      if (is.na(extraDLT[1]))	{
        pseudoy<-yall
        pseudod<-dall
      }

      #IDENTIFYING STARTING DOSE#
      if (combine01==FALSE)	{
        yall<-pseudoy
        dall<-pseudod
        y1<-rms::cr.setup(yall)
        dose<-dall[y1$subs]
        y<-y1$y
        cohort<-y1$cohort
        reg<-rms::lrm(y ~cohort+dose)
        fullmodel<-reg

        d<-(-((1-targetDLT)/(targetDLT)))
        c<-(exp(reg$coef['Intercept'])+exp(reg$coef['Intercept']+reg$coef['y1=yall>=1'])+
              exp(reg$coef['Intercept']+reg$coef['y1=yall>=2']))
        b<-(exp((2*reg$coef['Intercept'])+reg$coef['y1=yall>=1'])+exp((2*reg$coef['Intercept'])+reg$coef['y1=yall>=2'])+
              exp((2*reg$coef['Intercept'])+reg$coef['y1=yall>=1']+reg$coef['y1=yall>=2']))
        a<-(exp((3*reg$coef['Intercept'])+reg$coef['y1=yall>=1']+reg$coef['y1=yall>=2']))
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
                     (1+exp(reg$coef['Intercept']+reg$coef['y1=yall>=1']+(reg$coef['dose']*blah)))*
                     (1+exp(reg$coef['Intercept']+reg$coef['y1=yall>=2']+(reg$coef['dose']*blah))))
        possibledoses<-ifelse((dcheck>=targetDLT-0.01 & dcheck<=targetDLT+0.01),blah,NA)
        d.new<-blah[which(is.na(possibledoses)==FALSE)]
        startdose<-ifelse(length(d.new)==1, round(d.new,0), NA)
        if (discrete==TRUE)	{
          discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
          esttoxdiscdose<-diffdose<-absdiffdose<-c()
          for (i in 1:length(discretedoses))	{
            esttoxdiscdose[i]<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*discretedoses[i])))*
                                    (1+exp(reg$coef['Intercept']+reg$coef['y1=yall>=1']+(reg$coef['dose']*discretedoses[i])))*
                                    (1+exp(reg$coef['Intercept']+reg$coef['y1=yall>=2']+(reg$coef['dose']*discretedoses[i]))))
            diffdose[i]<-(startdose-discretedoses[i])
            absdiffdose[i]<-abs(startdose-discretedoses[i])
          }
          tt<-order(absdiffdose, decreasing=FALSE)
          ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
          if (rounddown==TRUE)	{
            discdose<-ifelse(ordernextdose[1,2]>startdose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
          }
          if (rounddown==FALSE)	{
            discdose<-ordernextdose[1,2]
          }
        }
        if (discrete==FALSE)	{
          esttoxdiscdose<-NA
          discdose<-NA
        }
      }

      if (combine01==TRUE)	{
        yall<-pseudoy
        dall<-pseudod
        y1<-rms::cr.setup(yall)
        dose<-dall[y1$subs]
        y<-y1$y
        cohort<-y1$cohort
        reg<-rms::lrm(y ~cohort+dose)
        fullmodel<-reg

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
        startdose<-ifelse(length(d.new)==1, round(d.new,0), NA)
        if (discrete==TRUE)	{
          discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
          esttoxdiscdose<-diffdose<-absdiffdose<-c()
          for (i in 1:length(discretedoses))	{
            esttoxdiscdose[i]<-dcheck<-1/((1+exp(reg$coef['Intercept']+(reg$coef['dose']*discretedoses[i])))*
                                            (1+exp(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*discretedoses[i]))))
            diffdose[i]<-(startdose-discretedoses[i])
            absdiffdose[i]<-abs(startdose-discretedoses[i])
          }
          tt<-order(absdiffdose, decreasing=FALSE)
          ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
          if (rounddown==TRUE)	{
            discdose<-ifelse(ordernextdose[1,2]>startdose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
          }
          if (rounddown==FALSE)	{
            discdose<-ordernextdose[1,2]
          }
        }
        if (discrete==FALSE)	{
          esttoxdiscdose<-NA
          discdose<-NA
        }
      }

      ########################################
      #Plotting Pseudodata####################
      ########################################
      if (plotit==TRUE)	{
        if (discrete==FALSE)	{
          if (combine01==FALSE)	{
            graphics::plot(c(0,round(dose90*1.05,0)),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp0<- reg$coef['Intercept']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))
            tmp1<- reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))
            tmp2<- reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))
            tmp3<- reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1/(1+exp(tmp0))),lty=1, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1/(1+exp(tmp3)))*(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            graphics::points(startdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x=(startdose*1.05), y=(targetDLT+.03), labels=startdose)
            graphics::legend(0,1,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Starting Dose","Dose 10 & Dose 90"),
                   lty=c(1,2,3,4,-1,-1),pch=c(-1,-1,-1,-1,15,1),lwd=c(2,2,2,2,1,1),
                   col=c(1,1,1,1,"red",1),bg="white")
            graphics::title("Continuation Ratio Model CRM Pseudodata")
          }
          if (combine01==TRUE)	{
            graphics::plot(c(0,round(dose90*1.05,0)),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))))))
            crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))))))
            crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,round(dose90*1.05,0),1))))))
            graphics::lines(seq(0,round(dose90*1.05,0),1),1-crequal01,lty=2, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1-crequal01)*(1-crequal2),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,round(dose90*1.05,0),1),(1-crequal01)*(1-crequal2)*(1-crequal3),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            graphics::points(startdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x=(startdose*1.05), y=(targetDLT+.03), labels=startdose)
            graphics::legend(0,1,c("p(y>=2)","p(y>=3)","p(y>=4)","Starting Dose","Dose 10 & Dose 90"),
                   lty=c(2,3,4,-1,-1),pch=c(-1,-1,-1,15,1),lwd=c(2,2,2,1,1),
                   col=c(1,1,1,"red",1),bg="white")
            graphics::title("Continuation Ratio Model CRM Pseudodata")
          }
        }
        if (discrete==TRUE)	{
          if (combine01==FALSE)	{
            lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
            graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0)))),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            tmp0<- reg$coef['Intercept']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp1<- reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp2<- reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            tmp3<- reg$coef['Intercept']+reg$coef['cohort=yall>=3']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1/(1+exp(tmp0))),lty=1, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=2, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1/(1+exp(tmp3)))*(1/(1+exp(tmp2)))*(1/(1+exp(tmp1)))*(1/(1+exp(tmp0))),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            for (i in 1:length(discretedoses))	{
              graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
              graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
            graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
            graphics::legend(0,1,c("p(y>=1)","p(y>=2)","p(y>=3)","p(y>=4)","Starting Discrete Dose","Dose 10 & Dose 90","Estimated Toxicity"),
                   lty=c(1,2,3,4,-1,-1,-1),pch=c(-1,-1,-1,-1,15,1,15),lwd=c(2,2,2,2,1,1,1),
                   col=c(1,1,1,1,"red",1,"gray90"),bg="white")
            graphics::title("Continuation Ratio Model CRM Pseudodata")
          }
          if (combine01==TRUE)	{
            lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
            graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0)))),c(0,1), type="n",ylab="Probability",xlab="Dose")
            graphics::points(c(dose10,dose90),c(0.10, 0.90))
            crequal01<-1/(1+(exp(-(reg$coef['Intercept']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))))))
            crequal2<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=1']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))))))
            crequal3<-1/(1+(exp(-(reg$coef['Intercept']+reg$coef['cohort=yall>=2']+(reg$coef['dose']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))))))
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1-crequal01,lty=2, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1-crequal01)*(1-crequal2),lty=3, col=1, lwd=2)
            graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),(1-crequal01)*(1-crequal2)*(1-crequal3),lty=4, col=1, lwd=2)
            graphics::abline(h=targetDLT,lty=3)
            graphics::abline(v=startdose,lty=3)
            graphics::abline(h=0,lty=1)
            for (i in 1:length(discretedoses))	{
              graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
              graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
            graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
            graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
            graphics::legend(0,1,c("p(y>=2)","p(y>=3)","p(y>=4)","Starting Discrete Dose","Dose 10 & Dose 90","Estimated Toxicity"),
                   lty=c(2,3,4,-1,-1,-1),pch=c(-1,-1,-1,15,1,15),lwd=c(2,2,2,1,1,1),
                   col=c(1,1,1,"red",1,"gray90"),bg="white")
            graphics::title("Continuation Ratio Model CRM Pseudodata")
          }
        }
      }
    }

    ####################################################################
    ###Binary CRM ######################################################
    ####################################################################
    if (design=='CRM')	{

      d<-c(rep(dose10,100),rep(dose90,100))
      y<-c(rep(c(0,1),c(90,10)),rep(c(0,1),c(10,90)))
      reg<-rms::lrm(y~d)

      extraDLT<-c()
      extraDLT[1]<-ifelse(stabilize==TRUE,.30,NA)
      extraDLT[2]<-ifelse(stabilize==TRUE,.50,NA)
      if (is.na(extraDLT[1])==FALSE)	{
        extra<-vector("list",2)
        for (i in 1:2)	{
          dose<-round((log(extraDLT[i]/(1-extraDLT[i]))-reg$coef['Intercept'])/reg$coef['d'])
          prob1<- reg$coef['Intercept']+ (reg$coef['d']*dose)
          prob1<-1/(1+exp(-prob1))
          prob0<-1-prob1
          extra[[i]][1]<-dose
          extra[[i]][2]<-prob0
          extra[[i]][3]<-prob1
        }

        dose30<-round(extra[[1]],2)
        dose50<-round(extra[[2]],2)
        percent30<-round(dose30[2:3]*100,0)
        percent50<-round(dose50[2:3]*100,0)

        pseudoy<-append(y,c(rep(0,percent30[1]),rep(1,percent30[2]),rep(0,percent50[1]),rep(1,percent50[2])),after=length(y))
        pseudod<-append(d,c(rep(round(dose30[1],0),100),rep(round(dose50[1],0),100)),after=length(d))
      }

      if (is.na(extraDLT[1]))	{
        pseudoy<-y
        pseudod<-d
      }

      fullmodel<-rms::lrm(pseudoy~pseudod)
      startdose<-round((log(targetDLT/(1-targetDLT))-fullmodel$coef['Intercept'])/fullmodel$coef['pseudod'])
      if (discrete==TRUE)	{
        discretedoses<-discretedoses[order(discretedoses, decreasing=FALSE)]
        esttoxdiscdose<-diffdose<-absdiffdose<-c()
        for (i in 1:length(discretedoses))	{
          esttoxdiscdose[i]<-1/(1+exp(-1*(fullmodel$coef['Intercept']+ (fullmodel$coef['pseudod']*discretedoses[i]))))
          diffdose[i]<-(startdose-discretedoses[i])
          absdiffdose[i]<-abs(startdose-discretedoses[i])
        }
        tt<-order(absdiffdose, decreasing=FALSE)
        ordernextdose<-cbind(absdiffdose,discretedoses)[tt,]
        if (rounddown==TRUE)	{
          discdose<-ifelse(ordernextdose[1,2]>startdose, discretedoses[which(discretedoses==ordernextdose[1,2])-1] ,ordernextdose[1,2])
        }
        if (rounddown==FALSE)	{
          discdose<-ordernextdose[1,2]
        }
      }
      if (discrete==FALSE)	{
        esttoxdiscdose<-NA
        discdose<-NA
      }

      if (plotit==TRUE)	{
        if (discrete==FALSE)	{
          graphics::plot(c(0,round(dose90*1.05,0)),c(0,1), type="n",ylab="Probability", xlab="Dose")
          graphics::points(c(dose10,dose90),c(0.10, 0.90))
          crm<- fullmodel$coef['Intercept']+(fullmodel$coef['pseudod']*seq(0,round(dose90*1.05,0),1))
          graphics::lines(seq(0,round(dose90*1.05,0),1),1/(1+exp(-crm)),lty=3, col=1, lwd=2)
          graphics::abline(h=targetDLT,lty=3)
          graphics::abline(v=startdose,lty=3)
          graphics::abline(h=0,lty=1)
          graphics::points(startdose,targetDLT,pch=15,cex=1.3,col="red")
          graphics::text(x=(startdose*1.05), y=(targetDLT+.03), labels=startdose)
          graphics::legend(0,1,c("Probability of a DLT","Starting Dose","Dose 10 & Dose 90"),
                 lty=c(3,-1,-1),pch=c(-1,15,1),lwd=c(2,1,1),col=c(1,"red",1),bg="white")
          graphics::title("Binary CRM Pseudodata")
        }
        if (discrete==TRUE)	{
          lengthbar<-(max(discretedoses)-min(discretedoses))/(4*(length(discretedoses)+1))
          graphics::plot(c(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0)))),c(0,1), type="n",ylab="Probability",xlab="Dose")
          graphics::points(c(dose10,dose90),c(0.10, 0.90))
          crm<-fullmodel$coef['Intercept']+(fullmodel$coef['pseudod']*seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1))
          graphics::lines(seq(0,max(round(discretedoses+lengthbar,0),(round(dose90*1.05,0))),1),1/(1+exp(-crm)),lty=3, col=1, lwd=2)
          graphics::abline(h=targetDLT,lty=3)
          graphics::abline(v=startdose,lty=3)
          graphics::abline(h=0,lty=1)
          for (i in 1:length(discretedoses))	{
            graphics::rect(xleft=discretedoses[i]-lengthbar, ybottom=0, xright=discretedoses[i]+lengthbar, ytop=esttoxdiscdose[i],col="gray90")
            graphics::text(x =(discretedoses[i]), y=esttoxdiscdose[i]/2, labels=round(esttoxdiscdose[i]*100,0))}
          graphics::points(discdose,targetDLT,pch=15,cex=1.3,col="red")
          graphics::text(x =(discdose*1.05), y=(targetDLT+.03), labels=discdose)
          graphics::legend(0,1,c("Probability of a DLT","Starting Discrete Dose", "Dose 10 & Dose 90", "Estimated Toxicity"),
                 lty=c(3,-1,-1,-1),pch=c(-1,15,1,15),lwd=c(2,1,1,1),col=c(1,"red",1,"gray90"),bg="white")
          graphics::title("Binary CRM Pseudodata")
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

    return(list('Pseudodata Toxicities'=pseudoy,
                'Pseudodata Doses'=pseudod,
                'Starting Dose'=startdose,
                'Starting Discrete Dose'=discdose,
                'Estimated Toxicity at Discrete Doses'=esttoxdiscdose2,
                'Discrete Doses'=discretedoses,
                'Regression Model'=fullmodel))
    }




