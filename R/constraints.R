constraints <-
  function(dose=dose, lastdose=lastdose, actualdltcohort=actualdltcohort, numberdltrule=NA,
           lowerlimitrule=NA, upperlimitrule=NA, dltrule=NA,
           increaserule=NA, minimum=NA, maximum=NA)
  {

    #WARNINGS#
    if ((is.na(numberdltrule)==FALSE) & (is.na(dltrule)==TRUE))	{
      stop('You must specify how many DLTs in a cohort are required to use the DLT rule
           (i.e. 2 patients experience a DLT in a cohort)')
    }
    if ((is.na(lowerlimitrule)==FALSE) & (lowerlimitrule>0 & lowerlimitrule<1) & (is.na(minimum)==TRUE|is.na(maximum)))	{
      stop('If you want the lower limit rule to be scaled percentage of the range of doses,
           you must specify a minimum and maximum dose')
    }
    if ((is.na(lowerlimitrule)==FALSE) & lowerlimitrule<0)	{
      stop('If you wish to use the lowerlimitrule you must specify either a number between 0 and 1 (indicating
           a minimum percentage of the range of doses) or a number equal to 0 or greater than 1 (indicating a minimum
           dose level to test)')
    }
    if ((is.na(upperlimitrule)==FALSE) & (upperlimitrule<1 & upperlimitrule>0) & (is.na(minimum)==TRUE|is.na(maximum)))	{
      stop('If you want the upper limit rule to be scaled percentage of the range of doses,
           you must specify a minimum and maximum dose')
    }
    if ((is.na(upperlimitrule)==FALSE) & upperlimitrule<=0)	{
      stop('If you wish to use the upperlimitrule you must specify either a number between 0 and 1 (indicating
           a maximum percentage of the range of doses) or a number greater than 1 (indicating a maximum
           dose level to test)')
    }
    if ((is.na(maximum)==FALSE) & maximum<=0)	{
      stop('If you wish to use the upperlimitrule in terms of a percentage you must specify a
           maximum dose level to test greater than 0')
    }
    if ((is.na(minimum)==FALSE) & minimum<=0)	{
      stop('If you wish to use the lowerlimitrule in terms of a percentage you must specify a
           minimum dose level to test greater than 0')
    }
    if ((is.na(minimum)==FALSE) & is.na(maximum)==FALSE & (minimum>maximum | minimum==maximum))	{
      stop('If you wish to use the lowerlimitrule and/or upperlimitrule in terms of a percentage
           you must specify a minimum dose level less than and not equal to the maximum dose you wish
           to consider')
    }
    if ((is.na(dltrule)==FALSE) & dltrule<=0)	{
      stop('If you wish to use the dltrule you must specify either a number between 0 and 1 (indicating
           the minimum percentage of the last dose tested to decrease by) or a number greater than 1 (indicating a
           minimum dose level to decrease with the next tested cohort)')
    }
    if ((is.na(increaserule)==FALSE) & increaserule<=0)	{
      stop('If you wish to use the increaserule you must specify either a number between 0 and 1 (indicating
           the maximum percentage of the last dose tested to increase by) or a number greater than 1 (indicating a
           maximum dose level to increase with the next tested cohort)')
    }

    dosehat.original<-dose
    dosehat.new<-dose
    dosehat<-lastdose

    ################
    #DLT CONSTRAINT#
    ################
    if (is.na(dltrule)==FALSE)	{

      if (dltrule>0 & dltrule<1 & is.na(lowerlimitrule)==FALSE)	{
        ###DLTrule is in terms of percentage only###
        if (lowerlimitrule>=1 | lowerlimitrule==0)	{
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-(dosehat*dltrule))) &
                                    ((dosehat-(dosehat*dltrule))<lowerlimitrule)),
                                 lowerlimitrule, dosehat.new)
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dltdosehat.new>(dosehat-(dosehat*dltrule)))&
                                    (dltdosehat.new!=lowerlimitrule)), dosehat-(dosehat*dltrule), dltdosehat.new)
        }
        if (lowerlimitrule<1 & lowerlimitrule>0)	{
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-(dosehat*dltrule))) &
                                    ((dosehat-(dosehat*dltrule))<(minimum+((maximum-minimum)*lowerlimitrule)))),
                                 (minimum+((maximum-minimum)*lowerlimitrule)), dosehat.new)
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dltdosehat.new>(dosehat-(dosehat*dltrule)))&
                                    (dltdosehat.new!=(minimum+((maximum-minimum)*lowerlimitrule)))),
                                 dosehat-(dosehat*dltrule), dltdosehat.new)
        }
      }

      if ((dltrule>=1 | dltrule==0) & is.na(lowerlimitrule)==FALSE)	{
        ###DLTrule is a dose level to decrease by###
        if (lowerlimitrule>=1 | lowerlimitrule==0)	{
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-dltrule)) &
                                    ((dosehat-dltrule)<lowerlimitrule)),
                                 lowerlimitrule, dosehat.new)
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dltdosehat.new>(dosehat-dltrule)) &
                                    (dltdosehat.new!=lowerlimitrule)),(dosehat-dltrule), dltdosehat.new)
        }
        if (lowerlimitrule<1 & lowerlimitrule>0)	{
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-dltrule)) &
                                    ((dosehat-dltrule)<(minimum+((maximum-minimum)*lowerlimitrule)))),
                                 (minimum+((maximum-minimum)*lowerlimitrule)), dosehat.new)
          dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dltdosehat.new>(dosehat-dltrule))&
                                    (dltdosehat.new!=((minimum+((maximum-minimum)*lowerlimitrule))))),
                                 (dosehat-dltrule), dltdosehat.new)
        }
      }

      ##No lowerlimitrule specifed##
      if (dltrule>0 & dltrule<1 & is.na(lowerlimitrule)==TRUE)	{
        ###DLTrule is in terms of percentage only###
        dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-(dosehat*dltrule)))),
                               dosehat-(dosehat*dltrule), dosehat.new)
      }
      if ((dltrule>=1 | dltrule==0) & is.na(lowerlimitrule)==TRUE)	{
        ###DLTrule is ONLY a raw number to decrease###
        dltdosehat.new<-ifelse(((actualdltcohort>=numberdltrule) & (dosehat.new>(dosehat-dltrule))),
                               (dosehat-dltrule), dosehat.new)
      }
      dltruleused<-ifelse(dltdosehat.new!=dosehat.original, 1, 0)
      dosehat.new<-ifelse(dltruleused==1,dltdosehat.new, dosehat.new)
    }
    dltruleused<-ifelse(is.na(dltrule),0, dltruleused)

    ###############
    #INCREASE RULE#
    ###############
    if (is.na(increaserule)==FALSE)	{

      ###Increase Rule is in terms of percentage###
      if (increaserule>0 & increaserule<1)	{
        incrdosehat.new<-ifelse(dosehat.new>dosehat+(dosehat*increaserule),
                                dosehat+(dosehat*increaserule),dosehat.new)
      }

      ###Increase Rule is a dose level to increase###
      if (increaserule==0 | increaserule>=1)	{
        incrdosehat.new<-ifelse(dosehat.new-increaserule>dosehat,
                                dosehat+increaserule,dosehat.new)
      }

      incrruleused<-ifelse(incrdosehat.new!=dosehat.new,1,0)
      dosehat.new<-ifelse(incrdosehat.new!=dosehat.new,incrdosehat.new,dosehat.new)
    }
    incrruleused<-ifelse(is.na(increaserule),0,incrruleused)

    ##################
    #LOWER LIMIT RULE#
    ##################
    if (is.na(lowerlimitrule)==FALSE)	{
      ###>=1 implies that the lowerlimit rule is a fixed number###
      if (lowerlimitrule >= 1 | lowerlimitrule==0)	{
        lowdosehat.new<-ifelse(dosehat.new<lowerlimitrule,lowerlimitrule,dosehat.new)
      }
      ###<1 implies that the lowerlimit rule is a scaled percentage based on max and min values###
      if (lowerlimitrule < 1 & lowerlimitrule>0)	{
        lowdosehat.new<-ifelse(dosehat.new<(minimum+((maximum-minimum)*lowerlimitrule)),
                               (minimum+((maximum-minimum)*lowerlimitrule)),dosehat.new)
      }
      lowstopruleused<-ifelse(lowdosehat.new!=dosehat.new,1,0)
      dosehat.new<-ifelse(lowdosehat.new!=dosehat.new,lowdosehat.new,dosehat.new)
    }
    lowstopruleused<-ifelse(is.na(lowerlimitrule),0,lowstopruleused)

    ##################
    #UPPER LIMIT RULE#
    ##################
    if (is.na(upperlimitrule)==FALSE)	{
      ###>=1 implies that the upperlimit rule is a fixed number###
      if (upperlimitrule >= 1 | upperlimitrule==0)	{
        updosehat.new<-ifelse(dosehat.new>upperlimitrule,upperlimitrule,dosehat.new)
      }
      ###<1 implies that the upperlimit rule is a scaled percentage based on max and min values###
      if (upperlimitrule < 1 & upperlimitrule>0)		{
        updosehat.new<-ifelse(dosehat.new>(maximum-((maximum-minimum)*upperlimitrule)),
                              (maximum-((maximum-minimum)*upperlimitrule)),dosehat.new)
      }
      upstopruleused<-ifelse(updosehat.new!=dosehat.new,1,0)
      dosehat.new<-ifelse(updosehat.new!=dosehat.new,updosehat.new,dosehat.new)
    }
    upstopruleused<-ifelse(is.na(upperlimitrule),0,upstopruleused)

    rulesused<-c(dltruleused, incrruleused, lowstopruleused, upstopruleused)
    constraintused<-"None Used"
    if (rulesused[1]==1) 	{
      constraintused<-ifelse(dosehat.new==lowerlimitrule | dosehat.new==(minimum+((maximum-minimum)*lowerlimitrule)),
                             "Lower Bound Stop Rule Used","DLT Rule Used")
    }
    constraintused<-ifelse(rulesused[2]==1,"Increase Rule Used", constraintused)
    constraintused<-ifelse(rulesused[3]==1,"Lower Bound Stop Rule Used", constraintused)
    constraintused<-ifelse(rulesused[4]==1,"Upper Bound Stop Rule Used", constraintused)
    dosehat.new<-ifelse(sum(rulesused)==0,dosehat.original,dosehat.new)

    return(list(dosehat.new=dosehat.new,rulesused=rulesused,constraintused=constraintused))
    }

