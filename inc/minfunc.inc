	integer evform_dat,evform_th,iD,iev
	common /minfunc0/ evform_dat,evform_th,iD,iev

        integer evform
        real*8 rnevent_ren
        common /getevent/ rnevent_ren,evform

        integer multi_flag
        common /multi/ multi_flag

	real*8 chisq_stat(eventdim,ndetect,-1:1),chisq_param_tot
	common /chisq/ chisq_stat,chisq_param_tot