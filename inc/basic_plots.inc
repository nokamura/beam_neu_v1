C -*- fortran -*-
         basic_binsize_tmp = basic_binsize
         binsize_factor_tmp = binsize_factor
	 MH_tmp = int(zz(111))
         L_tmp = zz(116)
         nu_mode_tmp = int(zz(120))
	 detect_nu_tmp = int(zz(118))
         ihfunc_tmp = int(zz(115))
         ihisto_tmp = ihisto

         ymin = 0d0
         ymax = 5d0
         basic_binsize = 0.01
         binsize_factor = 10
         evform = evform_dat
c
c	 Flux	 
c
c         z(116) = 1d0 ! L
c         z(116) = 295d0 ! L
         zz(116) = 0.28d0
         zz(115) = 1 ! ihfunc 
         ihisto = 2
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="flux.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
c
c	 Transition Probability
c
         binsize_factor = 1
         include 'bining_x.inc' ! set basic binning
	 zz(111) = 1 ! NH
         zz(116) = L_tmp ! L
         zz(115) = 2 
         ihisto = 0
	 zz(118) = 1
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_nh_ne.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = 2
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_nh_nm.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -1
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_nh_ae.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -2
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_nh_am.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)

	 zz(111) = -1 ! IH
         zz(116) = 1d3 ! L
         zz(115) = 2 
         ihisto = 0
         binsize_factor = 1
         include 'bining_x.inc' ! set basic binning
	 zz(118) = 1
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_ih_ne.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = 2
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_ih_nm.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -1
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_ih_ae.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -2
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="prob_ih_am.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
c
c	 Xsec
c
         zz(115) = 3
         ihisto = 1
	 zz(150) = 1 ! CC
	 zz(118) = 1 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_cc_ne.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = 2 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_cc_nm.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -1 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_cc_ae.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -2 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_cc_am.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)

	 zz(150) = 0 ! NC
	 zz(118) = 1 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_nc_ne.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = 2 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_nc_nm.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -1 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_nc_ae.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
	 zz(118) = -2 ! detected neutrino
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="xsec_nc_am.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)
   

c
c	 Flux*P*xsec
c
         zz(116) = 1000
         zz(120) = 2
	 zz(118) = 3
	 zz(111) = 1
         basic_binsize = 0.01
         binsize_factor = 5
         zz(115) = 6 ! ihfunc 
         ihisto = 2
         zz(150) = 0
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event1
     &        ,hevent1,nevent1,ierr) 
         open(1,file="flux_P_xsec.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)


c
c	 Events after 1pi0 cut
c
         zz(116) = 1000
         zz(120) = 2
	 zz(118) = 3
	 zz(111) = 1
         basic_binsize = 0.01
         binsize_factor = 5
         zz(115) = 6 ! ihfunc 
         ihisto = 2
         zz(150) = 0
         include 'bining_x.inc' ! set basic binning
         call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins
     &        ,x,evform,serror,snmax,ihisto,event2
     &        ,hevent2,nevent2,ierr) 
         do i = 1,nbins
            E = ( x(i) +x(i-1) )/2d0
            call get_xsecfrac3(E,int(zz(150)),0,5,int(zz(118)),frac)
            event1(i) = frac*event2(i)
         enddo
         open(1,file="pi0event.dat",status="replace")
         do i = 0,nbins-1
            write(1,*) x(i),event1(i+1)
         enddo
         close(1)

         zz(116) = L_tmp
         zz(120) = nu_mode_tmp
	 zz(118) = detect_nu_tmp
	 zz(111) = MH_tmp
         basic_binsize = basic_binsize_tmp
         binsize_factor = binsize_factor_tmp
         zz(115) = ihfunc_tmp
         ihisto = ihisto_tmp
