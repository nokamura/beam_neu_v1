      program smeardist_mono
      implicit none
      include 'inc/params.inc'
      integer i
      integer inutype,int_mode,evform,icc,ierr
      real*8 Ev,rievent,event_mode(1000),xmin,xmax,rnevent_mode
      real*8 rtotevent,frac,rievent2,event_total(1000),rnevent_total
      real*8 z(3),hevent_mode(1000)
      real*8 fErec_dist
      external fErec_dist

      evform = 2
      serror = 1d-2
      snmax = 10
      ihisto = 2

      Ev = 2d0 ! incoming neutrino energy in [GeV]
      ymin = 0d0
      ymax = 3d0
      rievent = 1d0

      basic_binsize = 0.05d0
      binsize_factor = 1
      include 'inc/bining_x.inc'
      icc = 1

      inutype = 1 ! detected neutrino flavor: 1:nu_e 2:nu_mu -1:anti-nu_e -2: anti-nu_mu
      int_mode = 1 ! 1:CCQE-H 2:CCQE-O 3:RES-H 4:RES-O
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ne_ccqeh_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 2
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ne_ccqeo_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 3
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ne_ccrsh_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 4
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ne_ccrso_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 1,2
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ne_ccqe_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 3,4
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ne_res_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 1,4
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ne_total_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)


      inutype = -1 ! detected neutrino flavor: 1:nu_e 2:nu_mu -1:anti-nu_e -2: anti-nu_mu
      int_mode = 1 ! 1:CCQE-H 2:CCQE-O 3:RES-H 4:RES-O
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ae_ccqeh_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 2
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ae_ccqeo_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 3
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ae_ccrsh_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      int_mode = 4
      include 'inc/get_nudist_mode.inc'
      open(1,file="smeardist_ae_ccrso_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_mode(i+1),rnevent_mode
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 1,2
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ae_ccqe_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 3,4
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ae_res_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)

      do i = 1,nbins
         event_total(i) = 0d0
      enddo
      rnevent_total = 0d0
      do int_mode = 1,4
         include 'inc/get_nudist_mode.inc'
         do i = 1,nbins
            event_total(i) = event_total(i) +event_mode(i)
         enddo
         rnevent_total = rnevent_total +rnevent_mode
      enddo
      open(1,file="smeardist_ae_total1_2gev.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),event_total(i+1),rnevent_total
      enddo
      close(1)

      end
