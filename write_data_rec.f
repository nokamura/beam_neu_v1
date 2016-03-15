      subroutine write_dat
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 27 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,inu,idet,iev
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ib = 1
      iD = 1
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_nb.sk.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.sk.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.sk.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_nb.sk.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.sk.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.sk.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_nb.sk.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.sk.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.sk.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_nb.sk.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.sk.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.sk.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_nb.sk.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.sk.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.sk.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_nb.sk.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.sk.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.sk.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.sk.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.sk.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.sk.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.sk.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.sk.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.sk.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 2
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_nb.oki.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.oki.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.oki.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_nb.oki.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.oki.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.oki.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_nb.oki.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.oki.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.oki.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_nb.oki.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.oki.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.oki.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_nb.oki.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.oki.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.oki.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_nb.oki.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.oki.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.oki.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.oki.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.oki.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.oki.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.oki.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.oki.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.oki.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 3
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_nb.kr.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.kr.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.kr.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_nb.kr.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.kr.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.kr.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_nb.kr.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_nb.kr.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_nb.kr.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_nb.kr.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.kr.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.kr.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_nb.kr.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.kr.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.kr.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_nb.kr.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_nb.kr.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_nb.kr.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_nb.kr.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_nb.kr.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_nb.kr.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_nb.kr.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_nb.kr.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_nb.kr.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      ib = -1
      iD = 1
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_ab.sk.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.sk.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.sk.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_ab.sk.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.sk.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.sk.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_ab.sk.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.sk.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.sk.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_ab.sk.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.sk.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.sk.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_ab.sk.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.sk.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.sk.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_ab.sk.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.sk.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.sk.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.sk.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.sk.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.sk.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.sk.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.sk.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.sk.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 2
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_ab.oki.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.oki.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.oki.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_ab.oki.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.oki.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.oki.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_ab.oki.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.oki.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.oki.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_ab.oki.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.oki.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.oki.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_ab.oki.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.oki.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.oki.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_ab.oki.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.oki.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.oki.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.oki.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.oki.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.oki.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.oki.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.oki.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.oki.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 3
      inu = 1
      idet = 1
      iev = 1
      open(1,file="dat_ab.kr.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.kr.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.kr.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="dat_ab.kr.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.kr.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.kr.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="dat_ab.kr.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="dat_ab.kr.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="dat_ab.kr.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="dat_ab.kr.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.kr.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.kr.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="dat_ab.kr.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.kr.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.kr.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="dat_ab.kr.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="dat_ab.kr.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="dat_ab.kr.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="dat_ab.kr.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="dat_ab.kr.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="dat_ab.kr.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="dat_ab.kr.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="dat_ab.kr.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="dat_ab.kr.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_dat(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      return
      end


      subroutine write_th
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 27 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,inu,idet,iev
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ib = 1
      iD = 1
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_nb.sk.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.sk.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.sk.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_nb.sk.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.sk.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.sk.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_nb.sk.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.sk.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.sk.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_nb.sk.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.sk.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.sk.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_nb.sk.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.sk.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.sk.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_nb.sk.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.sk.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.sk.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.sk.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.sk.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.sk.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.sk.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.sk.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.sk.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 2
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_nb.oki.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.oki.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.oki.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_nb.oki.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.oki.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.oki.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_nb.oki.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.oki.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.oki.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_nb.oki.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.oki.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.oki.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_nb.oki.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.oki.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.oki.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_nb.oki.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.oki.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.oki.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.oki.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.oki.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.oki.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.oki.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.oki.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.oki.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 3
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_nb.kr.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.kr.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.kr.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_nb.kr.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.kr.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.kr.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_nb.kr.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_nb.kr.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_nb.kr.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_nb.kr.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.kr.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.kr.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_nb.kr.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.kr.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.kr.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_nb.kr.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_nb.kr.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_nb.kr.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_nb.kr.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_nb.kr.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_nb.kr.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_nb.kr.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_nb.kr.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_nb.kr.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      ib = -1
      iD = 1
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_ab.sk.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.sk.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.sk.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_ab.sk.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.sk.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.sk.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_ab.sk.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.sk.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.sk.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_ab.sk.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.sk.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.sk.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_ab.sk.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.sk.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.sk.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_ab.sk.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.sk.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.sk.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.sk.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.sk.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.sk.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.sk.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.sk.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.sk.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 2
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_ab.oki.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.oki.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.oki.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_ab.oki.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.oki.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.oki.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_ab.oki.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.oki.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.oki.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_ab.oki.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.oki.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.oki.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_ab.oki.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.oki.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.oki.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_ab.oki.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.oki.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.oki.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.oki.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.oki.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.oki.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.oki.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.oki.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.oki.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      iD = 3
      inu = 1
      idet = 1
      iev = 1
      open(1,file="th_ab.kr.ne.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ne.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ne.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ne.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ne.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ne.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.kr.ne.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ne.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ne.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ne.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ne.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ne.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.kr.ne.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ne.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ne.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ne.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ne.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ne.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ne.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = 2
      idet = 1
      iev = 1
      open(1,file="th_ab.kr.nm.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nm.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nm.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nm.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nm.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nm.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.kr.nm.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nm.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nm.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nm.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nm.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nm.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.kr.nm.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nm.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nm.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nm.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nm.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nm.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nm.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      inu = 3
      idet = 1
      iev = 1
      open(1,file="th_ab.kr.nt.ne.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.ne.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nt.ne.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.ne.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nt.ne.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nt.ne.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nt.ne.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nt.ne.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 2
      iev = 1
      open(1,file="th_ab.kr.nt.nm.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.nm.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nt.nm.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.nm.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nt.nm.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nt.nm.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nt.nm.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nt.nm.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = 3
      iev = 1
      open(1,file="th_ab.kr.nt.nt.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.nt.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.nt.nt.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.nt.nt.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.nt.nt.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.nt.nt.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.nt.nt.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.nt.nt.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -1
      idet = -1
      iev = 1
      open(1,file="th_ab.kr.ae.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ae.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ae.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ae.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ae.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ae.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.kr.ae.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ae.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ae.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ae.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ae.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ae.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.kr.ae.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.ae.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.ae.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.ae.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.ae.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.ae.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.ae.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -2
      idet = -1
      iev = 1
      open(1,file="th_ab.kr.am.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.am.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.am.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.am.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.am.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.am.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.kr.am.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.am.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.am.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.am.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.am.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.am.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.kr.am.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.am.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.am.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.am.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.am.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.am.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.am.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)


      inu = -3
      idet = -1
      iev = 1
      open(1,file="th_ab.kr.at.ae.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.ae.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.at.ae.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.ae.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.at.ae.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.at.ae.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.at.ae.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.at.ae.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -2
      iev = 1
      open(1,file="th_ab.kr.at.am.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.am.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.at.am.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.am.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.at.am.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.at.am.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.at.am.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.at.am.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      idet = -3
      iev = 1
      open(1,file="th_ab.kr.at.at.e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.at.e2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 2
      open(1,file="th_ab.kr.at.at.m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_raw_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      open(1,file="th_ab.kr.at.at.m2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 3
      open(1,file="th_ab.kr.at.at.pi0.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 4
      open(1,file="th_ab.kr.at.at.gam.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 5
      open(1,file="th_ab.kr.at.at.m2e.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)
      iev = 6
      open(1,file="th_ab.kr.at.at.e2m.dat",status="replace")
      do i = 0,nbins-1
         write(1,*) x(i),hErec_th(i+1,iev,idet,inu,iD,ib)
      enddo
      close(1)

      return
      end


