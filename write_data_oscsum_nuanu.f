      subroutine write_data_oscsum_nuanu(itype)
C     ****************************************************
C     By Yoshitaro Takaesu @ Okayama U. Dec 9 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS
      integer itype
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,iev,inu
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do ib = -1,1,2  ! 1:nu_mu focusing beam -1:anti-nu_mu focusing beam
         do iD = 1,3  ! 1:SK 2:Oki 3:Korea
            do inu = -1,1,2  ! 1: neutrino is detected  
C                             -1: anti-neutrino is detected
               do iev = 1,6 ! 1:e2e 2:m2e 3:pi0 4:gam 5:e2m 6:m2e
                  call write_data_oscsum_nuanu_unit(itype,ib,iD,inu,iev)
               enddo
            enddo
         enddo
      enddo
            
      return
      end

      subroutine write_data_oscsum_nuanu_unit(itype,ib,iD,inu,iev)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
c      include 'inc/params.inc'
c      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer itype,ib,iD,inu,idet,iev
C     LOCAL VARIABLES 
      character*2 cb,cD,cnu,cdet
      character*5 file_mid
      character*9 file_app
c      integer i
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
CCC   Prepare output files
      if (ib.eq.1) cb = "nb"
      if (ib.eq.-1) cb = "ab"
      if (iD.eq.1) cD = "sk"                  
      if (iD.eq.2) cD = "ok"
      if (iD.eq.3) cD = "kr"
      file_mid = cb//"."//cD

      if (itype.eq.1) then
         file_app = "dat_"//file_mid
      elseif (itype.eq.2) then
         file_app = "the_"//file_mid
      endif
      if (inu.eq.1) then
         if (iev.eq.1) then
            open(1,file=file_app//".nu.e2e.dat",status="replace")
         elseif (iev.eq.2) then
            open(1,file=file_app//".nu.m2m.dat",status="replace")
         elseif (iev.eq.3) then
            open(1,file=file_app//".nu.pi0.dat",status="replace")
         elseif (iev.eq.4) then
            open(1,file=file_app//".nu.gam.dat",status="replace")
         elseif (iev.eq.5) then
            open(1,file=file_app//".nu.m2e.dat",status="replace")
         elseif (iev.eq.6) then
            open(1,file=file_app//".nu.e2m.dat",status="replace")
         endif
      elseif (inu.eq.-1) then
         if (iev.eq.1) then
            open(1,file=file_app//".anu.e2e.dat",status="replace")
         elseif (iev.eq.2) then
            open(1,file=file_app//".anu.m2m.dat",status="replace")
         elseif (iev.eq.3) then
            open(1,file=file_app//".anu.pi0.dat",status="replace")
         elseif (iev.eq.4) then
            open(1,file=file_app//".anu.gam.dat",status="replace")
         elseif (iev.eq.5) then
            open(1,file=file_app//".anu.m2e.dat",status="replace")
         elseif (iev.eq.6) then
            open(1,file=file_app//".anu.e2m.dat",status="replace")
         endif
      endif

CCC   Write results in output files
      call write_data_oscsum_nuanu_atm(1,itype,ib,iD,inu,iev)
      close(1)        

      return
      end


      subroutine write_data_oscsum_nuanu_atm(ilun,itype,ib,iD,inu,iev)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer ilun,itype,ib,iD,inu,idet,iev,ineut
C     LOCAL VARIABLES 
      integer i,j
      real*8 ans(imaxint),tot
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
CCC      write(ilun,*) "# E_nu, total, CCQE-H, CCQE-O,"
CCC     &            ,"CCRes-H, CCRes-O for CC events" 
CCC      write(ilun,*) "# E_nu, total, NCQE, NCRS,"
CCC     &            ," NCDI, NCCO+NCDF for NC events" 

      do i = 0,nbins-1

         tot = 0d0
         do j = 1,imaxint
            ans(j) = 0d0
         enddo

         if (inu.eq.1) then
            if (itype.eq.1) then
               do j = 1,imaxint
                  do idet = 1,3
                     do ineut = -3,3
                        ans(j) = ans(j) +hErec_dat(i+1,j,iev,idet,ineut,
     &                       iD,ib)
                     enddo
                  enddo
                  tot = tot +ans(j)
               enddo
            elseif (itype.eq.2) then
               do j = 1,imaxint
                  do idet = 1,3
                     do ineut = -3,3
                        ans(j) = ans(j) +hErec_th(i+1,j,iev,idet,ineut,
     &                       iD,ib)
                     enddo
                  enddo
                  tot = tot +ans(j)
               enddo
            endif      
         elseif (inu.eq.-1) then   
            if (itype.eq.1) then
               do j = 1,imaxint
                  do idet = -3,-1
                     do ineut = -3,3
                        ans(j) = ans(j) +hErec_dat(i+1,j,iev,idet,ineut,
     &                       iD,ib)
                     enddo
                  enddo
                  tot = tot +ans(j)
               enddo
            elseif (itype.eq.2) then
               do j = 1,imaxint
                  do idet = -3,-1
                     do ineut = -3,3
                        ans(j) = ans(j) +hErec_th(i+1,j,iev,idet,ineut,
     &                       iD,ib)
                     enddo
                  enddo
                  tot = tot +ans(j)
               enddo
            endif      
         endif

         write(ilun,'(f5.3,f10.5,f10.5,f10.5,f10.5,f10.5)') 
     &        (x(i+1)+x(i))/2d0,tot,ans(1),ans(2),ans(3),ans(4)

      enddo
      
      return
      end
