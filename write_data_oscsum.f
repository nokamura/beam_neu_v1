      subroutine write_data_oscsum(itype)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
c      include 'inc/params.inc'
c      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer itype
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,iev 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do ib = -1,1,2
         do iD = 1,3
            do iev = 1,6
               call write_data_oscsum_unit(itype,ib,iD,iev)
            enddo
         enddo
      enddo
            
      return
      end

      subroutine write_data_oscsum_unit(itype,ib,iD,iev)
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
      if (iev.eq.1) then
         open(1,file=file_app//".e2e.dat",status="replace")
      elseif (iev.eq.2) then
         open(1,file=file_app//".m2m.dat",status="replace")
      elseif (iev.eq.3) then
         open(1,file=file_app//".pi0.dat",status="replace")
      elseif (iev.eq.4) then
         open(1,file=file_app//".gam.dat",status="replace")
      elseif (iev.eq.5) then
         open(1,file=file_app//".m2e.dat",status="replace")
      elseif (iev.eq.6) then
         open(1,file=file_app//".e2m.dat",status="replace")
      endif
      call write_data_oscsum_atm(1,itype,ib,iD,iev)
      close(1)        

      return
      end


      subroutine write_data_oscsum_atm(ilun,itype,ib,iD,iev)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer ilun,itype,ib,iD,inu,idet,iev
C     LOCAL VARIABLES 
      integer i,j
      real*8 ans(imaxint),tot
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
c      write(ilun,*) "# E_nu, total, CCQE-H, CCQE-O,"
c     &            ,"CCRes-H, CCRes-O for CC events" 
c      write(ilun,*) "# E_nu, total, NCQE, NCRS,"
c     &            ," NCDI, NCCO+NCDF for NC events" 
      do i = 0,nbins-1
         tot = 0d0
         do j = 1,imaxint
            ans(j) = 0d0
         enddo
         if (itype.eq.1) then
            do j = 1,imaxint
               do idet = -3,3
                  do inu = -3,3
                    ans(j) = ans(j) +hErec_dat(i+1,j,iev,idet,inu,iD,ib)
                  enddo
               enddo
               tot = tot +ans(j)
            enddo
         elseif (itype.eq.2) then
            do j = 1,imaxint
               do idet = -3,3
                  do inu = -3,3
                     ans(j) = ans(j) +hErec_th(i+1,j,iev,idet,inu,iD,ib)
                  enddo
               enddo
               tot = tot +ans(j)
            enddo
         endif         
         write(ilun,'(f5.3,f10.5,f10.5,f10.5,f10.5,f10.5)') 
     &        x(i),tot,ans(1),ans(2),ans(3),ans(4)
      enddo
      
      return
      end
