      subroutine write_chi2event_dat
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Feb 13 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,ilike
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ib = 1
      iD = 1
      ilike = 1
      open(1,file="chi2dat_nb.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_nb.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2dat_nb.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_nb.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2dat_nb.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_nb.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      ib = -1
      iD = 1
      ilike = 1
      open(1,file="chi2dat_ab.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_ab.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2dat_ab.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_ab.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2dat_ab.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2dat_ab.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_dat(i+1,ilike,iD,ib)
      enddo
      close(1)

      return
      end


      subroutine write_chi2event_th
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Feb 13 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,ilike
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ib = 1
      iD = 1
      ilike = 1
      open(1,file="chi2th_nb.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_nb.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2th_nb.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_nb.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2th_nb.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_nb.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      ib = -1
      iD = 1
      ilike = 1
      open(1,file="chi2th_ab.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_ab.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2th_ab.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_ab.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2th_ab.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2th_ab.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),hErec_chi2_th(i+1,ilike,iD,ib)
      enddo
      close(1)

      return
      end


      subroutine write_chi2_stat
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Feb 13 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i,j,k
      integer ib,iD,ilike
      real*8 chi2_stat(nmode_like,ndetect,-1:1),chi2_stat_tot
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do i = -1,1,2
         do j = 1,ndetect
            do k = 1,nmode_like
               chi2_stat(k,j,i) = 0d0
            enddo
         enddo
      enddo

      ib = 1
      iD = 1
      ilike = 1
      open(1,file="chi2stat_nb.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_nb.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2stat_nb.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_nb.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2stat_nb.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_nb.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      ib = -1
      iD = 1
      ilike = 1
      open(1,file="chi2stat_ab.sk.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_ab.sk.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 2
      ilike = 1
      open(1,file="chi2stat_ab.ok.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_ab.ok.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      iD = 3
      ilike = 1
      open(1,file="chi2stat_ab.kr.elike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)
      ilike = 2
      open(1,file="chi2stat_ab.kr.mlike.dat",status="replace")
      do i = 0,nbins_chi2(ilike,iD)-1
         write(1,*) x_chi2(i,ilike,iD),chi2_stat_bin(i+1,ilike,iD,ib)
         chi2_stat(ilike,iD,ib) = chi2_stat(ilike,iD,ib)
     &        +chi2_stat_bin(i+1,ilike,iD,ib)
      enddo
      close(1)

      chi2_stat_tot = 0d0
      do i = -1,1,2
         do j = 1,ndetect
            do k = 1,nmode_like
               chi2_stat_tot = chi2_stat_tot +chi2_stat(k,j,i)
            enddo
         enddo
      enddo

      open(1,file="chi2stat_all.dat",status="replace")
      do i = -1,1,2
         do j = 1,ndetect
            write(1,*) chi2_stat(1,j,i), chi2_stat(2,j,i)
         enddo
      enddo
      write(1,*) chi2_stat_tot
      close(1)

      return
      end
