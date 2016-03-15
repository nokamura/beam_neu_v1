C -*- fortran -*-
            xmin = Emin
            xmax = Emax
            nbins_basic = int( ( xmax -xmin ) / basic_binsize ) ! nbins should be less than 100000
            do i = 0,nbins_basic
               x_basic(i) = xmin +basic_binsize*i
c               yy(i) = x(i)
            enddo
            if (xmax -x_basic(nbins_basic).gt.1d-7) then
               nbins_basic = nbins_basic +1
               x_basic(nbins_basic) = xmax
c               yy(nbins) = x(nbins)
            endif
