C -*- fortran -*-
            xmin = dsqrt(ymin)
            xmax = dsqrt(ymax)
            nbins = int( ( xmax -xmin ) / rdx ) ! nbins should be less than 100000
            do i = 0,nbins
               x(i) = xmin +rdx*i
               yy(i) = x(i)**2
            enddo
