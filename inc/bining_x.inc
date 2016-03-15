C -*- fortran -*-
            xmin = ymin
            xmax = ymax
            binsize = basic_binsize*binsize_factor
            nbins = int( ( xmax -xmin ) / binsize ) ! nbins should be less than 100000
            if (nbins.eq.0) then
               nbins = 1
               x(0) = xmin
               x(1) = xmax
               yy(0) = x(0)
               yy(1) = x(1)
            else
               do i = 0,nbins
                  x(i) = xmin +binsize*i
                  yy(i) = x(i)
               enddo
               if (xmax-x(nbins).gt.1d-7) then
                  nbins = nbins+1
                  x(nbins) = xmax
                  yy(nbins) = x(nbins)
               endif
            endif
