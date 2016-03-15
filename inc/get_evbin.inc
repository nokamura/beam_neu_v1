C -*- fortran -*-
               ymin = yymin(ilike,iD)
               ymax = yymax(ilike,iD)
               include 'bining_x.inc'               
               evnbins(ilike,iD) = nbins
               do i = 0,nbins
                  evx(i,ilike,iD) = x(i)
               enddo
