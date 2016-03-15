* HEigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th
************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.
      subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
      implicit none
      integer n, ldA, ldU, sort, MAXDIM
      parameter (MAXDIM=16)
      complex*16 A(ldA,n), U(ldU,n)
      real*8 d(n),SYM_EPS
      parameter (SYM_EPS=2D0**(-103))
      
      integer p, q, j
      real*8 red, off, thresh
      real*8 t, delta, invc, s
      complex*16 x, y, Apq
      real*8 ev(2,MAXDIM)
      
      integer sweep
      common /nsweeps/ sweep
      
      real*8 Sq
      complex*16 c
      Sq(c) = dble(c*dconjg(c))
      
      if( n .gt. MAXDIM ) then
         print *, "Dimension too large"
         d(1) = -999
         return
      endif
      
      do p = 1, n
         ev(1,p) = 0
         ev(2,p) = dble(A(p,p))
         d(p) = ev(2,p)
      enddo
      
      do p = 1, n
         do q = 1, n
            U(q,p) = 0
         enddo
         U(p,p) = 1
      enddo
      
      red = .04D0/n**4
      
      do sweep = 1, 50
         off = 0
         do q = 2, n
            do p = 1, q - 1
               off = off + Sq(A(p,q))
            enddo
         enddo
         if( .not. off .gt. SYM_EPS ) goto 100
         
         thresh = 0
         if( sweep .lt. 4 ) thresh = off*red
         
         do q = 2, n
            do p = 1, q - 1
               Apq = A(p,q)
               off = Sq(Apq)
               if( sweep .gt. 4 .and. off .lt.
     &              SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
                  A(p,q) = 0
               else if( off .gt. thresh ) then
                  t = .5D0*(ev(2,p) - ev(2,q))
                  t = 1/(t + sign(sqrt(t**2 + off), t))
                  
                  delta = t*off
                  ev(1,p) = ev(1,p) + delta
                  ev(2,p) = d(p) + ev(1,p)
                  ev(1,q) = ev(1,q) - delta
                  ev(2,q) = d(q) + ev(1,q)
                  
                  invc = sqrt(delta*t + 1)
                  s = t/invc
                  t = delta/(invc + 1)
                  
                  do j = 1, p - 1
                     x = A(j,p)
                     y = A(j,q)
                     A(j,p) = x + s*(dconjg(Apq)*y - t*x)
                     A(j,q) = y - s*(Apq*x + t*y)
                  enddo
                  
                  do j = p + 1, q - 1
                     x = A(p,j)
                     y = A(j,q)
                     A(p,j) = x + s*(Apq*dconjg(y) - t*x)
                     A(j,q) = y - s*(Apq*dconjg(x) + t*y)
                  enddo
                  
                  do j = q + 1, n
                     x = A(p,j)
                     y = A(q,j)
                     A(p,j) = x + s*(Apq*y - t*x)
                     A(q,j) = y - s*(dconjg(Apq)*x + t*y)
                  enddo
                  
                  A(p,q) = 0
                  
                  do j = 1, n
                     x = U(p,j)
                     y = U(q,j)
                     U(p,j) = x + s*(Apq*y - t*x)
                     U(q,j) = y - s*(dconjg(Apq)*x + t*y)
                  enddo
               endif
            enddo
         enddo
         
         do p = 1, n
            ev(1,p) = 0
            d(p) = ev(2,p)
         enddo
      enddo
      
      print *, "Bad convergence in HEigensystem"
      
 100  if( sort .eq. 0 ) return
      
*     sort the eigenvalues
      
      do p = 1, n - 1
         j = p
         t = d(p)
         do q = p + 1, n
            if( sort*(t - d(q)) .gt. 0 ) then
               j = q
               t = d(q)
            endif
         enddo
         
         if( j .ne. p ) then
            d(j) = d(p)
            d(p) = t
            do q = 1, n
               x = U(p,q)
               U(p,q) = U(j,q)
               U(j,q) = x
            enddo
         endif
      enddo
      end
      


* CEigensystem.F
* diagonalization of a complex n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th
************************************************************************
** CEigensystem diagonalizes a general complex n-by-n matrix.
** Input: n, A = n-by-n matrix
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^-1.
      subroutine CEigensystem(n, A,ldA, d, U,ldU, sort)
      implicit none
      integer n, ldA, ldU, sort
      complex*16 A(ldA,n), U(ldU,n), d(n)

      integer p, q, j, MAXDIM
      real*8 red, off, thresh,EPS
      parameter (EPS=2D0**(-102))
      complex*16 delta, t, s, invc, sx, sy, tx, ty
      complex*16 x, y
      parameter (MAXDIM=16)
      complex*16 ev(2,MAXDIM)

      integer sweep
      common /nsweeps/ sweep

      real*8 Sq
      complex*16 c
      Sq(c) = dble(c*dconjg(c))
      
      if( n .gt. MAXDIM ) then
         print *, "Dimension too large"
         d(1) = -999
         return
      endif
      
      do p = 1, n
         ev(1,p) = 0
         ev(2,p) = A(p,p)
         d(p) = ev(2,p)
      enddo
      
      do p = 1, n
         do q = 1, n
            U(q,p) = 0
         enddo
         U(p,p) = 1
      enddo
      
      red = .01D0/n**4
      
      do sweep = 1, 100
         off = 0
         do q = 2, n
            do p = 1, q - 1
               off = off + Sq(A(p,q)) + Sq(A(q,p))
            enddo
         enddo
         if( .not. off .gt. EPS ) goto 100
         
         thresh = 0
         if( sweep .lt. 4 ) thresh = off*red
         
         do q = 2, n
            do p = 1, q - 1
               off = Sq(A(p,q)) + Sq(A(q,p))
               if( sweep .gt. 4 .and. off .lt.
     &              EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
                  A(p,q) = 0
                  A(q,p) = 0
               else if( off .gt. thresh ) then
                  delta = A(p,q)*A(q,p)
                  x = .5D0*(ev(2,p) - ev(2,q))
                  y = sqrt(x**2 + delta)
                  t = x - y
                  s = x + y
                  if( Sq(t) .lt. Sq(s) ) t = s
                  
                  t = 1/t
                  delta = delta*t
                  ev(1,p) = ev(1,p) + delta
                  ev(2,p) = d(p) + ev(1,p)
                  ev(1,q) = ev(1,q) - delta
                  ev(2,q) = d(q) + ev(1,q)
                  
                  invc = sqrt(delta*t + 1)
                  s = t/invc
                  t = t/(invc + 1)
                  sx = s*A(p,q)
                  ty = t*A(p,q)
                  sy = s*A(q,p)
                  tx = t*A(q,p)
                  
                  do j = 1, n
                     x = A(j,p)
                     y = A(j,q)
                     A(j,p) = x + sy*(y - ty*x)
                     A(j,q) = y - sx*(x + tx*y)
                     x = A(p,j)
                     y = A(q,j)
                     A(p,j) = x + sx*(y - tx*x)
                     A(q,j) = y - sy*(x + ty*y)
                  enddo
                  
                  A(p,q) = 0
                  A(q,p) = 0
                  
                  do j = 1, n
                     x = U(p,j)
                     y = U(q,j)
                     U(p,j) = x + sx*(y - tx*x)
                     U(q,j) = y - sy*(x + ty*y)
                  enddo
               endif
            enddo
         enddo
         
         do p = 1, n
            ev(1,p) = 0
            d(p) = ev(2,p)
         enddo
      enddo
      
      print *, "Bad convergence in CEigensystem"
      
 100  if( sort .eq. 0 ) return
      
*     sort the eigenvalues by their real part
      
      do p = 1, n - 1
         j = p
         t = d(p)
         do q = p + 1, n
            if( sort*(dble(t) - dble(d(q))) .gt. 0 ) then
               j = q
               t = d(q)
            endif
         enddo
         
         if( j .ne. p ) then
            d(j) = d(p)
            d(p) = t
            do q = 1, n
               x = U(p,q)
               U(p,q) = U(j,q)
               U(j,q) = x
            enddo
         endif
      enddo
      end
