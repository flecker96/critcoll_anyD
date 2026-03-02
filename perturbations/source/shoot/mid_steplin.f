C     nu-stage Implicit Runge-Kutta (2nu-order), solved by iteration

      subroutine mid_step(n3, ny, d, Delta, lamb, 
     $                    xin, y3in, xout, y3out,
     $        nu, uBcoll, vBcoll, fBcoll, ia2Bcoll)

      implicit none
      integer ny, n3, nu, lamb
      double precision xin, y3in(n3), xout, y3out(n3), d
      external derivslin

      integer nymax, numax
      include '../nymax.inc'

      parameter(numax=3)
     
      double precision onemin(nymax,nymax), oneplus(nymax,nymax), 
     $    uBcoll(nymax), vBcoll(nymax), ddy3dx(nymax), dy3(nymax),
     $    fBcoll(nymax), ia2Bcoll(nymax), rhs(nymax),
     $    test(nymax,nymax), ddy3dxtest(nymax), rhstest(nymax)
     
      double precision dx, Delta, x
      double precision a(numax,numax), b(numax), c(numax)

      integer i, ii, j, its, indx(nymax)

      double precision y3(n3), par

      


      if (ny .gt. nymax) stop 'nymax too small in implicit_step'
      
C     Start new calculation

      do j=1,n3
         y3(j) = y3in(j)
      end do


C     Collocation points
      dx = ( xout - xin ) 
C      do i=1,nu
C         x(i) = xin + dx * c(i)
C      end do
      x = xin + dx / 2.d0

      do i = 1,n3
         dy3(i) = 0.d0
      end do

      do i = 1, n3
            dy3(i) = 1.d0           
            call derivslin(n3, ny, d, Delta, lamb, dy3, x, ddy3dx,
     $              uBcoll, vBcoll, fBcoll, ia2Bcoll) 
            dy3(i) = 0.d0
            do j = 1, n3
               oneplus(j,i) = 0.5d0 * dx * ddy3dx(j)
               onemin(j,i) = - 0.5d0 * dx * ddy3dx(j)
               test(j,i) = ddy3dx(j)
            end do
            oneplus(i,i) = oneplus(i,i) + 1.d0
            onemin(i,i) = onemin(i,i) + 1.d0
      end do

      do i=1,n3
           rhstest(i) = 0.d0
         do j = 1, n3
            rhstest(i) = rhstest(i) + test(i,j) * y3(j)
         end do
      end do

      call derivslin(n3, ny, d, Delta, lamb, y3, x, ddy3dxtest,
     $              uBcoll, vBcoll, fBcoll, ia2Bcoll) 

      do i=1,n3
         write(6,*) rhstest(i) - ddy3dxtest(i)
      end do

C     Compute rhs for step equation      
      do i = 1,n3
         rhs(i) = 0.d0
         do j = 1, n3
            rhs(i) = rhs(i) + oneplus(i,j) * y3(j)
         end do
      end do
      
C     Calculate suggested change in onemin.y3out = rhs
C     LU decomposition of onemin, stored in the same space
      call ludcmp(onemin, n3, nymax, indx, par)
C     Feed rhs to the subroutine using y3out
      do j=1,n3
          y3out(j) = rhs(j)
      end do
      call lubksb(onemin, n3, nymax, indx, y3out)
C     Improve result
C      call mprove(test, lhsmat, ny, nymax, indx, rhs, yout)
C      call mprove(test, lhsmat, ny, nymax, indx, rhs, yout)
C      call mprove(test, lhsmat, ny, nymax, indx, rhs, yout)
      




      return

      end

