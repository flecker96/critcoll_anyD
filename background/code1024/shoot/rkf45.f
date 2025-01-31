      subroutine rkf45_step(ny, d, xin, yin, xout, yout,
     $     derivs, crit, dxtry, dxnext)

      implicit none
      integer ny, nu, maxits, itsreach
      double precision xin, yin(ny), xout, yout(ny), crit, d
      external derivs
      double precision junkc, junkf1(ny), junkf2(ny), junkf3(ny)

      integer nymax
      include '../nymax.inc'

      parameter(nu=6)
      double precision xi(nymax,nu), xiold(nymax,nu), 
     $   errmax, pointdiff, yscal, sfty, pgrow, pshrink,
     $   errcon
      double precision dx, dxtry, dxnext, dxtemp, 
     $   x(nu), f(nymax,nu)
      double precision A(nu,nu), B(nu), BHAT(nu), C(nu)

      integer i, ii, j, its

      integer ndivs, div
      double precision y4(nymax), y5(nymax)

      if (ny .gt. nymax) stop 'nymax too small in implicit_step'

C     Define RKF45 coefficients
      ! Define the C vector (collocation points)
      C(1) = 0.0D0
      C(2) = 1.0D0 / 4.0D0
      C(3) = 3.0D0 / 8.0D0
      C(4) = 12.0D0 / 13.0D0
      C(5) = 1.0D0
      C(6) = 1.0D0 / 2.0D0

      ! Define the A matrix (Butcher tableau coefficients)
      A(1,1) = 0.0D0
      A(2,1) = 1.0D0 / 4.0D0
      A(3,1) = 3.0D0 / 32.0D0
      A(3,2) = 9.0D0 / 32.0D0
      A(4,1) = 1932.0D0 / 2197.0D0
      A(4,2) = -7200.0D0 / 2197.0D0
      A(4,3) = 7296.0D0 / 2197.0D0
      A(5,1) = 439.0D0 / 216.0D0
      A(5,2) = -8.0D0
      A(5,3) = 3680.0D0 / 513.0D0
      A(5,4) = -845.0D0 / 4104.0D0
      A(6,1) = -8.0D0 / 27.0D0
      A(6,2) = 2.0D0
      A(6,3) = -3544.0D0 / 2565.0D0
      A(6,4) = 1859.0D0 / 4104.0D0
      A(6,5) = -11.0D0 / 40.0D0

      ! Define the B vector (for the 4th order solution)
      B(1) = 25.0D0 / 216.0D0
      B(2) = 0.0D0
      B(3) = 1408.0D0 / 2565.0D0
      B(4) = 2197.0D0 / 4104.0D0
      B(5) = -1.0D0 / 5.0D0
      B(6) = 0.0D0

      ! Define the BHAT vector (for the 5th order solution)
      BHAT(1) = 16.0D0 / 135.0D0
      BHAT(2) = 0.0D0
      BHAT(3) = 6656.0D0 / 12825.0D0
      BHAT(4) = 28561.0D0 / 56430.0D0
      BHAT(5) = -9.0D0 / 50.0D0
      BHAT(6) = 2.0D0 / 55.0D0

C     Parameters for stepper
      sfty = 0.9d0
      pgrow = - 0.2d0
      pshrink = - 0.25d0
      errcon = (5.d0 / sfty)**(1.d0 / pgrow)
      
C     Start new calculation

      do j=1,ny
         y4(j) = yin(j)
         y5(j) = yin(j) 
      end do
       
C     Collocation points
C      dx = ( xout - xin )
      
      dx = dxtry
      
 1    continue
 
C     Collocation points
      do i=1,nu
	 x(i) = xin + dx * C(i)
      end do
      
C     Stages
      do i=1,nu
         do j=1,ny
            xi(j,i) = yin(j)
            do ii=1,i-1
               xi(j,i) = xi(j,i) + dx * A(i,ii) * f(j,ii)
            end do
         end do
         call derivs(ny, d, xi(1,i), x(i), f(1,i), 
     $                   junkc, junkf1,junkf2) 
      end do
      
      do i=1,ny
         do j=1,nu
            y4(i) = y4(i) + dx * B(j) * f(i,j)
            y5(i) = y5(i) + dx * BHAT(j) * f(i,j)
         end do
      end do
     
     
C     Compute maximum pointwise error relative to |y(i)| + |dx * dydx(i)|      
      errmax = 0.d0
      do j=1,ny
           yscal = abs(yin(j)) + abs(dx * xi(j,1))
           pointdiff = abs(y4(j)-y5(j)) / yscal
           errmax = max(errmax, pointdiff)
      end do
      
      errmax = errmax / crit
      
      write(6,*) errmax, dx
      
C     Decide whether to grow or shrink stepsize
      if (errmax.gt.1) then
         dxtemp = sfty * dx * errmax**pshrink
C        Shrink no more than by factor of 10         
         dx = sign(max(abs(dxtemp),0.1d0 * abs(dx)) , dx)
         
         xout = xin + dx  
         if (xout.eq.xin) stop 'stepsize underflow in rkf45'
         goto 1
      else
C        Compute guess for next step      
C        Grow no more than by factor of 5
         if (errmax.gt.errcon) then
            dxnext = sfty * dx * errmax**pgrow
         else
            dxnext = 5.d0 * dx
         end if
         
         xout = xin + dx 
         
         do j=1,ny
            yout(j) = y4(j)    
         end do
      
      end if    
      
      
      return
      
      end
     
      
