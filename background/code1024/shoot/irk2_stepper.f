      subroutine irk2_stepper(ny, d, xin, yin, xout, yout,
     $      prec_irk, tol, dxin, dxguess, maxits, derivs)

      implicit none
      integer ny, nu, maxits, itsreach
      double precision xin, yin(ny), xout, yout(ny), prec_irk, d, tol
      external derivs
      integer nymax
      include '../nymax.inc'

      double precision errmax, sfty, 
     $   pgrow, pshrink,
     $   xnow, resnorm2, ynorm2
      double precision dxnow, dxin, dxguess
      integer i, j

      parameter(nu=2, sfty = 0.9d0, pgrow = - 0.2d0, 
     $          pshrink = - 0.2d0)
      logical repeat
      double precision res(nymax), ynow(nymax), ynowhalf(nymax)
      
C     Start new calculation
      
      dxnow = dxin
      xnow = xin + dxnow
       
 1    continue
      
      call implicit_step(ny, d, xin, yin, xnow, ynow,
     $     derivs, nu, prec_irk, maxits, itsreach, repeat)

      if (repeat) then
C     If the implicit step did not converge, restart with dx/2            
            dxnow = dxnow / 2.d0
            xnow = xin + dxnow
            write(6,*) 'halved dx at x=',xin 
            goto 1
      else 
            continue
      end if
      
      call implicit_step(ny, d, xin, yin, xin + dxnow / 2.d0,
     $     ynowhalf,derivs, nu, prec_irk, maxits, itsreach, 
     $     repeat)
      call implicit_step(ny, d, xin + dxnow / 2.d0, ynowhalf, 
     $     xnow, ynowhalf, derivs, nu, prec_irk, maxits, itsreach, 
     $     repeat)
 
C     Compute error (Richardson extrapolation)
      errmax = 0.d0
      do j=1,ny
            res(j) = (ynow(j) - ynowhalf(j)) * 16.d0 / 15.d0
      end do

      resnorm2 = 0.d0
      ynorm2 = 0.d0
      do j=1,ny
C           Leaves out the entry in y where Delta is stored            
            if (j.eq.5) cycle
            resnorm2 = resnorm2 + res(j)**2
            ynorm2 = ynorm2 + yin(j)**2
      end do
         
      resnorm2 = sqrt( resnorm2 )
      ynorm2 = sqrt( ynorm2 )
      
      errmax = resnorm2 / (ynorm2 * tol)
C      errmax = resnorm2 / tol
C      write(6,*) errmax, ynorm2

C     Decide whether to grow or shrink stepsize
      if (errmax.gt.1) then
C        Shrink no more than by factor of 10
         dxnow = dxnow * max(0.1d0, sfty * errmax**pshrink)
         xnow = xin + dxnow
         if (xnow.eq.xin) stop 'stepsize underflow in irk2_stepper'
         goto 1
      else
C        Compute guess for next step      
C        Grow not more than by factor of 5
         dxguess = dxnow * 
     $          min(5.d0, max(sfty * errmax**pgrow, sfty))
C            write(6,*) xin, sfty * errmax**pgrow
      end if
C      write(6,*) errmax, dxnow, errmax**pgrow
      xout = xnow
      
      do i=1,ny
         yout(i) = ynow(i)
      end do
      
      return
      
      end
     
      
