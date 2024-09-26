      subroutine derivslin(ny, d, Delta, lamb, y, x, dydx, 
     $                      uB, vB, fB, ia2B)

      implicit none
      integer ny
      double precision y(ny), x, d, lamb,
     $        dydx(ny)

      integer nymax
      include '../nymax.inc'

      integer i
      double precision Delta, ia2(nymax),
     $     u(nymax), v(nymax), f(nymax), 
     $     dudtau(nymax), dvdtau(nymax), 
     $     dfdtau(nymax), 
     $     dudx(nymax), dvdx(nymax), djunkxidx(nymax), dfdx(nymax),
     $     uB(nymax), vB(nymax), fB(nymax), ia2B(nymax),
     $     dudtauB(nymax), dvdtauB(nymax)
      
      
C     Construct required background derivatives
      call fourdiff1(ny, uB, dudtauB, Delta)
      call fourdiff1(ny, vB, dvdtauB, Delta)
      
      
C     Reconstruct fields in tau-space.
      call fieldsfromylin(ny, d, y, x,
     $     u, v, f, ia2, Delta, lamb,
     $     dudtau, dvdtau, dfdtau,
     $     uB, vB, fB, ia2B)
    

C     Calculate the derivatives. 
      do i=1,ny
        dudx(i) = (- 2.d0*x**2*dudtau(i)*ia2B(i)**2 + 
     $       2.d0*x*dudtauB(i)*f(i)*ia2B(i)**2 - 
     $       2.d0*x*dudtau(i)*fB(i)*ia2B(i)**2 - 
     $       (x + fB(i))*ia2B(i)*
     $        (2.d0*lamb*x*ia2B(i) + 
     $       fB(i)*(2.d0*(d - 3.d0) - (d - 2.d0)*ia2B(i)))*u(i) - 
     $       6.d0*x*fB(i)*ia2(i)*uB(i) + 
     $       2.d0*d*x*fB(i)*ia2(i)*uB(i) - 
     $       6.d0*fB(i)**2*ia2(i)*uB(i) + 
     $       2.d0*d*fB(i)**2*ia2(i)*uB(i) + 
     $       6.d0*x*f(i)*ia2B(i)*uB(i) - 
     $       2.d0*d*x*f(i)*ia2B(i)*uB(i) - 
     $       2.d0*x*f(i)*ia2B(i)**2*uB(i) + 
     $       d*x*f(i)*ia2B(i)**2*uB(i) + 
     $       (d - 2.d0)*fB(i)*(x + fB(i))*ia2B(i)**2*v(i) - 
     $       2.d0*x*f(i)*ia2B(i)**2*vB(i) + 
     $       d*x*f(i)*ia2B(i)**2*vB(i))/
     $     (2.d0*x*(x + fB(i))**2*ia2B(i)**2)
      
       dvdx(i) = -(2.d0*x**2*dvdtau(i)*ia2B(i)**2 + 
     $        2.d0*x*dvdtauB(i)*f(i)*ia2B(i)**2 - 
     $        2.d0*x*dvdtau(i)*fB(i)*ia2B(i)**2 + 
     $        (d - 2.d0)*(x - fB(i))*fB(i)*ia2B(i)**2*u(i) - 
     $        2.d0*x*f(i)*ia2B(i)**2*uB(i) + 
     $        d*x*f(i)*ia2B(i)**2*uB(i) + 
     $        (x - fB(i))*ia2B(i)*
     $         (2.d0*lamb*x*ia2B(i) + 
     $       fB(i)*(6.d0 - 2.d0*d + (d - 2.d0)*ia2B(i)))*v(i) - 
     $        6.d0*x*fB(i)*ia2(i)*vB(i) + 
     $        2.d0*d*x*fB(i)*ia2(i)*vB(i) + 
     $        6.d0*fB(i)**2*ia2(i)*vB(i) - 
     $        2.d0*d*fB(i)**2*ia2(i)*vB(i) + 
     $        6.d0*x*f(i)*ia2B(i)*vB(i) - 
     $        2.d0*d*x*f(i)*ia2B(i)*vB(i) - 
     $        2.d0*x*f(i)*ia2B(i)**2*vB(i) + 
     $        d*x*f(i)*ia2B(i)**2*vB(i))/
     $      (2.d0*x*(x - fB(i))**2*ia2B(i)**2)
      
        djunkxidx(i) = 0.d0

        dfdx(i) = -((d - 3.d0)*(ia2(i)*fB(i) - f(i)*ia2B(i) 
     $               + f(i)*ia2B(i)**2))/(x*ia2B(i)**2)
     
      end do
      
C     Pack up dydx in the inverse way of unpacking y.
      call yfromfields(ny, dudx, dvdx, dfdx, dydx)
      
      return
      end 
