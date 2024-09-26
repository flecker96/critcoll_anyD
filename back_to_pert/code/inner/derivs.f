      subroutine derivs(ny, d, y, x, dydx,
     $     junkc, junkf1, junkf2)

      implicit none
      integer ny
      double precision y(ny), x, d, 
     $       dydx(ny), junkc, junkf1(ny), junkf2(ny)

      integer nymax
      include '../nymax.inc'

      integer i
      double precision Delta, ia2(nymax),
     $     u(nymax), v(nymax), junkxi(nymax), f(nymax), 
     $     dudtau(nymax), dvdtau(nymax), djunkxidtau(nymax), 
     $     dfdtau(nymax), 
     $     dudx(nymax), dvdx(nymax), djunkxidx(nymax), dfdx(nymax)
      
C     Reconstruct fields in tau-space.
      call fieldsfromy(ny, d, y, x,
     $     u, v, junkxi, f, ia2, Delta,
     $     dudtau, dvdtau, djunkxidtau, dfdtau)


C     Calculate the derivatives. 
      do i=1,ny
        dudx(i) = (f(i) * ((d - 2.d0) * v(i)
     $             - ( 2.d0 * (d - 3.d0) / ia2(i) - (d - 2.d0)) * u(i))
     $             - 2.d0 * x * dudtau(i)) / (2.d0 * x * (f(i) + x))
      
        dvdx(i) = (f(i) * ((d - 2.d0) * u(i)
     $             - ( 2.d0 * (d - 3.d0) / ia2(i) - (d - 2.d0)) * v(i))
     $             + 2.d0 * x * dvdtau(i)) / (2.d0 * x * (f(i) - x))

        djunkxidx(i) = 0.d0

        dfdx(i) = (d - 3.d0) * f(i) * (1.d0 - ia2(i)) / (ia2 (i) * x)
      end do
      
      
C     Pack up dydx in the inverse way of unpacking y.
      call yfromfields(ny, dudx, dvdx, dfdx, dydx)

      dydx(5) = 0.d0
      
      return
      end 
