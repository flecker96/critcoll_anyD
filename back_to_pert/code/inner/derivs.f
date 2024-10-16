      subroutine derivs(ny, d, y, x, dydx)

      implicit none
      integer ny
      double precision y(ny), x, d, dydx(ny)

      integer nymax
      include '../nymax.inc'

      integer i
      double precision Delta, ia2(nymax),
     $     u(nymax), v(nymax), f(nymax), 
     $     dudtau(nymax), dvdtau(nymax), 
     $     dfdtau(nymax), 
     $     dudx(nymax), dvdx(nymax), dfdx(nymax)
      
C     Reconstruct fields in tau-space.
      call fieldsfromy(ny, d, y, x,
     $     u, v, f, ia2, Delta,
     $     dudtau, dvdtau, dfdtau)


C     Calculate the derivatives. 
      do i=1,ny
        dudx(i) = (f(i) * ((d - 2.d0) * v(i)
     $             - ( 2.d0 * (d - 3.d0) / ia2(i) - (d - 2.d0)) * u(i))
     $             - 2.d0 * x * dudtau(i)) / (2.d0 * x * (f(i) + x))
      
        dvdx(i) = (f(i) * ((d - 2.d0) * u(i)
     $             - ( 2.d0 * (d - 3.d0) / ia2(i) - (d - 2.d0)) * v(i))
     $             + 2.d0 * x * dvdtau(i)) / (2.d0 * x * (f(i) - x))

        dfdx(i) = (d - 3.d0) * f(i) * (1.d0 - ia2(i)) / (ia2 (i) * x)
      end do
      
      
C     Pack up dydx in the inverse way of unpacking y.
      call yfromfields(ny, dudx, dvdx, dfdx, dydx)

      dydx(5) = 0.d0
      
      return
      end 
