      subroutine fieldsfromylin(ny, d, y, x, 
     $     u, v, f, ia2, Delta, lamb,
     $     dudtau, dvdtau, dfdtau,
     $     uB, vB, fB, ia2B)

      implicit none
      integer ny
      double precision y(ny), x, d, lamb,
     $     u(ny), v(ny), f(ny), ia2(ny), Delta,
     $     dudtau(ny), dvdtau(ny), dfdtau(ny),
     $     uB(ny), vB(ny), fB(ny), ia2B(ny)
      integer nymax
      include '../nymax.inc'

      integer i
      double precision caux(2*nymax), caux1(2*nymax), 
     $     coeff1(nymax), coeff2(nymax)

C     y contains the odd Fourier components of delU+idelV and the even ones
C     of idelf. 
C     To reduce aliasing error, double the number of Fourier components,
C     padding the high frequencies in the middle with zeros. This means we
C     go from ny/2 points in frequency-space to ny.

      call double(y, caux, ny/2)   
      
C     Before Fourier transforming caux  in place, make a copy 
C     multiplied by -ik, to calculate derivatives later.
      call mik(caux, caux1, ny, Delta)

C     Fourier transform to real space with ny points.
      call four1(caux, ny, -1)

C     Separate real and imaginary, odd and even for 4 real functions 
C     on ny points each. i labels points in real-space.
      do i=1,ny/2
         u(i) = 0.5d0 * (caux(2*i-1) - caux(2*(i+ny/2)-1))
         v(i) = 0.5d0 * (caux(2*i) - caux(2*(i+ny/2)))
         f(i) = 0.5d0 * (caux(2*i) + caux(2*(i+ny/2)))
      end do

C     Complete by shifting symmetries.
      do i=1,ny/2
         u(i+ny/2) = - u(i)
         v(i+ny/2) = - v(i)
         f(i+ny/2) = f(i)
      end do

C     Now do the whole job again for the derivatives.
      call four1(caux1, ny, -1)
      do i=1,ny/2
         dudtau(i) = 0.5d0 * (caux1(2*i-1) - caux1(2*(i+ny/2)-1))
         dvdtau(i) = 0.5d0 * (caux1(2*i) - caux1(2*(i+ny/2)))
         dfdtau(i) = 0.5d0 * (caux1(2*i) + caux1(2*(i+ny/2)))
      end do
      do i=1,ny/2
         dudtau(i+ny/2) = - dudtau(i)
         dvdtau(i+ny/2) = - dvdtau(i)
         dfdtau(i+ny/2) = dfdtau(i)
      end do

C     Calculate ia2 from linearized constraint: ia2,tau + coeff1 ia2 + coeff2 = 0.
      do i=1,ny
        coeff1(i) = - ((- 2.d0 + d)**3.d0 * ((x + fB(i))*uB(i)**2 
     $                  + (x - fB(i))*vB(i)**2))/(8.d0 * x)
     $              + 3.d0 - d + lamb
        coeff2(i) = -((d - 2.d0)**3 * ia2B(i)
     $             *(2.d0 * u(i)*(x + fB(i))*uB(i) 
     $               + 2.d0 * v(i)*(x - fB(i))*vB(i) 
     $               + f(i)*(uB(i)**2 - vB(i)**2)))/(8.d0 * x)  
      end do

      call inhom(ia2, coeff1, coeff2, Delta, ny)
      
     
      return
      end







