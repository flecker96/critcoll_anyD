      subroutine fieldsfromy(ny, d, y, x, 
     $     u, v, junkxi, f, ia2, Delta,
     $     dudtau, dvdtau, djunkxidtau, dfdtau)

      implicit none
      integer ny
      double precision y(ny), x, d,
     $     u(ny), v(ny), junkxi(ny), f(ny), ia2(ny), Delta,
     $     dudtau(ny), dvdtau(ny), djunkxidtau(ny), dfdtau(ny)

      integer nymax
      include '../nymax.inc'

      integer i
      double precision caux(2*nymax), caux1(2*nymax), 
     $     coeff1(nymax), coeff2(nymax)

C     y contains the odd Fourier components of U+iV and the even ones
C     of junkxi+if. By definition junkxi has no fundamental frequency
C     cosine. 
C
C     It is in fact identically zero in this new version of the inner patch,
C     but in old code for the inner patch we used it for xi0.
C
C     So in the Fourier transform of junkxi+if we set the real part of 
C     k=2, y(5), equal to minus the the real part of k=ny-2, y(2*ny-3). 
C     In y(5) we store Delta, which therefore must be saved first.
C     To reduce aliasing error, double the number of Fourier components,
C     padding the high frequencies in the middle with zeros. This means we
C     go from ny/2 points in frequency-space to ny.

      Delta = y(5)
      call double(y, caux, ny/2)
      caux(5) = - caux(2*ny-3)
      
      
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
         junkxi(i) = 0.5d0 * (caux(2*i-1) + caux(2*(i+ny/2)-1))
         f(i) = 0.5d0 * (caux(2*i) + caux(2*(i+ny/2)))
      end do

C     Complete by shifting symmetries.
      do i=1,ny/2
         u(i+ny/2) = - u(i)
         v(i+ny/2) = - v(i)
         junkxi(i+ny/2) = junkxi(i)
         f(i+ny/2) = f(i)
      end do

C     Now do the whole job again for the derivatives.
      call four1(caux1, ny, -1)
      do i=1,ny/2
         dudtau(i) = 0.5d0 * (caux1(2*i-1) - caux1(2*(i+ny/2)-1))
         dvdtau(i) = 0.5d0 * (caux1(2*i) - caux1(2*(i+ny/2)))
         djunkxidtau(i)= 0.5d0 * (caux1(2*i-1) + caux1(2*(i+ny/2)-1))
         dfdtau(i) = 0.5d0 * (caux1(2*i) + caux1(2*(i+ny/2)))
      end do
      do i=1,ny/2
         dudtau(i+ny/2) = - dudtau(i)
         dvdtau(i+ny/2) = - dvdtau(i)
         djunkxidtau(i+ny/2) = djunkxidtau(i)
         dfdtau(i+ny/2) = dfdtau(i)
      end do

C     Calculate ia2 from its constraint: ia2,tau + coeff1 ia2 + coeff2 = 0.

      do i=1,ny
        coeff1(i) = - ((x + f(i)) * u(i)**2 + (x - f(i)) * v(i)**2) 
     $                 * (d - 2.d0)**3 / (8.d0 * x) 
     $               - (d - 3.d0)
     
        coeff2(i) = (d - 3.d0)   
      end do

      call inhom(ia2, coeff1, coeff2, Delta, ny)
      
     
      return
      end







