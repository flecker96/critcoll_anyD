      subroutine compmismatch(ny, n3, x, yl, yr, vec)

      implicit none
      integer ny, n3
      double precision yl(ny), yr(ny), vec(n3), x,
     $     ul(ny/2), vl(ny/2), fl(ny/2), psil(ny/2),
     $     ur(ny/2), vr(ny/2), fr(ny/2), psir(ny/2)

      integer nymax
      include '../nymax.inc'

      integer i
     

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

      if (yl(5).ne.yr(5)) stop 'Error in Delta'

C     Replace the slot where Delta was stored with the correct value -Im(f_2) again
      yl(5) = - yl(ny-3)
      yr(5) = - yr(ny-3)
      

C     Fourier transform to real space with ny/2 points.
      call four1(yl, ny/2, -1)
      call four1(yr, ny/2, -1)

C     Separate real and imaginary, odd and even for 3 real functions 
C     on ny/2 points each. i labels points in real-space.
      do i=1,ny/4
         ul(i) = 0.5d0 * (yl(2*i-1) - yl(2*(i+ny/4)-1))
         vl(i) = 0.5d0 * (yl(2*i) - yl(2*(i+ny/4)))
         fl(i) = 0.5d0 * (yl(2*i) + yl(2*(i+ny/4)))
      end do

C     Complete by shifting symmetries.
      do i=1,ny/4
         ul(i+ny/4) = - ul(i)
         vl(i+ny/4) = - vl(i)
         fl(i+ny/4) = fl(i)
      end do

C     Do the same for the variables from the right
      do i=1,ny/4
         ur(i) = 0.5d0 * (yr(2*i-1) - yr(2*(i+ny/4)-1))
         vr(i) = 0.5d0 * (yr(2*i) - yr(2*(i+ny/4)))
         fr(i) = 0.5d0 * (yr(2*i) + yr(2*(i+ny/4)))
      end do

C     Complete by shifting symmetries.
      do i=1,ny/4
         ur(i+ny/4) = - ur(i)
         vr(i+ny/4) = - vr(i)
         fr(i+ny/4) = fr(i)
      end do

      do i=1, ny/2
         psil(i)  = ( vl(i) - ul(i) ) / 2.d0 / x**2
         psir(i)  = ( vr(i) - ur(i) ) / 2.d0 / x**2
      end do

      do i=1,n3/3
         vec(i       ) = fl(i) - fr(i)
         vec(i   +  n3/3) = ul(i) - ur(i)
         vec(i  + 2*n3/3) = psil(i) - psir(i)
      end do
     
      return
      end







