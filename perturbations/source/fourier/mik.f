C     Multiplies Fourier array x by -ik, thus giving the derivative.
C     The high frequency is suppressed, as being a cosine.
C     x and kx can be the same arrray, as I use a temp variable for
C     the swap that corresponds to multiplying by i.

      subroutine mik(x, kx, n, period)
      implicit none
      integer n, i
      double precision x(2*n), kx(2*n), period, k1, temp

      k1 = 6.28318530717959d0 / period

C     Positive frequencies.
      do i=0,n/2-1
         temp = x(2*i+2) * i *k1
         kx(2*i+2) = - x(2*i+1) * i * k1
         kx(2*i+1) = temp
      end do

C     High frequency is suppressed, because cos -> sin -> 0.
      kx(n+1) = 0.d0
      kx(n+2) = 0.d0

C     Negative frequencies (excluding high frequency).
      do i=n/2+1,n-1
         temp = x(2*i+2) * (i-n) * k1
         kx(2*i+2) = - x(2*i+1) * (i-n) * k1
         kx(2*i+1) = temp
      end do

      return
      end
