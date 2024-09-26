C     Takes f(x) in n points x=0, period/n 2*period/n, ... (n-1)*period/n
C     and calculates f(x) at arbitrary x.

      subroutine fourinterpol(n, fn, period, x, fx)

      implicit none
      integer n
      double precision fn(n), period, x, fx

      integer nymax, i
      include '../nymax.inc'
      double precision fk(2*nymax), twopi, phi1, phi, refx, imfx, c, s
      parameter (twopi = 6.28318530717959d0)


      if (n .gt. nymax) stop 'nymax too small in fourinterpol'

C     Calculate complex Fourier coefficients
      do i=1,n
         fk(2*i-1) = fn(i) / n
         fk(2*i) = 0.d0
      end do
      call four1(fk, n, 1)

C     Calculate f(x) as sum over Fourier coefficients.
      refx = 0.d0
      imfx = 0.d0
      phi1 = twopi * x / period 
      do i=1,n
         if (i .le. n/2) then
            phi = (i-1) * phi1
         else
            phi = (i-1-n) * phi1
         end if
         c = cos(phi)
         s = sin(phi)
         refx = refx + c * fk(2*i-1) + s * fk(2*i)
         imfx = imfx - s * fk(2*i-1) + c * fk(2*i) 
      end do
c      write(6,*) imfx
      fx = refx

      return
      end 
