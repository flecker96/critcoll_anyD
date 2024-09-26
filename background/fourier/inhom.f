C     Solves the equation x'+fx+g=0 where f and g are periodic
C     functions. The integration constant is fixed by demanding that 
C     x be also periodic. This version does it in closed form,
C     see notebook p.92.
C--------------------------------------------------------------------------
      subroutine inhom(x, f, g, period, n)
C--------------------------------------------------------------------------
      implicit none
      integer n
      double precision period, x(n), f(n), g(n)
C-------------------------------------------------------------------------
C     Uses lamint 
C--------------------------------------------------------------------------
      integer i, nymax
      include '../nymax.inc'
      double precision f0, bigf(nymax),
     $     h(nymax), h0, bigh(nymax)
C--------------------------------------------------------------------------

      call lamint(n, f, period, 0.d0, f0, bigf)

      do i=1,n
         h(i) = - g(i) * exp(bigf(i))
      end do

      call lamint(n ,h, period, f0, h0, bigh)

      do i=1,n
         x(i) = exp(- bigf(i)) * bigh(i)
      end do
      
      return
      end

