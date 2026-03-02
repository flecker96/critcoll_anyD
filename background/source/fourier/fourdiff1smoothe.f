      subroutine fourdiff1smoothe(n, x, xp, period)

      integer n, nymax
      include '../nymax.inc'
      double precision x(n), xp(n), caux(2*nymax), period

      do i=1,n
         caux(2*i-1) = x(i) / n
         caux(2*i) = 0.d0
      end do
      call four1(caux, n, 1)
      call halve(caux, caux, n/2)
      call halve(caux, caux, n/4)
      call double(caux, caux, n/4)
      call double(caux, caux, n/2)
      call mik(caux, caux, n, period)
      call four1(caux, n, -1)
      do i=1,n
         xp(i) = caux(2*i-1)
      end do

      return
      end
