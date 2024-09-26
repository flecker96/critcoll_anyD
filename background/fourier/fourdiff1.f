      subroutine fourdiff1(n, x, xp, period)

      implicit none
      integer i, n, nymax
      include '../nymax.inc'
      double precision x(n), xp(n), caux(2*nymax), caux2(2*nymax),
     $   period

      do i=1,n
         caux(2*i-1) = x(i) / n
         caux(2*i) = 0.d0
      end do
      call four1(caux, n, 1)
      call mik(caux, caux2, n, period)
      call four1(caux2, n, -1)
      do i=1,n
         xp(i) = caux2(2*i-1)
      end do

      return
      end
