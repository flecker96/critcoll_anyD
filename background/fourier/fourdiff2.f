      subroutine fourdiff2(n, x, xp, y, yp, period)

      integer n, nymax
      include '../nymax.inc'
      double precision x(n), xp(n), y(n), yp(n), caux(2*nymax), period

      if (n .eq. 1) then
         xp(1) = 0.d0
         yp(1) = 0.d0
         return
      end if

      do i=1,n
         caux(2*i-1) = x(i) / n
         caux(2*i) = y(i) / n
      end do
      call four1(caux, n, 1)
      call mik(caux, caux, n, period)
      call four1(caux, n, -1)
      do i=1,n
         xp(i) = caux(2*i-1)
         yp(i) = caux(2*i)
      end do

      return
      end
