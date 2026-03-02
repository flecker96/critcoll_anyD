      subroutine smoothe(n, x, xp, period)

      integer n, nmax
      parameter (nmax=1024)
      double precision x(n), xp(n), caux(2*nmax), period

      do i=1,n
         caux(2*i-1) = x(i) / n
         caux(2*i) = 0.d0
      end do
      call four1(caux, n, 1)
      call halve(caux, caux, n/2)
      call halve(caux, caux, n/4)
c      call halve(caux, caux, n/8)
c      call halve(caux, caux, n/16)
c      call double(caux, caux, n/16)
c      call double(caux, caux, n/8)
      call double(caux, caux, n/4)
      call double(caux, caux, n/2)
      call four1(caux, n, -1)
      do i=1,n
         xp(i) = caux(2*i-1)
      end do

      return
      end
